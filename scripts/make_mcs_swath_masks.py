import xarray as xr
import numpy as np
import pandas as pd
import cftime
import yaml
import calendar
import os, glob, re
import time
import argparse
import logging
import traceback
import sys
import gc
import intake
import easygems.healpix as egh
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.zarr_tools import setup_dask_client, initialize_zarr_store, append_chunk_to_zarr
from src.cof_paths import data_root
from src.utilities import convert_cftime_to_standard_calendar
from src.mcs_tc_filter import filter_mcs_tc_overlaps, MCS_TC_FILTER_THRESHOLD
from pyflextrkr.ft_utilities import load_config
from pyflextrkr.ftfunctions import olr_to_tb

# Import for parallel processing
try:
    from distributed import as_completed
except ImportError:
    as_completed = None  # Will only be needed if parallel=True

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Suppress verbose Dask logging
    logging.getLogger('distributed').setLevel(logging.WARNING)
    logging.getLogger('distributed.worker').setLevel(logging.WARNING)
    logging.getLogger('distributed.core').setLevel(logging.WARNING)
    logging.getLogger('distributed.comm').setLevel(logging.WARNING)
    logging.getLogger('distributed.nanny').setLevel(logging.WARNING)
    logging.getLogger('distributed.scheduler').setLevel(logging.WARNING)

#--------------------------------------------------------------------------------------------------
def create_track_swaths_and_coverage(tracknumber):
    """
    Create track swaths and coverage count arrays from a 3D track number array.
    
    Parameters:
    -----------
    tracknumber : numpy.ndarray
        3D array with dimensions (time, y, x) containing track numbers at each timestep.
        Background/no-track pixels should be 0.
    
    Returns:
    --------
    track_swaths_dict : dict
        Dictionary with track_id as keys and 2D swath arrays as values.
        Each swath shows all pixels covered by that track (labeled with track_id).
    track_coverage_dict : dict
        Dictionary with track_id as keys and 2D coverage count arrays as values.
        Each array shows the number of timesteps each pixel was covered by that track.
    
    Example:
    --------
    >>> tracknumber = np.zeros((24, 200, 150), dtype=int)
    >>> # ... populate tracknumber with track data ...
    >>> swaths, coverage = create_track_swaths_and_coverage(tracknumber)
    """
    # Get unique track numbers (excluding 0 which is background)
    unique_tracks = np.unique(tracknumber)
    unique_tracks = unique_tracks[unique_tracks > 0]
    
    # Initialize dictionaries
    track_swaths_dict = {}
    track_coverage_dict = {}
    
    # Process each unique track using vectorized operations
    for track_id in unique_tracks:
        # Create binary mask for this track across all timesteps
        track_mask = (tracknumber == track_id).astype(int)
        
        # Track swath: all pixels covered by this track labeled with track number
        swath = np.any(track_mask, axis=0).astype(int) * track_id
        track_swaths_dict[track_id] = swath
        
        # Coverage count: number of times each pixel is covered
        coverage = np.sum(track_mask, axis=0)
        track_coverage_dict[track_id] = coverage
    
    return track_swaths_dict, track_coverage_dict

#--------------------------------------------------------------------------------------------------
def combine_swaths_with_priority(track_swaths_dict, track_coverage_dict):
    """
    Combine multiple track swaths into a single 2D array, resolving overlaps
    by assigning the track number with the largest coverage count.
    
    Parameters:
    -----------
    track_swaths_dict : dict
        Dictionary with track_id as keys and 2D swath arrays as values
    track_coverage_dict : dict
        Dictionary with track_id as keys and 2D coverage count arrays as values
    
    Returns:
    --------
    combined_swath : numpy.ndarray
        2D array with combined swaths, where overlapping pixels are assigned
        to the track with the highest coverage count
    """
    # Get the shape from the first swath
    first_key = list(track_swaths_dict.keys())[0]
    shape = track_swaths_dict[first_key].shape
    
    # Initialize output array
    combined_swath = np.zeros(shape, dtype=int)
    
    # Initialize array to track maximum coverage at each pixel
    max_coverage = np.zeros(shape, dtype=int)
    
    # Process each track using vectorized operations
    for track_id in track_swaths_dict.keys():
        # Get swath and coverage for this track
        swath = track_swaths_dict[track_id]
        coverage = track_coverage_dict[track_id]
        
        # Create mask for pixels belonging to this track's swath
        track_mask = swath > 0
        
        # Update combined_swath where:
        # 1. Current pixel has no assignment yet (combined_swath == 0), OR
        # 2. Current track has higher coverage than previous assignment
        update_mask = track_mask & ((combined_swath == 0) | (coverage > max_coverage))
        
        # Apply updates using vectorized operations
        combined_swath = np.where(update_mask, track_id, combined_swath)
        max_coverage = np.where(update_mask, coverage, max_coverage)
    
    return combined_swath

#--------------------------------------------------------------------------------------------------
def create_latitude_dependent_tb_threshold(lat_values):
    """
    Create latitude-dependent brightness temperature thresholds.
    
    Parameters:
    -----------
    lat_values : numpy.ndarray
        Array of latitude values
        
    Returns:
    --------
    tb_thresh : numpy.ndarray
        Array of temperature thresholds (K) matching lat_values shape
    """
    abs_lat = np.abs(lat_values)
    tb_thresh = np.zeros_like(lat_values, dtype=np.float32)
    
    # Tropics: |lat| <= 30, tb_thresh = 250 K
    tropical_mask = abs_lat <= 30
    tb_thresh[tropical_mask] = 250.0
    
    # Mid-latitudes: 30 < |lat| <= 60, linearly decrease from 250 to 230 K
    midlat_mask = (abs_lat > 30) & (abs_lat <= 60)
    tb_thresh[midlat_mask] = 250.0 - 20.0 * ((abs_lat[midlat_mask] - 30.0) / 30.0)
    
    # High latitudes: |lat| > 60, tb_thresh = 230 K
    highlat_mask = abs_lat > 60
    tb_thresh[highlat_mask] = 230.0
    
    return tb_thresh

#--------------------------------------------------------------------------------------------------
def classify_cloud_types(tb, pr, tb_thresh, mcs_mask, pr_threshold=0.5):
    """
    Classify cloud types based on brightness temperature and precipitation.
    
    Parameters:
    -----------
    tb : numpy.ndarray
        Brightness temperature (K)
    pr : numpy.ndarray
        Precipitation rate (mm/h)
    tb_thresh : numpy.ndarray
        Latitude-dependent brightness temperature threshold (K)
    mcs_mask : numpy.ndarray
        MCS mask array (values > 0 indicate MCS pixels)
        Classification only occurs where mcs_mask == 0 (non-MCS areas)
    pr_threshold : float, optional
        Precipitation threshold for classification (default: 0.5 mm/h)
    
    Returns:
    --------
    cloud_type : numpy.ndarray
        Cloud type classification:
        0 = Unclassified (inside MCS regions)
        1 = Deep convective (tb < tb_thresh & pr >= pr_threshold)
        2 = Stratiform (tb < tb_thresh & pr < pr_threshold)
        3 = Non-deep convective (tb >= tb_thresh & pr >= pr_threshold)
        4 = Drizzle (tb >= tb_thresh & pr < pr_threshold)
    """
    # Initialize cloud type array with zeros
    cloud_type = np.zeros_like(tb, dtype=np.int8)
    
    # Only classify where mcs_mask == 0 (non-MCS areas)
    valid_mask = mcs_mask == 0
    
    # Classification conditions (only applied where valid_mask is True)
    # 1. Deep convective: tb < tb_thresh & pr >= pr_threshold
    deep_conv = (tb < tb_thresh) & (pr >= pr_threshold) & valid_mask
    cloud_type[deep_conv] = 1
    
    # 2. Stratiform: tb < tb_thresh & pr < pr_threshold
    stratiform = (tb < tb_thresh) & (pr < pr_threshold) & valid_mask
    cloud_type[stratiform] = 2
    
    # 3. Non-deep convective: tb >= tb_thresh & pr >= pr_threshold
    nondeep_conv = (tb >= tb_thresh) & (pr >= pr_threshold) & valid_mask
    cloud_type[nondeep_conv] = 3
    
    # 4. Drizzle: tb >= tb_thresh & pr < pr_threshold
    drizzle = (tb >= tb_thresh) & (pr < pr_threshold) & valid_mask
    cloud_type[drizzle] = 4
    
    return cloud_type

#--------------------------------------------------------------------------------------------------
def find_priority_based_cloud_type(cloud_types_time_series):
    """
    Find cloud type based on priority ranking: 1 > 2 > 3 > 4 (VECTORIZED).
    
    If a cell has any occurrence of type 1, assign type 1.
    If no type 1 but has type 2, assign type 2.
    If no type 1 or 2 but has type 3, assign type 3.
    If no type 1, 2, or 3 but has type 4, assign type 4.
    Otherwise, assign 0 (unclassified).
    
    Parameters:
    -----------
    cloud_types_time_series : numpy.ndarray
        Array of shape (n_times, n_cells) with cloud type values starting from 0
        
    Returns:
    --------
    priority_based : numpy.ndarray
        Array of shape (n_cells,) with priority-based cloud type
    """
    n_times, n_cells = cloud_types_time_series.shape
    
    # Initialize result array with zeros
    priority_based = np.zeros(n_cells, dtype=np.int8)
    
    # Check for presence of each cloud type (any occurrence over time)
    has_type_1 = np.any(cloud_types_time_series == 1, axis=0)
    has_type_2 = np.any(cloud_types_time_series == 2, axis=0)
    has_type_3 = np.any(cloud_types_time_series == 3, axis=0)
    has_type_4 = np.any(cloud_types_time_series == 4, axis=0)
    
    # Apply priority: 1 > 2 > 3 > 4
    # Start from lowest priority and work up (so higher priorities overwrite)
    priority_based[has_type_4] = 4
    priority_based[has_type_3] = 3
    priority_based[has_type_2] = 2
    priority_based[has_type_1] = 1
    
    return priority_based

#--------------------------------------------------------------------------------------------------
def compute_mean_precip_by_cloud_type(cloud_types_timeseries, pr_timeseries):
    """
    Compute frequency-weighted mean precipitation for each cloud type over the aggregation window.
    
    This computes: (conditional mean) × (frequency of that type)
    When summed across all 4 types, this equals the simple mean of all precipitation.

    Parameters
    ----------
    cloud_types_timeseries : np.ndarray
        Array of shape (n_times, n_cells) with cloud type values (0-4).
    pr_timeseries : np.ndarray
        Array of shape (n_times, n_cells) with precipitation values (mm/h).

    Returns
    -------
    dc_pr : np.ndarray
        Frequency-weighted deep convective precipitation (cloud_type==1), shape (n_cells,)
    st_pr : np.ndarray
        Frequency-weighted stratiform precipitation (cloud_type==2), shape (n_cells,)
    nd_pr : np.ndarray
        Frequency-weighted non-deep convective precipitation (cloud_type==3), shape (n_cells,)
    dz_pr : np.ndarray
        Frequency-weighted drizzle precipitation (cloud_type==4), shape (n_cells,)
        
    Note
    ----
    dc_pr + st_pr + nd_pr + dz_pr = mean(pr_timeseries, axis=0)
    """
    import warnings
    
    n_times = cloud_types_timeseries.shape[0]
    
    # Mask precipitation by cloud type for each time
    dc_mask = cloud_types_timeseries == 1
    st_mask = cloud_types_timeseries == 2
    nd_mask = cloud_types_timeseries == 3
    dz_mask = cloud_types_timeseries == 4

    # Compute conditional means (average precipitation WHEN cell is this type)
    # Suppress RuntimeWarnings for empty slices
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', r'Mean of empty slice')
        dc_pr_cond = np.nanmean(np.where(dc_mask, pr_timeseries, np.nan), axis=0)
        st_pr_cond = np.nanmean(np.where(st_mask, pr_timeseries, np.nan), axis=0)
        nd_pr_cond = np.nanmean(np.where(nd_mask, pr_timeseries, np.nan), axis=0)
        dz_pr_cond = np.nanmean(np.where(dz_mask, pr_timeseries, np.nan), axis=0)
    
    # Compute frequency of each type (fraction of time steps)
    dc_freq = dc_mask.sum(axis=0) / n_times
    st_freq = st_mask.sum(axis=0) / n_times
    nd_freq = nd_mask.sum(axis=0) / n_times
    dz_freq = dz_mask.sum(axis=0) / n_times
    
    # Weight conditional means by frequency
    # Set NaN conditional means to 0 before weighting (cells that never had this type)
    dc_pr = np.nan_to_num(dc_pr_cond, nan=0.0) * dc_freq
    st_pr = np.nan_to_num(st_pr_cond, nan=0.0) * st_freq
    nd_pr = np.nan_to_num(nd_pr_cond, nan=0.0) * nd_freq
    dz_pr = np.nan_to_num(dz_pr_cond, nan=0.0) * dz_freq

    return dc_pr, st_pr, nd_pr, dz_pr

#--------------------------------------------------------------------------------------------------
def remove_pr_where_tb_missing(pr, tb):
    """
    Set the precipitation to NaN at every pixel and hour where Tb is missing.

    Without Tb a pixel cannot be classified (there is no cloud type), so its rain would be counted in tot_pr but
    in no category and would show up as a residual. Removing it here, before anything is derived from pr,
    makes tot_pr and dc/st/nd/dz_pr see the same NaN (a missing hour counts as zero with the same divisor),
    so the budget closes. Where Tb is valid nothing changes. The sources whose Tb comes from OLR are missing Tb
    only where pr is missing too, so this changes nothing for them; the IR-IMERG Tb has gaps (about 0.3% of
    the cell-hours within 60S-60N, and about 0.3% of the rain fell there).

    Parameters:
    -----------
    pr : xarray.DataArray
        Precipitation (time, cell), mm/h
    tb : xarray.DataArray
        Brightness temperature on the same (time, cell) grid

    Returns:
    --------
    (pr_masked, stats) : the masked precipitation and a dict with
        'px_hours'     : pixel-hours that had rain (> 0) and no Tb
        'rain_removed' : that rain, summed over pixels and hours (mm/h; HEALPix cells have equal area)
        'rain_total'   : all rain before the removal, summed the same way
    """
    no_tb = tb.isnull()
    stats = {
        'px_hours': int(((pr > 0) & no_tb).sum()),
        'rain_removed': float(pr.where(no_tb).sum()),
        'rain_total': float(pr.sum()),
    }
    return pr.where(~no_tb), stats

#--------------------------------------------------------------------------------------------------
def add_tb_pr_to_dataset(_ds, config):
    """
    Add brightness temperature (tb) and precipitation (pr) variables to a dataset chunk.
    
    This function reads OLR and precipitation data from a catalog, handles calendar conversions
    if necessary, and adds tb and pr variables to the input dataset.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with time coordinate
    config : dict
        Configuration dictionary containing catalog information and variable names:
        - catalog_file: Path to the intake catalog file
        - catalog_location: Location within the catalog
        - catalog_source: Source name within the catalog
        - catalog_params: Parameters for the catalog source
        - varname_precip_liq: Variable name of the precipitation field used as 'pr'
        - varname_precip_ice: Optional variable name of a separate frozen-precipitation field to add to
          'pr'. main() passes None: frozen precipitation is deliberately ignored so every source uses its
          own 'pr' field only.
        - pcp_convert_factor: Conversion factor for precipitation (e.g., to convert to mm/h)
    
    Returns:
    --------
    _ds : xarray.Dataset
        Dataset with added 'tb' and 'pr' variables. 'pr' is NaN wherever 'tb' is missing (see
        remove_pr_where_tb_missing); what that removed is in _ds.attrs['tb_gap'].
    """
    # Extract config values
    catalog_file = config['catalog_file']
    catalog_location = config['catalog_location']
    catalog_source = config['catalog_source']
    catalog_params = config['catalog_params']
    varname_olr = config['varname_olr']
    varname_precip_liq = config['varname_precip_liq']
    varname_precip_ice = config['varname_precip_ice']
    pr_convert_factor = config['pcp_convert_factor']
    
    # Special handling for IR_IMERG: tb is already in the dataset, no OLR conversion needed
    if catalog_source == 'IR_IMERG':
        varname_olr = 'Tb'
        varname_precip_liq = 'precipitation'

    # Read OLR/precipitation data from catalog
    if catalog_location is None:
        cat = intake.open_catalog(catalog_file)
    else:
        cat = intake.open_catalog(catalog_file)[catalog_location]
    ds_p = cat[catalog_source](**catalog_params).to_dask()

    # List of variables to keep
    vars_keep = [varname_olr, varname_precip_liq, varname_precip_ice]
    # Identify variables to drop (those not in the keep list)
    dropvars_list = [v for v in ds_p.data_vars if v not in vars_keep]
    # Drop unwanted variables
    ds_p = ds_p.drop_vars(dropvars_list, errors='ignore')

    # 1. Determine calendar types and convert if needed
    ds_p_calendar_type = type(ds_p.time.values[0]).__name__
    ds_calendar_type = type(_ds.time.values[0]).__name__

    # Convert ds_p time to match _ds if they differ (convert non-standard to standard)
    if ds_p_calendar_type != ds_calendar_type:
        # Convert cftime objects to numpy datetime64 (standard calendar)
        converted_times = convert_cftime_to_standard_calendar(ds_p.time.values)
        
        # Create new dataset with converted time coordinate
        ds_p = ds_p.assign_coords(time=converted_times)

    # 2. Find common time range between datasets
    common_times = sorted(set(ds_p['time'].values)
                            .intersection(set(_ds['time'].values)))
    
    if not common_times:
        raise ValueError("No common time values between mask dataset and catalog dataset!")
    
    # 3. Select ds_p to only the common times (BEFORE expensive operations).
    # Eagerly load the small time slice into memory here so the full multi-year
    # catalog task graph (with its native 262144-cell chunks) is released before
    # we derive pr/tb. Without .load(), the 4x chunk size mismatch between the
    # catalog (cell=262144) and the MCS mask zarr (cell=65536) causes a large
    # intermediate memory spike when _ds.load() is called later in the worker.
    ds_p = ds_p.sel(time=common_times).load()
    
    # 4. Convert precipitation units (only for common times)
    # Check liquid precipitation variable
    if varname_precip_liq in list(ds_p.keys()):
        # Convert liquid precipitation to mm/h (if conversion factor is provided)
        if pr_convert_factor is not None:
            pr = ds_p[varname_precip_liq] * pr_convert_factor
        else:
            pr = ds_p[varname_precip_liq]
    # Check if the ice precipitation variable exist in the dataset
    if varname_precip_ice in list(ds_p.keys()):
        # Convert ice precipitation to liquid equivalent (if conversion factor is provided)
        if pr_convert_factor is not None:
            prs = ds_p[varname_precip_ice] * pr_convert_factor
        else:
            prs = ds_p[varname_precip_ice]
        # Add ice precipitation to get total precipitation
        pr = pr + prs
    
    # Convert OLR to Tb (only for common times)
    # For IR_IMERG, tb is already in the dataset, skip conversion
    if catalog_source == 'IR_IMERG':
        tb = ds_p[varname_olr]  # varname_olr is already set to 'tb' for IR_IMERG
    else:
        tb = olr_to_tb(ds_p[varname_olr])
    
    # Rain at pixels whose Tb is missing cannot be classified, so it is removed before anything is derived from pr
    pr, tb_gap = remove_pr_where_tb_missing(pr, tb)

    # 5. Add precipitation & tb to the dataset
    _ds = _ds.sel(time=common_times)
    _ds["pr"] = pr
    _ds["tb"] = tb
    _ds.attrs['tb_gap'] = tb_gap

    return _ds

#--------------------------------------------------------------------------------------------------
def add_tc_mask_to_dataset(_ds, config):
    """
    Add a TC track mask ('tc_mask') variable to a dataset chunk, already aligned to
    _ds's native time coordinate and HEALPix cell grid.

    This is needed so process_timechunk_swath() can exclude TC-contaminated MCS tracks
    (via the shared filter_mcs_tc_overlaps(), see src/mcs_tc_filter.py) *before* running
    cloud-type classification, instead of Step 3 (make_cooccurrence_masks.py) discovering
    the same overlap after this script has already zeroed cloud_types/dc_pr/st_pr/nd_pr/
    dz_pr for the (soon to be partly invalid) MCS swath. See src/mcs_tc_filter.py's
    docstring for the full explanation of why this matters.

    Two loading modes, selected by which config keys are present (checked in this order):
      - config['tc_source_zarr']: a single zarr store already on the target grid with a
        'tc_mask' or 'TC_int_tag' variable (e.g. the ERA5-derived tracking used for the
        observational IMERGv7 source, which has no per-source TempestExtremes output of
        its own - see scripts/remap_era5_masks_healpix.py).
      - config['dir_te'] + config['source_te'] + config['source_res']: one or more
        per-source TempestExtremes NetCDF track files, following the same
        TC_test_tracks_{source_te}_{source_res}.*.nc convention combine_tracking_masks.py
        (Step 2) already uses. config['tc_glob_pattern'] may override the filename
        template entirely (glued onto dir_te) for sources that don't follow it (e.g. the
        bare "TC_test_tracks_era5_*.nc" pattern used for the ERA5 tracking files).

    Missing TC data for a given native timestep defaults tc_mask to 0 ("no TC") rather
    than dropping the timestep - unlike Tb/pr, a missing TC file must not shrink this
    chunk, and "assume no TC" is the conservative, already-status-quo behavior for any
    time/source this function cannot resolve.

    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with a time coordinate and 'cell' dimension (same grid as
        mcs_mask/tb/pr).
    config : dict
        See loading modes above.

    Returns:
    --------
    _ds : xarray.Dataset
        Dataset with an added 'tc_mask' variable, shape (time, cell) matching _ds.
    """
    tc_source_zarr = config.get('tc_source_zarr')
    dir_te = config.get('dir_te')
    source_te = config.get('source_te')
    source_res = config.get('source_res')
    tc_glob_pattern = config.get('tc_glob_pattern')

    if tc_source_zarr:
        ds_tc = xr.open_dataset(tc_source_zarr, engine='zarr', chunks={}, mask_and_scale=False)
    elif dir_te and (tc_glob_pattern or (source_te and source_res)):
        pattern = tc_glob_pattern if tc_glob_pattern else f"TC_test_tracks_{source_te}_{source_res}.*.nc"
        files_tc = sorted(glob.glob(os.path.join(dir_te, pattern)))
        if not files_tc:
            raise FileNotFoundError(f"No TC track files found: {os.path.join(dir_te, pattern)}")
        datasets = []
        for f in files_tc:
            try:
                datasets.append(xr.open_dataset(f, chunks={}, mask_and_scale=False))
            except Exception as e:
                print(f"  Warning: could not open TC file {f}: {e}")
        if not datasets:
            raise ValueError("Could not open any TC track files")
        ds_tc = xr.concat(datasets, dim='time').sortby('time')
        _, index = np.unique(ds_tc['time'].values, return_index=True)
        if len(index) < len(ds_tc['time']):
            ds_tc = ds_tc.isel(time=sorted(index))
    else:
        raise ValueError(
            "add_tc_mask_to_dataset requires either config['tc_source_zarr'], or "
            "config['dir_te'] plus config['tc_glob_pattern'] (or config['source_te'] and "
            "config['source_res']) to locate this source's TC track data."
        )

    # Standardize dims/variable name to match the MCS mask grid
    if 'ncol' in ds_tc.dims:
        ds_tc = ds_tc.rename({'ncol': 'cell'})
    if 'TC_int_tag' in ds_tc.variables:
        ds_tc = ds_tc.rename({'TC_int_tag': 'tc_mask'})
    for v in ['lat', 'lon']:
        if v in ds_tc.variables:
            ds_tc = ds_tc.drop_vars(v)
    if 'tc_mask' not in ds_tc.variables:
        raise KeyError(
            f"TC dataset has neither 'tc_mask' nor 'TC_int_tag' variable; found: "
            f"{list(ds_tc.data_vars)}"
        )

    # Calendar alignment, same approach as add_tb_pr_to_dataset above
    tc_calendar_type = type(ds_tc.time.values[0]).__name__
    ds_calendar_type = type(_ds.time.values[0]).__name__
    if tc_calendar_type != ds_calendar_type:
        converted_times = convert_cftime_to_standard_calendar(ds_tc.time.values)
        ds_tc = ds_tc.assign_coords(time=converted_times)

    # Reindex (not intersect) onto _ds's native time steps. TC track data is only
    # available 6-hourly for every source (confirmed: all 5 GSRM TempestExtremes outputs
    # and the ERA5-derived IMERGv7 tracks), coarser than the hourly-or-finer native
    # cadence _ds is at here. A plain reindex with a zero fill_value would leave every
    # native hour that doesn't exactly match a 6-hourly label at "no TC" by default,
    # making the whole MCS-TC exclusion nearly a no-op. Forward-fill instead: each native
    # hour picks up the most recent (<=) 6-hourly TC snapshot, i.e. "TC was present in
    # this cell for the 6-h window this hour falls in" - matching the granularity the TC
    # data actually has. tolerance caps the fill at just under one 6-h step so a genuine
    # gap in TC data (a missing window) still falls through to NaN -> 0 ("no TC") below,
    # rather than incorrectly stretching a stale value across a real gap.
    tc_mask = ds_tc['tc_mask'].reindex(
        time=_ds['time'].values, method='ffill', tolerance=pd.Timedelta(hours=5, minutes=59)
    ).load()
    _ds = _ds.assign(tc_mask=tc_mask)

    return _ds

#--------------------------------------------------------------------------------------------------
def process_timechunk_swath(_ds, tb_thresh=None, verbose=False):
    """
    Process a time chunk of dataset to create MCS swath masks and cloud type classification.
    
    Parameters:
    -----------
    _ds : xarray.Dataset
        Input dataset chunk with dimensions (time, cell) containing 'mcs_mask', 'ccs_mask', 
        'tb', 'pr', and coordinates 'lat', 'lon'.
    tb_thresh : numpy.ndarray, optional
        Pre-computed latitude-dependent brightness temperature threshold.
        If None, will be computed from lat coordinate.
    verbose : bool
        If True, print progress information.
    
    Returns:
    --------
    dict : Dictionary containing:
        'mcs_mask': MCS swath mask (1D array, cell dimension)
        'cloud_types': Aggregated cloud type classification (1D array, cell dimension)
        'dc_pr', 'st_pr', 'nd_pr', 'dz_pr': Frequency-weighted mean precipitation by cloud type (mm/h)
        'tot_pr': Window-mean total precipitation (mm/h), from the same hourly pr, not masked by the swath
        'input_all_nan': True when the window has no valid pr and no valid Tb (not written to the store)
    """
    if verbose:
        print(f"Processing time chunk with {len(_ds.time)} time steps...")
    
    # Extract track number arrays
    mcs_mask = _ds['mcs_mask'].values  # shape (time, cell)
    # Replace NaN with 0 (when mask_and_scale=True, fill_value becomes NaN)
    mcs_mask = np.nan_to_num(mcs_mask, nan=0.0).astype(int)

    # ===== MCS-TC exclusion (must happen BEFORE cloud-type classification below) =====
    # Exclude MCS tracks that significantly overlap a TC in this window, using the same
    # test (src/mcs_tc_filter.filter_mcs_tc_overlaps, MCS_TC_FILTER_THRESHOLD) that Step 3
    # (make_cooccurrence_masks.py) applies downstream. Doing this here, before
    # classify_cloud_types() runs on each native timestep, is what lets pixels freed by
    # the exclusion get real DC/ND/ST/DZ classification from Tb/pr instead of being stuck
    # at the zero values Step 3's later, isolated exclusion could not undo. See
    # src/mcs_tc_filter.py's module docstring for the full explanation.
    # 'tc_mask' is only present when the caller has provided TC track data (see
    # add_tc_mask_to_dataset); skip gracefully (old behavior) if it hasn't, e.g. small
    # manual/test invocations of this function.
    if 'tc_mask' in _ds.variables:
        tc_mask_arr = np.nan_to_num(_ds['tc_mask'].values, nan=0.0).astype(int)
        mcs_tc_filter_result = filter_mcs_tc_overlaps(
            mcs_mask=xr.DataArray(mcs_mask, dims=['time', 'cell']),
            tc_mask=xr.DataArray(tc_mask_arr, dims=['time', 'cell']),
            overlap_threshold=MCS_TC_FILTER_THRESHOLD,
            verbose=verbose,
        )
        mcs_mask = mcs_tc_filter_result['mcs_filtered'].values.astype(int)
        # Write the TC-filtered mask back so every later read of _ds['mcs_mask'] in this
        # function (the per-native-hour classify_cloud_types loop below) sees the same,
        # already-excluded mask rather than the pre-filter one.
        _ds = _ds.assign(mcs_mask=(_ds['mcs_mask'].dims, mcs_mask))
    elif verbose:
        print("  No 'tc_mask' found on this chunk; skipping MCS-TC exclusion "
              "(cloud_types/dc_pr/etc. may undercount precip from MCS tracks that "
              "overlap a TC - see src/mcs_tc_filter.py).")

    # # Sum CCS mask over time and convert to binary
    # # First replace NaN with 0, then sum
    # ccs_mask_values = _ds['ccs_mask'].values
    # ccs_mask_values = np.nan_to_num(ccs_mask_values, nan=0.0)
    # ccs_mask_sum = ((ccs_mask_values > 0).sum(axis=0) > 0).astype(int)  # shape (cell)

    # Create swaths and coverage for MCS
    mcs_swaths_dict, mcs_coverage_dict = create_track_swaths_and_coverage(mcs_mask)
    if not mcs_swaths_dict:
        # No MCS present anywhere in this aggregation window (e.g., MCS-free window or short
        # test subset). Fall back to an all-zero swath instead of indexing into an empty dict.
        combined_mcs_swath = np.zeros(mcs_mask.shape[1:], dtype=int)
    else:
        combined_mcs_swath = combine_swaths_with_priority(mcs_swaths_dict, mcs_coverage_dict)

    # Filter out CCS that overlap with MCS swaths
    # ccs_mask_sum[combined_mcs_swath > 0] = 0
    
    # ===== Cloud Type Classification =====
    # Create latitude-dependent tb threshold if not provided
    if tb_thresh is None:
        lat_values = _ds['lat'].values
        tb_thresh = create_latitude_dependent_tb_threshold(lat_values)
    
    # A window with no valid pr and no valid Tb at all (e.g. the first steps of a model run that are
    # missing because of spin-up, depending on how the data were post-processed) has nothing to classify.
    # Flag it so the caller can tell it apart from a failed task: all-NaN input is expected, not an error.
    input_all_nan = not (np.isfinite(_ds['tb'].values).any() or np.isfinite(_ds['pr'].values).any())

    # Classify cloud types for each time step
    cloud_types_timeseries = []
    for t in range(len(_ds.time)):
        tb_t = _ds['tb'].isel(time=t).values
        pr_t = _ds['pr'].isel(time=t).values
        mcs_mask_t = _ds['mcs_mask'].isel(time=t).values
        
        # Replace NaN with 0 (when mask_and_scale=True, fill_value becomes NaN)
        mcs_mask_t = np.nan_to_num(mcs_mask_t, nan=0.0).astype(int)
        
        # Classify cloud types (only in non-MCS areas)
        cloud_type_t = classify_cloud_types(tb_t, pr_t, tb_thresh, mcs_mask_t, pr_threshold=0.5)
        cloud_types_timeseries.append(cloud_type_t)
    
    # Stack into array (time, cell)
    cloud_types_timeseries = np.stack(cloud_types_timeseries, axis=0)
    
    # Extract precipitation time series for computing precipitation by cloud type
    pr_timeseries = _ds['pr'].values  # shape (time, cell)
    
    # Compute frequency-weighted mean precipitation for each cloud type
    dc_pr, st_pr, nd_pr, dz_pr = compute_mean_precip_by_cloud_type(cloud_types_timeseries, pr_timeseries)

    # Window-mean total precipitation from the same hourly pr as the classification above (not masked by
    # the swath). Missing steps count as zero and the denominator is n_steps, exactly as for dc_pr..dz_pr,
    # so those four sum to tot_pr outside the swath and any unclassified rain stays visible as residual.
    tot_pr = np.nan_to_num(pr_timeseries, nan=0.0).sum(axis=0) / pr_timeseries.shape[0]

    # Find cloud type based on priority ranking (1 > 2 > 3 > 4)
    cloud_types_aggregated = find_priority_based_cloud_type(cloud_types_timeseries)
    
    if verbose:
        n_classified_before = (cloud_types_aggregated > 0).sum()
        print(f"  Cloud types before MCS filter: {n_classified_before:,} cells")
    
    # Apply MCS priority: Set cloud type to 0 where MCS swath exists
    # Classification hierarchy: MCS > Cloud Types > Unclassified
    # If a cell has an MCS at any time during the window, MCS takes priority
    cloud_types_aggregated = np.where(combined_mcs_swath > 0, 0, cloud_types_aggregated)
    
    # Filter precipitation by cloud type: set to 0 where MCS swath exists
    # This makes the 4 cloud types mutually exclusive with MCS precipitation
    mcs_mask = combined_mcs_swath > 0
    dc_pr = np.where(mcs_mask, 0, dc_pr)
    st_pr = np.where(mcs_mask, 0, st_pr)
    nd_pr = np.where(mcs_mask, 0, nd_pr)
    dz_pr = np.where(mcs_mask, 0, dz_pr)

    # import matplotlib.pyplot as plt
    # cloud_types_da = xr.DataArray(cloud_types_aggregated, dims=['cell'], coords={'lat': ('cell', _ds['lat'].values), 'lon': ('cell', _ds['lon'].values), 'cell': _ds['cell'].values})
    # egh.healpix_show(cloud_types_da.where(cloud_types_da > 0), cmap='tab10', vmin=0, vmax=3)
    # import pdb; pdb.set_trace()
    if verbose:
        n_classified = (cloud_types_aggregated > 0).sum()
        n_mcs_overlap = n_classified_before - n_classified
        print(f"  Cloud types after MCS filter: {n_classified:,} cells")
        print(f"  Removed {n_mcs_overlap:,} cells due to MCS overlap")
        print(f"  ✅ Completed processing for this time chunk")

    return {
        'mcs_mask': combined_mcs_swath,
        'cloud_types': cloud_types_aggregated,
        'dc_pr': dc_pr,
        'st_pr': st_pr,
        'nd_pr': nd_pr,
        'dz_pr': dz_pr,
        'tot_pr': tot_pr,
        # Not a mask variable (the zarr writer only takes the names in mask_variables): tells the caller that
        # the input window had no valid pr or Tb, so an all-zero result is expected rather than a failure.
        'input_all_nan': input_all_nan,
    }

def process_timechunk_wrapper_zarr(start_idx, end_idx, zarr_path, config, verbose=False,
                                    store_offsets=None):
    """
    Wrapper function for processing a time chunk by reading from zarr file.
    
    This approach avoids serialization issues by:
    1. Having each worker read the zarr file directly
    2. Using integer indices instead of passing time arrays (reduces graph size)
    
    Args:
        start_idx: Start index in the time dimension
        end_idx: End index in the time dimension (exclusive)
        zarr_path: Path (str) or list of paths to the zarr store(s).
            When a list is provided the stores are opened and concatenated
            along the time dimension before slicing.
        config: Configuration dictionary with catalog and variable information
        verbose: Whether to print verbose output
        store_offsets: Optional list of (global_start, global_end_exclusive, path) tuples.
            When provided, each worker opens only the store(s) that contain the
            requested time range instead of opening and concatenating all stores.
        
    Returns:
        tuple: (time_str, results_dict) on success, or (None, {'error': 'ExcType: message'}) if the
        chunk raised, so the caller can report why (stream_process_to_zarr / resolve_chunk_result)
    """
    try:
        # Open only the zarr store(s) needed for this time slice.
        # store_offsets lets each worker target a single annual store rather than
        # opening and concatenating all stores (avoids per-worker memory overhead).
        open_stores = []
        if store_offsets is not None:
            relevant = [(gs, ge, p) for gs, ge, p in store_offsets
                        if gs < end_idx and ge > start_idx]
            if len(relevant) == 1:
                gs, _ge, path = relevant[0]
                _s = xr.open_dataset(path, engine='zarr', chunks=None)
                open_stores.append(_s)
                _ds = _s.isel(time=slice(start_idx - gs, end_idx - gs))
            else:
                parts = []
                for gs, ge, path in relevant:
                    _s = xr.open_dataset(path, engine='zarr', chunks=None)
                    open_stores.append(_s)
                    local_start = max(0, start_idx - gs)
                    local_end = min(ge - gs, end_idx - gs)
                    parts.append(_s.isel(time=slice(local_start, local_end)))
                _ds = xr.concat(parts, dim='time')
        elif isinstance(zarr_path, list):
            for p in zarr_path:
                _s = xr.open_dataset(p, engine='zarr', chunks=None)
                open_stores.append(_s)
            _ds = xr.concat(open_stores, dim='time').isel(time=slice(start_idx, end_idx))
        else:
            _s = xr.open_dataset(zarr_path, engine='zarr', chunks=None)
            open_stores.append(_s)
            _ds = _s.isel(time=slice(start_idx, end_idx))
        
        # Get the first time value for output
        out_time_val = _ds.time.values[0]
        
        # Load zarr data into memory and close the open stores to free memory
        _ds = _ds.load()
        _ds = _ds.pipe(egh.attach_coords)
        for _s in open_stores:
            _s.close()
        
        # Add tb and pr variables to the dataset chunk (pr is removed where Tb is missing)
        _ds = add_tb_pr_to_dataset(_ds, config)
        tb_gap = _ds.attrs.pop('tb_gap', None)
        # Add the TC track mask so process_timechunk_swath can exclude TC-contaminated
        # MCS tracks before cloud-type classification (see add_tc_mask_to_dataset's
        # docstring and src/mcs_tc_filter.py). Config controls this per-source; a source
        # with no TC config keys set simply skips exclusion (logged inside
        # process_timechunk_swath), matching pre-fix behavior for that source only.
        if config.get('tc_source_zarr') or config.get('dir_te'):
            _ds = add_tc_mask_to_dataset(_ds, config)
        # Eagerly load tb/pr into numpy to release the catalog's dask task graphs.
        # Without this, each .values call in process_timechunk_swath triggers a
        # separate catalog read and keeps the full task graph alive in worker memory.
        _ds = _ds.load()
        gc.collect()
        
        # Process this time chunk (all times in the chunk)
        timestep_results = process_timechunk_swath(_ds, verbose=verbose)
        # Not a mask variable (the zarr writer only takes the names in mask_variables): what the removal of rain
        # at missing-Tb pixels took out of this window, summed up by stream_process_to_zarr().
        timestep_results['tb_gap'] = tb_gap

        # Explicit cleanup
        del _ds
        gc.collect()
        
        return str(out_time_val), timestep_results
        
    except Exception as e:
        print(f"Error processing time chunk {start_idx}-{end_idx}: {e}")
        traceback.print_exc()
        return None, {'error': f"{type(e).__name__}: {e}"}

#--------------------------------------------------------------------------------------------------
def is_result_degenerate(result):
    """
    Detect a "silently failed" chunk result: mcs_mask AND cloud_types both entirely
    zero across every cell.

    This is a defensive, symptom-level check rather than a root-cause fix. Audit found
    that a small, non-deterministic fraction of chunks in a full parallel production run
    come back from process_timechunk_wrapper_zarr() as a "successful" future (no
    exception, no None result - future.result() returns cleanly) but with every array
    zeroed out, even though the exact same computation, re-run standalone outside the
    live Dask cluster, always produces correct, real coverage for the same inputs. Two
    distinct failure signatures were observed empirically: a confirmed Dask worker
    restart (distributed.nanny "Restarting worker") correlating with a large *contiguous*
    block of bad frames in one run, and, in a separate rerun with zero worker restarts,
    two isolated bad frames (including the very first chunk processed, ruling out any
    "accumulated over a long run" explanation on its own) - i.e. more than one underlying
    trigger can produce this same symptom. Rather than chase every possible transient
    cause individually, this check catches the shared symptom directly: real precipitation
    covers the globe somewhere in any genuine 6-hourly window, so mcs_mask and
    cloud_types both being zero at every one of the ~786k cells (in dc_pr/st_pr/nd_pr/
    dz_pr terms: the whole globe classified as having had literally zero precipitation
    for the entire window) is not a physically plausible result - it is the fingerprint
    of a lost/corrupted task, not real data.

    Parameters:
    -----------
    result : dict
        A results dict as returned by process_timechunk_swath() / the second element of
        process_timechunk_wrapper_zarr()'s return tuple. Expected keys include 'mcs_mask'
        and 'cloud_types' (numpy arrays).

    Returns:
    --------
    bool : True if the result looks like a lost/corrupted chunk (retry-worthy), False if
        it looks like genuine data (including a legitimate rare case with no MCS/cloud-
        type coverage in just one of the two fields, which alone is not suspicious).
    """
    if not result:
        return False
    mcs_mask = result.get('mcs_mask')
    cloud_types = result.get('cloud_types')
    if mcs_mask is None or cloud_types is None:
        return False
    return bool(np.all(mcs_mask == 0) and np.all(cloud_types == 0))

# Retries after the first attempt at a window, and the wait before retry n (n * this many seconds).
CHUNK_MAX_RETRIES = 3
CHUNK_RETRY_BACKOFF_SECONDS = 2.0

def resolve_chunk_result(first_attempt, run_again, label, logger,
                         max_retries=CHUNK_MAX_RETRIES, backoff_seconds=CHUNK_RETRY_BACKOFF_SECONDS):
    """
    Get a usable result for one aggregation window, retrying failures that may be transient.

    A window can fail in several ways: the task raises (unreadable zarr chunk, a netCDF/HDF5 handle
    error while reading the TC tracks, out of memory), the worker is killed, or the task "succeeds" with
    every array zero (see is_result_degenerate). A window whose input has no valid pr and no valid Tb at
    all is not a failure: model output can be missing at the start of a run (spin-up), depending on how
    it was post-processed, so it is reported as 'empty' without retries.

    Parameters:
    -----------
    first_attempt : callable
        Returns (time_str, result) for the attempt already under way (in parallel mode: waits on that
        task's future). May raise, e.g. when the worker died.
    run_again : callable
        Starts a fresh attempt and returns (time_str, result). It must really recompute. With Dask that
        means client.submit(..., pure=False): the default pure=True gives an identical call the same key,
        so the scheduler hands back the first attempt's cached result or exception and nothing is retried.
    label : str
        Names the window in log messages.
    logger : logging.Logger

    Returns:
    --------
    (result, status, detail) with status
        'ok'     : result is usable; detail says how many attempts it took when more than one
        'empty'  : the window has no valid pr or Tb; nothing to write, not a failure
        'failed' : every attempt raised, returned nothing, or stayed all-zero although the input has
                   data; result is None and detail gives the last reason
    """
    reason = "no attempt was made"
    for attempt in range(1 + max_retries):
        try:
            time_str, result = first_attempt() if attempt == 0 else run_again()
        except Exception as exc:
            result, reason = None, f"{type(exc).__name__}: {exc}"
        else:
            if result is None or time_str is None:
                reason = (result or {}).get('error') or "the task returned no result"
            elif result.get('input_all_nan'):
                return result, 'empty', "no valid pr or Tb in this window"
            elif is_result_degenerate(result):
                reason = "mcs_mask and cloud_types both entirely zero although the input has valid data"
            else:
                return result, 'ok', (f"needed {attempt + 1} attempts" if attempt else "")
        if attempt < max_retries:
            logger.warning(f"{label}: attempt {attempt + 1} failed ({reason}); retrying ({attempt + 1}/{max_retries})...")
            time.sleep(backoff_seconds * (attempt + 1))
    return None, 'failed', reason

#--------------------------------------------------------------------------------------------------
def stream_process_to_zarr(time_coords, mask_variables, output_path,
                          time_groups, output_time_coords,
                          client=None, logger=None, parallel=True,
                          input_zarr_path=None, batch_size=100,
                          config=None, store_offsets=None):
    """
    Stream process time chunks and write results to zarr with optional parallel processing.

    Uses a batched submission approach to avoid overwhelming the Dask scheduler with
    too many tasks at once. Instead of submitting all chunks at once, submits them
    in batches (super-chunks), waits for each batch to complete, then moves to the next.

    Each window is resolved by resolve_chunk_result(): failures that may be transient are retried,
    and every window ends up in exactly one of written / empty / failed, which is reported at the end.

    Parameters:
    -----------
    time_coords : array-like
        Array of input time coordinates to process
    mask_variables : list
        List of mask variable names to create
    output_path : str
        Path to output zarr file
    time_groups : dict
        Dictionary mapping aligned output times to lists of input time indices
    output_time_coords : array-like
        Array of aligned output time coordinates
    client : dask.distributed.Client, optional
        Dask client for parallel processing
    logger : logging.Logger, optional
        Logger instance
    parallel : bool
        Whether to use parallel processing
    input_zarr_path : str
        Path to input zarr file for workers to read from
    batch_size : int
        Number of output chunks to submit per batch (default: 100)
    config : dict
        Configuration dictionary with catalog and variable information

    Returns:
    --------
    dict : Summary of the run, with
        'written' : number of output frames written
        'total'   : number of output frames
        'failed'  : list of {'chunk', 'time', 'n_steps', 'reason'} for frames that could not be computed
                    or written even after retries; these read back as NaN and need a rerun
        'empty'   : same records for windows with no valid pr or Tb (expected, e.g. model spin-up);
                    they are left as NaN
        'partial' : same records for windows written from fewer hourly steps than a full window
                    (the mean is over the steps present)
        'tb_gap'  : rain removed at pixels with missing Tb, summed over the written windows
                    ('px_hours', 'rain_removed', 'rain_total', 'windows'; see remove_pr_where_tb_missing)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Time-aligned processing using time groups
    total_chunks = len(output_time_coords)
    output_times_list = list(output_time_coords)
    chunk_indices = list(range(total_chunks))
    full_steps = max(len(v) for v in time_groups.values())  # hourly steps in a complete window

    logger.info(f"Processing {len(time_coords)} input time steps into {total_chunks} aligned time groups")
    logger.info(f"Each chunk aggregates multiple hourly time steps into 1 swath mask at standard hours")

    written = np.zeros(total_chunks, dtype=bool)
    failed, empty, partial = [], [], []
    # Rain removed at pixels with missing Tb, summed over the written windows (see remove_pr_where_tb_missing)
    tb_gap_total = {'px_hours': 0, 'rain_removed': 0.0, 'rain_total': 0.0, 'windows': 0}

    def chunk_meta(chunk_idx):
        aligned_time = output_times_list[chunk_idx]
        # Convert numpy datetime64 to pandas Timestamp for dictionary lookup
        time_indices = time_groups[pd.Timestamp(aligned_time)]
        return {
            'chunk_idx': chunk_idx,
            'start_idx': int(time_indices[0]),
            'end_idx': int(time_indices[-1]) + 1,
            'n_steps': len(time_indices),
            'output_time': aligned_time,
        }

    def finish_chunk(meta, first_attempt, run_again):
        """Resolve one window (retrying transient failures), write it, and record what happened."""
        chunk_idx = meta['chunk_idx']
        label = f"Chunk {chunk_idx + 1} (time steps {meta['start_idx']}-{meta['end_idx'] - 1})"
        when = str(pd.Timestamp(meta['output_time']))[:16]
        record = {'chunk': chunk_idx, 'time': when, 'n_steps': meta['n_steps']}

        result, status, detail = resolve_chunk_result(first_attempt, run_again, label, logger)
        if status == 'empty':
            empty.append({**record, 'reason': detail})
            logger.info(f"{label} at {when}: {detail}; left as NaN (expected, e.g. model spin-up or missing input)")
            return
        if status == 'failed':
            failed.append({**record, 'reason': detail})
            logger.error(f"{label} at {when}: FAILED after {1 + CHUNK_MAX_RETRIES} attempts ({detail}); "
                         f"the frame stays NaN")
            return

        # The store's frame for this window is identified by the aligned output time, not by the first
        # hourly step the worker saw: a partial window (missing or late first hour) has a different first
        # step, and a result keyed by it would not be found by the writer and would be dropped as "missing".
        output_time = np.array([meta['output_time']], dtype='datetime64[ns]')
        try:
            logger.info(f"Writing chunk {chunk_idx + 1} to zarr...")
            # output_time holds exactly this chunk's single aggregated timestamp,
            # so its index into the store is chunk_idx itself (one chunk -> one
            # output frame here, unlike Step 3's multi-frame chunks).
            missing = append_chunk_to_zarr(
                chunk_results={str(output_time[0]): result},
                chunk_times=output_time,
                time_start=chunk_idx,
                mask_variables=mask_variables,
                output_path=output_path,
                logger=logger
            )
        except Exception as exc:
            traceback.print_exc()
            failed.append({**record, 'reason': f"write error: {type(exc).__name__}: {exc}"})
            logger.error(f"Error writing chunk {chunk_idx + 1} to zarr: {exc}")
            return
        if missing:
            failed.append({**record, 'reason': "result was not found by the writer"})
            logger.error(f"{label} at {when}: the writer found no result for {missing}")
            return

        written[chunk_idx] = True
        gap = result.get('tb_gap')
        if gap:
            tb_gap_total['px_hours'] += gap['px_hours']
            tb_gap_total['rain_removed'] += gap['rain_removed']
            tb_gap_total['rain_total'] += gap['rain_total']
            tb_gap_total['windows'] += int(gap['px_hours'] > 0)
        if meta['n_steps'] < full_steps:
            partial.append(record)
        logger.info(f"Chunk {chunk_idx + 1} complete: 1 swath mask(s) written" + (f" ({detail})" if detail else ""))

    if parallel and client is not None:
        # PARALLEL MODE: Submit chunks in batches to avoid overwhelming scheduler
        total_batches = (len(chunk_indices) + batch_size - 1) // batch_size
        logger.info(f"Using batched submission: {total_batches} batches of up to {batch_size} chunks each")

        # Prepare chunk metadata for chunks we're actually processing
        chunk_metadata = [chunk_meta(chunk_idx) for chunk_idx in chunk_indices]

        def submit(meta, **kwargs):
            return client.submit(
                process_timechunk_wrapper_zarr,
                meta['start_idx'],
                meta['end_idx'],
                input_zarr_path,
                config,
                verbose=False,
                store_offsets=store_offsets,
                **kwargs
            )

        # Process chunks in batches
        for batch_idx in range(total_batches):
            batch_start = batch_idx * batch_size
            batch_end = min((batch_idx + 1) * batch_size, len(chunk_metadata))
            batch_chunks = chunk_metadata[batch_start:batch_end]

            logger.info(f"")
            logger.info(f"{'='*80}")
            logger.info(f"BATCH {batch_idx + 1}/{total_batches}: Processing {len(batch_chunks)} chunks")
            logger.info(f"{'='*80}")

            # Submit all chunks in this batch to Dask workers (using indices only)
            futures = {}
            for meta in batch_chunks:
                futures[submit(meta)] = meta

            logger.info(f"Submitted {len(futures)} chunks to workers for this batch...")

            # Process results as they complete within this batch
            for future in as_completed(futures):
                meta = futures[future]
                logger.info(f"Processing chunk {meta['chunk_idx'] + 1}/{total_chunks}: "
                            f"time steps {meta['start_idx']}-{meta['end_idx'] - 1}")
                finish_chunk(
                    meta,
                    first_attempt=future.result,
                    # pure=False so a retry is a new task: with the default pure=True an identical
                    # submit() returns the first attempt's cached result or exception (see
                    # resolve_chunk_result).
                    run_again=lambda meta=meta: submit(meta, pure=False).result(),
                )

            logger.info(f"Batch {batch_idx + 1}/{total_batches} complete: "
                        f"{int(written.sum())}/{len(chunk_indices)} frames written so far")

    else:
        # SERIAL MODE: Process chunks one at a time
        logger.info("Processing chunks in serial mode...")

        for chunk_idx in chunk_indices:
            meta = chunk_meta(chunk_idx)
            logger.info(f"Processing chunk {chunk_idx + 1}/{total_chunks}: "
                        f"time steps {meta['start_idx']}-{meta['end_idx'] - 1}")

            def run(meta=meta):
                return process_timechunk_wrapper_zarr(
                    meta['start_idx'], meta['end_idx'], input_zarr_path, config,
                    verbose=False, store_offsets=store_offsets
                )
            finish_chunk(meta, first_attempt=run, run_again=run)

    n_written = int(written.sum())
    logger.info(f"Stream processing complete: {n_written}/{len(chunk_indices)} chunks written successfully")
    share = 100.0 * tb_gap_total['rain_removed'] / tb_gap_total['rain_total'] if tb_gap_total['rain_total'] > 0 else 0.0
    logger.info(f"Rain at pixels with missing Tb was removed before classification: {tb_gap_total['px_hours']:,} pixel-hours with "
                f"rain in {tb_gap_total['windows']} window(s), {share:.3f}% of the total rain (nothing to remove is expected "
                f"for sources whose Tb comes from OLR)")
    if empty:
        logger.info(f"{len(empty)} window(s) had no valid pr or Tb and were left as NaN (expected, e.g. model "
                    f"spin-up): {', '.join(r['time'] for r in empty[:10])}" + (" ..." if len(empty) > 10 else ""))
    if partial:
        logger.info(f"{len(partial)} window(s) had fewer than {full_steps} hourly steps and were averaged over "
                    f"the steps present: {', '.join(r['time'] + ' (' + str(r['n_steps']) + ')' for r in partial[:10])}"
                    + (" ..." if len(partial) > 10 else ""))
    if failed:
        logger.error(f"{len(failed)}/{len(chunk_indices)} frame(s) were NOT written and read back as NaN: "
                     + "; ".join(f"{r['time']} [{r['reason']}]" for r in failed[:10]) + (" ..." if len(failed) > 10 else "")
                     + ". Rerun Step 1 (or those windows) before using this store.")
    return {'written': n_written, 'total': len(chunk_indices), 'failed': failed, 'empty': empty, 'partial': partial,
            'tb_gap': tb_gap_total}

#--------------------------------------------------------------------------------------------------
def main():
    """Main function to run the make MCS swath process"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process MCS swath masks')
    parser.add_argument("-c", "--config", help="yaml config file for PyFLEXTRKR MCS tracking", required=True)
    parser.add_argument('--parallel', action='store_true', default=True,
                       help='Use parallel processing with Dask (default: True)')
    parser.add_argument('--no-parallel', action='store_false', dest='parallel',
                       help='Disable parallel processing')
    parser.add_argument('--workers', type=int, default=32,
                       help='Number of Dask workers (default: 32)')
    parser.add_argument('--threads-per-worker', type=int, default=1,
                       help='Number of threads per worker (default: 1)')
    parser.add_argument('--aggregation-window', type=int, default=6,
                       help='Number of hourly time steps to aggregate into one swath mask (default: 6)')
    parser.add_argument('--batch-size', type=int, default=100,
                       help='Number of chunks to submit per batch to avoid overwhelming scheduler (default: 100)')
    parser.add_argument('--test-steps', type=int, default=None,
                       help='Number of time steps to process for testing (default: all)')
    
    args = parser.parse_args()

    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)

    start_time = time.time()
    logger.info("Starting making MCS mask swath...")

    # Get aggregation window from input argument 
    # This determines the time window size to aggregate MCS masks (e.g., 6 = aggregate 6 hourly steps into 1 swath)
    aggregation_window = args.aggregation_window
    
    # Zarr storage chunk size (number of output time steps per zarr chunk)
    # 28 = 1 week of 6-hourly data (7 days * 4 swaths/day)
    zarr_chunk_size_time = 28
    
    # Define output variables
    mask_variables = ['mcs_mask', 'cloud_types', 'dc_pr', 'st_pr', 'nd_pr', 'dz_pr', 'tot_pr']
    
    # Define attributes for each variable
    var_attrs = {
        'mcs_mask': {
            'long_name': 'MCS swath mask',
            'description': 'MCS track swath mask aggregated over time window. Each pixel labeled with MCS track number.',
            'units': '1',
            'valid_range': [0, 1000000],
            'comment': 'Track numbers are assigned by PyFLEXTRKR. Overlapping tracks resolved by highest coverage count.'
        },
        'cloud_types': {
            'long_name': 'Cloud type classification',
            'description': 'Priority-based cloud type within non-MCS cold cloud shield over time window',
            'units': '1',
            'valid_range': [0, 4],
            'flag_values': [0, 1, 2, 3, 4],
            'flag_meanings': 'no_cloud deep_convective stratiform non_deep_convective drizzle',
            'comment': 'Classification based on brightness temperature and precipitation rate. Priority ranking: 1>2>3>4. Mutually exclusive with MCS.'
        },
        'dc_pr': {
            'long_name': 'Deep convective precipitation',
            'description': 'Frequency-weighted mean precipitation for deep convective clouds',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as deep convective (type 1), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'st_pr': {
            'long_name': 'Stratiform precipitation',
            'description': 'Frequency-weighted mean precipitation for stratiform clouds',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as stratiform (type 2), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'nd_pr': {
            'long_name': 'Non-deep convective precipitation',
            'description': 'Frequency-weighted mean precipitation for non-deep convective clouds',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as non-deep convective (type 3), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'dz_pr': {
            'long_name': 'Drizzle precipitation',
            'description': 'Frequency-weighted mean precipitation for drizzle',
            'units': 'mm h-1',
            'comment': 'Mean precipitation when cell is classified as drizzle (type 4), weighted by frequency. Mutually exclusive with MCS swath.'
        },
        'tot_pr': {
            'long_name': 'Window-mean total precipitation',
            'description': 'Mean over the aggregation window of the hourly precipitation used for the cloud-type classification',
            'units': 'mm h-1',
            'comment': ('Same variable, units and time alignment as dc_pr/st_pr/nd_pr/dz_pr; not masked by the MCS swath. '
                        'Each source uses its own pr field only (frozen precipitation is ignored). Missing hourly values count as zero. '
                        'dc_pr+st_pr+nd_pr+dz_pr equals tot_pr outside the MCS swath when every hourly step is classified. '
                        'Downstream steps (monthly and extreme-precipitation analyses) use this as the total precipitation.')
        }
    }
    
    # Parallel processing configuration
    parallel = args.parallel
    n_workers = args.workers
    threads_per_worker = args.threads_per_worker
    batch_size = args.batch_size
    config_file = args.config

    # Load MCS tracking configuration for the specified source
    zoom = 8
    config = load_config(config_file)
    root_path = config.get("root_path")
    pixel_path_name = config.get("pixel_path_name")
    pixel_path = f"{root_path}{pixel_path_name}/"
    in_basename = config.get("zarr_output_presets", {}).get("healpix").get("out_filebase")
    # # Prefer dated zarr stores (e.g., *_hp8_v1_20190101.0000_20200101.0100.zarr).
    # # Pattern matches *_v1_{date}_*.zarr but not *_v1.zarr.
    # # Falls back to the canonical *_hp8_v1.zarr when no dated stores exist.
    # zarr_pattern = f"{pixel_path}{in_basename}hp{zoom}_v1_[0-9]*.zarr"
    # dated_zarr_stores = sorted(glob.glob(zarr_pattern))
    # if len(dated_zarr_stores) >= 1:
    #     logger.info(f"Found {len(dated_zarr_stores)} dated zarr store(s) to use:")
    #     for p in dated_zarr_stores:
    #         logger.info(f"  {os.path.basename(p)}")
    #     in_zarr = dated_zarr_stores if len(dated_zarr_stores) > 1 else dated_zarr_stores[0]
    # else:
    # Fall back to the canonical v1 store
    in_zarr = f"{pixel_path}{in_basename}hp{zoom}_v1.zarr"
    logger.info(f"No dated zarr stores found, using: {os.path.basename(in_zarr)}")

    # Catalog information
    catalog_file = config.get('catalog_file')
    catalog_location = config.get('catalog_location')
    catalog_source = config.get('catalog_source')
    catalog_params = config.get('catalog_params', {})
    # Update the zoom in the catalog_params to match the zoom in the mask file
    catalog_params.update({'zoom': zoom})
    # Precipitation, OLR variable names and conversion factor
    varname_precip_liq = 'pr'
    # Frozen precipitation is ignored for all sources (each model's 'pr' only, as in the 1-hourly
    # extreme-precipitation configs), so no separate snow field such as UM's 'prs' is added.
    varname_precip_ice = None
    varname_olr = 'rlut'
    pr_convert_factor = config.get('pcp_convert_factor')

    # Get source name from root path (e.g., /pscratch/sd/w/wcmca1/hackathon/mcs/scream/)
    source_name = os.path.basename(os.path.normpath(root_path))
    # Strip trailing year-range suffix (e.g., IMERGv7_2019_2021 -> IMERGv7)
    source_name = re.sub(r'_20\d{2}_20\d{2}$', '', source_name)

    # Output paths (under the pipeline data root, see src/cof_paths.py)
    out_dir = f"{data_root(logger)}mcs_masks/"
    out_basename = f"{source_name}_mcs_masks_hp{zoom}.zarr"
    out_zarr = f"{out_dir}{out_basename}"
    os.makedirs(out_dir, exist_ok=True)

    # Prepare config dictionary for processing
    processing_config = {
        'catalog_file': catalog_file,
        'catalog_location': catalog_location,
        'catalog_source': catalog_source,
        'catalog_params': catalog_params,
        'varname_olr': varname_olr,
        'varname_precip_liq': varname_precip_liq,
        'varname_precip_ice': varname_precip_ice,
        'pcp_convert_factor': pr_convert_factor,
        # TC track location for MCS-TC exclusion (see add_tc_mask_to_dataset). These keys
        # mirror combine_tracking_masks.py's (Step 2) config so both steps read the same
        # TC data for a given source. Sources without any of these set (e.g. not yet
        # backfilled into config_mcs_tbpf_*.yml) simply skip MCS-TC exclusion here and
        # fall back to Step 3's safety net, same as pre-fix behavior for that source only.
        'tc_source_zarr': config.get('tc_source_zarr'),
        'dir_te': config.get('dir_te'),
        'source_te': config.get('source_te'),
        'source_res': config.get('source_res'),
        'tc_glob_pattern': config.get('tc_glob_pattern'),
    }

    print("="*80)
    print("MAKE MCS SWATH MASK PROCESSING")
    print("="*80)
    print(f"Source: {source_name}")
    if isinstance(in_zarr, list):
        print(f"Input: {len(in_zarr)} zarr stores (concatenated by time):")
        for p in in_zarr:
            print(f"  {p}")
    else:
        print(f"Input: {in_zarr}")
    print(f"Output: {out_zarr}")
    print(f"Parallel processing: {parallel}")
    if parallel:
        print(f"Workers: {n_workers}, Threads per worker: {threads_per_worker}")
        print(f"Batch size: {batch_size} chunks per batch")

    # Setup Dask client
    client = setup_dask_client(
        parallel=parallel, 
        n_workers=n_workers, 
        threads_per_worker=threads_per_worker, 
        logger=logger,
    )

    failed_frames = []  # windows that could not be computed or written; decides the exit status
    try:
        # Load the full dataset
        print(f"\nLoading full dataset...")
        store_offsets = None  # per-store offset table so workers open only what they need
        try:
            # Read MCS mask dataset from zarr (single store or multiple dated stores).
            # Build store_offsets alongside so workers can target individual stores.
            if isinstance(in_zarr, list):
                ds_list = []
                store_offsets = []
                offset = 0
                for p in in_zarr:
                    _s = xr.open_zarr(p, consolidated=True, mask_and_scale=True)
                    n = len(_s.time)
                    store_offsets.append((offset, offset + n, p))
                    offset += n
                    ds_list.append(_s)
                ds = xr.concat(ds_list, dim='time')
            else:
                ds = xr.open_zarr(in_zarr, consolidated=True, mask_and_scale=True)
            
            # ds = ds.pipe(egh.attach_coords)  # Commented out for testing
            print(f"  ✅ Dataset loaded successfully")
            print(f"  Time steps: {len(ds.time)}")
            print(f"  Data variables: {list(ds.data_vars)}")
            print(f"  Spatial dimensions: {dict(ds.dims)}")
        
        except Exception as e:
            print(f"  ❌ Error loading dataset: {e}")
            sys.exit(1)
        # Limit time steps for testing if requested
        if args.test_steps is not None:
            ds = ds.isel(time=slice(0, args.test_steps))
            print(f"  📋 Limited to {args.test_steps} time steps for testing")
            
        # Initialize streaming processing configuration
        print(f"\nInitializing streaming zarr processing...")
        time_coords = ds["time"].values

        # Add processing metadata
        attrs = ds.attrs.copy()
        attrs.update({
            'processing_info': 'MCS swath mask creation',
            'processing_date': str(np.datetime64('today')),
            'processing_script': os.path.basename(__file__),
        })

        # Stream processing and writing to zarr
        print(f"\nStreaming processing and writing {len(time_coords)} time steps to zarr...")
        
        print(f"Using aggregation_window={aggregation_window} hours for swath computation")

        # Create output time coordinates aligned to standard hours (00, 06, 12, 18)
        # Convert to pandas for easier time manipulation
        time_df = pd.DataFrame({'time': pd.to_datetime(time_coords)})
        
        # Round times to nearest aggregation_window hour boundary
        # For 6-hour aggregation: rounds to 00, 06, 12, 18
        time_df['time_aligned'] = time_df['time'].dt.floor(f'{aggregation_window}h')
        
        # Get unique aligned times (these are the output time coordinates)
        output_time_coords = time_df['time_aligned'].unique()
        output_time_coords = np.array(output_time_coords, dtype='datetime64[ns]')
        
        # Create mapping from aligned times to input time indices
        # This tells us which input times belong to each output time
        time_groups = time_df.groupby('time_aligned').groups
        
        logger.info(f"Input: {len(time_coords)} hourly time steps")
        logger.info(f"Output: {len(output_time_coords)} {aggregation_window}-hourly swath times (aligned to standard hours)")
        logger.info(f"First input time: {time_coords[0]}, First output time: {output_time_coords[0]}")
        logger.info(f"Last input time: {time_coords[-1]}, Last output time: {output_time_coords[-1]}")
        logger.info(f"Zarr time chunking: {zarr_chunk_size_time} time steps per chunk (optimized for weekly access)")

        try:
            # Initialize the zarr store structure (once only)
            logger.info("Initializing zarr store...")
            initialize_zarr_store(
                output_path=out_zarr,
                time_coords=output_time_coords,  # Use aggregated time coordinates
                mask_variables=mask_variables,
                template_coords=ds.coords,
                attrs=attrs,
                chunk_size_time=zarr_chunk_size_time,  # Zarr storage chunks (28 = 1 week)
                var_attrs=var_attrs  # Variable attributes
            )

            # Stream process with chunked zarr writing
            summary = stream_process_to_zarr(
                time_coords=time_coords,  # Pass input time coords for processing
                mask_variables=mask_variables,
                output_path=out_zarr,
                time_groups=time_groups,  # Mapping of aligned times to input indices
                output_time_coords=output_time_coords,  # Aligned output time coordinates
                client=client,
                logger=logger,
                parallel=parallel,
                input_zarr_path=in_zarr,  # Pass the input zarr path for workers
                batch_size=batch_size,  # Number of chunks per batch
                config=processing_config,  # Configuration dictionary
                store_offsets=store_offsets,  # Per-store offset table for targeted reads
            )
            
            total_processed = summary['written']
            failed_frames = summary['failed']
            logger.info(f"Processing finished: {total_processed}/{summary['total']} frames written to {out_zarr}"
                        + (f", {len(failed_frames)} NOT written" if failed_frames else ""))

        except Exception as e:
            logger.error(f"Error writing chunked zarr: {e}")
            traceback.print_exc()
            print(f"  ❌ Error writing zarr: {e}")
            sys.exit(1)

        # Calculate total processing time
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Store success info to print after Dask cleanup
        success_info = {
            'elapsed_time': elapsed_time,
            'out_zarr': out_zarr,
            'total_processed': total_processed,
            'n_failed': len(failed_frames),
            'n_empty': len(summary['empty']),
        }


    finally:
        # Always cleanup client
        if client and parallel:
            # Suppress Dask shutdown messages by temporarily raising log level
            logging.getLogger('distributed').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.worker').setLevel(logging.CRITICAL)
            logging.getLogger('distributed.nanny').setLevel(logging.CRITICAL)
            
            client.close()
        
        # Print success message after Dask cleanup (so it's always visible at the end)
        if 'success_info' in locals():
            elapsed = success_info['elapsed_time']
            print(f"\n{'='*80}")
            if success_info['n_failed']:
                print(f"❌ PROCESSING FINISHED WITH {success_info['n_failed']} FRAME(S) NOT WRITTEN (they read back as NaN)")
            else:
                print(f"✅ PROCESSING COMPLETE!")
            print(f"{'='*80}")
            print(f"Output: {success_info['out_zarr']}")
            print(f"Chunks processed: {success_info['total_processed']}")
            if success_info['n_empty']:
                print(f"Windows with no valid pr or Tb (left as NaN, expected): {success_info['n_empty']}")
            print(f"Total time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
            print(f"{'='*80}\n")

    # Non-zero exit status if any frame is missing, so a job script or a chained pipeline notices.
    if failed_frames:
        sys.exit(1)

if __name__ == "__main__":
    main()