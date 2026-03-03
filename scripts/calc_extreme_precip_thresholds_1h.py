"""
Calculate extreme precipitation percentiles from HEALPix grid data (1-hourly).

This script computes precipitation percentiles (e.g., 90th, 95th) at different
time scales (1-hourly, 6-hourly, daily, etc.) from 1-hourly data on HEALPix grid.
All input data are sourced from the catalog, except IMERG which uses local 1-hourly
Zarr files.

Author: Zhe Feng, zhe.feng@pnnl.gov
"""
import numpy as np
import sys
import os
from pathlib import Path
import yaml
import xarray as xr
import pandas as pd
import time
import argparse
import intake
import logging
import easygems.healpix as egh


def parse_cmd_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description='Calculate extreme precipitation percentiles (1-hourly input)')
    parser.add_argument('--catalog_source', type=str, required=True,
                        help='Catalog source name (e.g., scream_ne120, IR_IMERG)')
    parser.add_argument('--config_file', type=str,
                        default='/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources_1h.yaml',
                        help='Path to configuration YAML file '
                             '(default: /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/config/config_sources_1h.yaml)')
    parser.add_argument('--zoom', type=int, default=8,
                        help='HEALPix zoom level (default: 8)')
    parser.add_argument('--percentiles', type=float, nargs='+', default=[95, 99],
                        help='Percentiles to compute (default: 95 99)')
    parser.add_argument('--time_durations', type=str, nargs='+', default=['1h'],
                        help='Time durations for resampling (default: 1h)')
    parser.add_argument('--method', type=str, default='linear',
                        choices=['linear', 'lower', 'higher', 'midpoint', 'nearest'],
                        help='Quantile interpolation method (default: linear)')
    parser.add_argument('--min_precip_threshold', type=float, default=0.1,
                        help='Minimum precipitation threshold (mm/h) to exclude from percentile calculation. '
                             'Values below this threshold are treated as missing. (default: 0.1 mm/h)')
    parser.add_argument('--output_dir', type=str, default='/pscratch/sd/w/wcmca1/hackathon/extreme_precip_1h/',
                        help='Output directory for NetCDF files '
                             '(default: /pscratch/sd/w/wcmca1/hackathon/extreme_precip_1h/)')
    parser.add_argument('--version', type=str, default='v1',
                        help='Version string for output files (default: v1)')

    args = parser.parse_args()
    args_dict = vars(args)
    return args_dict


def setup_logging():
    """
    Set up logging configuration.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)


def load_config(config_file, catalog_source):
    """
    Load configuration from YAML file for a specific catalog source.

    Args:
        config_file: str
            Path to the YAML configuration file
        catalog_source: str
            The catalog source key to load configuration for

    Returns:
        dict: Configuration dictionary for the specified source
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    if catalog_source not in config:
        raise ValueError(f"Catalog source '{catalog_source}' not found in config file. "
                         f"Available sources: {list(config.keys())}")

    return config[catalog_source]


def calc_precip_percentiles(pr, percentiles=[90, 95], time_duration='1h', method='linear',
                            min_precip_threshold=None, logger=None):
    """
    Calculate precipitation percentiles at specified time scales.

    Args:
        pr: xr.DataArray
            Precipitation data array with dimensions (time, cell) in mm/h
        percentiles: list of float
            Percentiles to compute (e.g., [90, 95] for 90th and 95th percentiles)
        time_duration: str
            Time duration for resampling. Options include:
            - '1h' or '1H': 1-hourly (no resampling if already 1-hourly)
            - '6h' or '6H': 6-hourly
            - '1D' or '24h': daily
            - '12h' or '12H': 12-hourly
            Default is '1h' (keeps original 1-hourly resolution)
        method: str
            Interpolation method for quantile calculation. Options include:
            - 'linear': linear interpolation (default)
            - 'lower': lower value
            - 'higher': higher value
            - 'midpoint': midpoint of two nearest values
            - 'nearest': nearest value
        min_precip_threshold: float, optional
            Minimum precipitation threshold (mm/h). Values below this threshold
            are excluded from percentile calculations (treated as NaN).
            If None, all precipitation values are included. Default is None.
        logger: logging.Logger
            Logger instance

    Returns:
        dict: Dictionary with percentile values as keys and DataArrays as values
              e.g., {90: pr_p90, 95: pr_p95}
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Apply minimum precipitation threshold if specified
    if min_precip_threshold is not None:
        logger.info(f"Applying minimum precipitation threshold: {min_precip_threshold} mm/h")
        logger.info(f"Values below threshold will be excluded from percentile calculation")
        # Replace values below threshold with NaN
        pr_filtered = pr.where(pr >= min_precip_threshold)
    else:
        logger.info("No minimum precipitation threshold applied")
        pr_filtered = pr

    # Resample to specified time duration if needed.
    # Input data are 1-hourly; skip resampling when time_duration is also 1h.
    if time_duration.lower() not in ['1h', '1hr']:
        # Convert time duration to pandas-compatible frequency string
        freq = time_duration.upper().replace('HR', 'H')
        logger.info(f"Resampling precipitation from 1-hourly to {freq}...")
        # Resample by taking mean precipitation rate (ignoring NaN)
        pr_resampled = pr_filtered.resample(time=freq).mean(keep_attrs=True)
    else:
        logger.info(f"Using original 1-hourly precipitation data...")
        pr_resampled = pr_filtered

    # Calculate percentiles
    results = {}
    for p in percentiles:
        logger.info(f"Computing {p}th percentile...")
        quantile_value = p / 100.0
        # skipna=True by default in quantile(), so NaN values are ignored
        pr_percentile = pr_resampled.quantile(quantile_value, dim='time', method=method, skipna=True)
        pr_percentile.name = f'pr_p{p}'
        pr_percentile.attrs['long_name'] = f'{p}th percentile precipitation'
        pr_percentile.attrs['units'] = 'mm/h'
        pr_percentile.attrs['time_duration'] = time_duration
        pr_percentile.attrs['method'] = method
        if min_precip_threshold is not None:
            pr_percentile.attrs['min_precip_threshold'] = f'{min_precip_threshold} mm/h'
            pr_percentile.attrs['note'] = f'Precipitation below {min_precip_threshold} mm/h excluded from calculation'
        results[p] = pr_percentile

    return results


def write_netcdf(results_dict, ds_p, output_filename, zoom, source_name,
                 start_datetime, end_datetime, time_duration, method,
                 min_precip_threshold=None, logger=None):
    """
    Write precipitation percentiles to a NetCDF file.

    Args:
        results_dict: dict
            Dictionary with percentile values as keys and DataArrays as values
        ds_p: xr.Dataset
            Original precipitation dataset (for coordinates)
        output_filename: str
            Output NetCDF filename
        zoom: int
            HEALPix zoom level
        source_name: str
            Data source name
        start_datetime: str
            Start date/time of data
        end_datetime: str
            End date/time of data
        time_duration: str
            Time duration used for resampling
        method: str
            Quantile interpolation method used
        min_precip_threshold: float, optional
            Minimum precipitation threshold used
        logger: logging.Logger
            Logger instance
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info(f'Preparing data for output file: {output_filename}')

    # Create variables dictionary
    var_dict = {}
    for p, data in results_dict.items():
        var_name = f'pr_p{int(p)}'
        var_dict[var_name] = (['cell'], data.values)

    # Create coordinates
    coord_dict = {
        'cell': (['cell'], ds_p['cell'].values),
        'lat': (['cell'], ds_p['lat'].values),
        'lon': (['cell'], ds_p['lon'].values),
    }

    # Add crs if available
    if 'crs' in ds_p:
        coord_dict['crs'] = ds_p['crs'].values

    # Create global attributes
    gattr_dict = {
        'Title': 'Precipitation percentiles at different time scales',
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov',
        'source_name': source_name,
        'start_date': start_datetime,
        'end_date': end_datetime,
        'created_on': time.ctime(time.time()),
        'grid_type': 'HEALPix',
        'zoom_level': zoom,
        'time_duration': time_duration,
        'quantile_method': method,
    }

    if min_precip_threshold is not None:
        gattr_dict['min_precip_threshold'] = f'{min_precip_threshold} mm/h'
        gattr_dict['threshold_note'] = f'Precipitation below {min_precip_threshold} mm/h excluded from percentile calculation'

    # Create output dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Add coordinate attributes
    dsout['cell'].attrs['long_name'] = 'HEALPix cell index'
    dsout['lon'].attrs['long_name'] = 'Longitude'
    dsout['lon'].attrs['units'] = 'degree'
    dsout['lat'].attrs['long_name'] = 'Latitude'
    dsout['lat'].attrs['units'] = 'degree'

    # Add variable attributes
    for p in results_dict.keys():
        var_name = f'pr_p{int(p)}'
        dsout[var_name].attrs['long_name'] = f'{int(p)}th percentile precipitation'
        dsout[var_name].attrs['units'] = 'mm/h'
        dsout[var_name].attrs['time_duration'] = time_duration
        dsout[var_name].attrs['percentile'] = p
        dsout[var_name].attrs['method'] = method
        if min_precip_threshold is not None:
            dsout[var_name].attrs['min_precip_threshold'] = f'{min_precip_threshold} mm/h'

    # Save the output file
    fillvalue = np.nan
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    logger.info(f'Writing output file: {output_filename}')
    dsout.to_netcdf(path=output_filename, mode='w', format='NETCDF4', encoding=encoding)
    logger.info(f'Successfully wrote: {output_filename}')

    return dsout


def load_precipitation_data(config_file, catalog_source, zoom, logger=None):
    """
    Load 1-hourly precipitation data from catalog or direct Zarr file.

    All data sources are loaded from the HEALPix catalog, with the exception
    of IMERG which uses local 1-hourly Zarr files.

    Args:
        config_file: str
            Path to configuration YAML file
        catalog_source: str
            Catalog source name
        zoom: int
            HEALPix zoom level
        logger: logging.Logger
            Logger instance

    Returns:
        tuple: (pr DataArray, ds_p Dataset, config dict)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    # Load configuration for the specified source
    config = load_config(config_file, catalog_source)

    # Extract configuration variables
    source_name = config.get('source_name')
    catalog_location = config.get('catalog_location', 'NERSC')
    catalog_params = config.get('catalog_params', {}).copy()
    varname_precip_liq = config.get('varname_precip_liq')
    varname_precip_ice = config.get('varname_precip_ice')
    pr_convert_factor = config.get('pr_convert_factor')
    start_datetime = config.get('start_datetime')
    end_datetime = config.get('end_datetime')

    # Catalog file
    catalog_file = "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"

    if catalog_source == "IMERG":
        # IMERG 1-hourly data is not in the catalog; use local Zarr files
        dir_healpix = "/pscratch/sd/w/wcmca1/GPM/healpix/"
        in_zarr = f"{dir_healpix}IMERG_V7_1H_zoom{zoom}_20190101_20211231.zarr"
        logger.info(f"Loading IMERG 1-hourly dataset (NOT from catalog): {in_zarr}")
        ds_p = xr.open_zarr(in_zarr, consolidated=True)
        ds_p = ds_p.pipe(egh.attach_coords)

    else:
        # All other sources: load from the HEALPix catalog
        # Update the zoom level in catalog_params
        catalog_params['zoom'] = zoom
        logger.info(f"Loading HEALPix catalog: {catalog_file}")
        in_catalog = intake.open_catalog(catalog_file)
        if catalog_location:
            in_catalog = in_catalog[catalog_location]

        logger.info(f"Loading {catalog_source} from catalog with params: {catalog_params}")
        ds_p = in_catalog[catalog_source](**catalog_params).to_dask()
        # Add lat/lon coordinates to the HEALPix DataSet
        ds_p = ds_p.pipe(egh.attach_coords)
    # import pdb; pdb.set_trace()

    # Check liquid precipitation variable
    if varname_precip_liq in list(ds_p.keys()):
        # Convert liquid precipitation to mm/h
        pr = ds_p[varname_precip_liq] * pr_convert_factor

    # Check if the ice precipitation variable exists in the dataset
    if varname_precip_ice and varname_precip_ice in list(ds_p.keys()):
        # Convert ice precipitation to liquid equivalent
        prs = ds_p[varname_precip_ice] * pr_convert_factor
        # Add ice precipitation to get total precipitation
        pr = pr + prs

    logger.info(f"Precipitation data loaded: {pr.shape}")

    return pr, ds_p, config


def main():
    """
    Main function to calculate extreme precipitation percentiles from 1-hourly data.
    """
    # Parse command-line arguments
    args_dict = parse_cmd_args()

    # Setup logging
    logger = setup_logging()
    logger.info("=" * 60)
    logger.info("Starting extreme precipitation percentile calculation (1-hourly)")
    logger.info("=" * 60)

    # Extract arguments
    config_file = args_dict['config_file']
    catalog_source = args_dict['catalog_source']
    zoom = args_dict['zoom']
    percentiles = args_dict['percentiles']
    time_durations = args_dict['time_durations']
    method = args_dict['method']
    min_precip_threshold = args_dict['min_precip_threshold']
    output_dir = args_dict['output_dir']
    version = args_dict['version']

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load precipitation data
    logger.info(f"Loading 1-hourly precipitation data for {catalog_source}...")
    pr, ds_p, config = load_precipitation_data(config_file, catalog_source, zoom, logger)

    source_name = config.get('source_name')
    start_datetime = config.get('start_datetime')
    end_datetime = config.get('end_datetime')

    # Process each time duration
    for time_duration in time_durations:
        logger.info("=" * 60)
        logger.info(f"Processing time duration: {time_duration}")
        logger.info("=" * 60)

        # Calculate percentiles
        results = calc_precip_percentiles(pr, percentiles=percentiles,
                                          time_duration=time_duration,
                                          method=method,
                                          min_precip_threshold=min_precip_threshold,
                                          logger=logger)

        # Compute results (convert from dask to numpy)
        logger.info("Computing results...")
        results_computed = {}
        for p, data in results.items():
            results_computed[p] = data.compute()
            logger.info(f"P{int(p)}: min={float(results_computed[p].min()):.3f}, "
                        f"max={float(results_computed[p].max()):.3f}, "
                        f"mean={float(results_computed[p].mean()):.3f} mm/h")

        # Create output filename
        time_str = time_duration.lower().replace('h', 'h').replace('d', 'd')
        out_basename = f"{source_name}_precip_percentiles_{time_str}_hp{zoom}_{version}.nc"
        output_filename = os.path.join(output_dir, out_basename)

        # Write to NetCDF
        write_netcdf(results_computed, ds_p, output_filename, zoom, source_name,
                     start_datetime, end_datetime, time_duration, method,
                     min_precip_threshold, logger)

    logger.info("=" * 60)
    logger.info("Extreme precipitation percentile calculation complete!")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
