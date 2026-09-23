"""Latitude domain of the ETC composites and statistics (Analysis 3).

The ETC tracks are global (the extraction keeps every point with |latitude| <= 90 - radius, so up to 70 degrees), but the COF products exist only
equatorward of 60 degrees: the MCS masks are not available poleward of that, and with them the overlap flags and the MCS/AR masks of the Step 3 products
(an ETC there is "isolated" or "AR only" by construction, whatever the weather). An ETC whose storm-relative box lies mostly poleward of 60 degrees
therefore has no reliable flag or mask, and it must not enter the composites or the statistics.

The rule is on the coverage of the box: the fraction of its rows (latitudes) that lie within |latitude| <= lat_limit. The box is +-20 degrees around the
ETC centre (161 rows at 0.25 degrees), so for a centre at latitude L in the northern hemisphere the coverage is (60 - (L - 20)) / 40:

    coverage >= 1.0  <=>  |L| <= 40      (the whole box is inside)
    coverage >= 0.8  <=>  |L| <= 48
    coverage >= 0.7  <=>  |L| <= 52      (default of the composites)
    coverage >= 0.5  <=>  |L| <= 60      (the centre is inside: the "centroid" rule, default of the statistics)

Both scripts take --lat-limit and --min-lat-coverage; a minimum coverage of 0 switches the rule off.
"""
import numpy as np

LAT_LIMIT = 60.0


def get_grid_resolution(ds):
    """
    Grid spacing (degrees) of the storm-relative box, from the global attributes lon_res and lat_res.

    The coordinates x and y of the ETC 2D stores are offsets in grid points (-80 ... 80 for a box of +-20 degrees at
    0.25 degrees), so they have to be multiplied by the spacing to get degrees. A store without these attributes is an
    error: assuming a spacing would silently give a wrong radius or a wrong box.
    """
    try:
        lon_res = float(ds.attrs['lon_res'])
        lat_res = float(ds.attrs['lat_res'])
    except KeyError as missing:
        raise ValueError(f"the dataset has no global attribute {missing}: x and y are offsets in grid points and "
                         f"need lon_res and lat_res (degrees) to be converted to degrees") from None
    if not (lon_res > 0 and lat_res > 0):
        raise ValueError(f"lon_res and lat_res must be positive, got {lon_res} and {lat_res}")
    return lon_res, lat_res


def lat_box_coverage(center_lat, y_offsets_deg, lat_res, lat_limit=LAT_LIMIT):
    """
    Fraction of the rows of the storm-relative box with |latitude| <= lat_limit.

    center_lat     : array (n,), latitude of the ETC centre in degrees
    y_offsets_deg  : array (ny,), latitude offsets of the box rows in degrees (y * lat_res)
    lat_res        : grid spacing in degrees; the centre is rounded to the grid as the extraction does (the box rows are centre + offsets)
    A NaN centre has coverage 0.
    """
    c = np.round(np.asarray(center_lat, dtype='float64') / lat_res) * lat_res
    rows = c[:, None] + np.asarray(y_offsets_deg, dtype='float64')[None, :]
    with np.errstate(invalid='ignore'):
        inside = np.abs(rows) <= lat_limit + 1e-9
    return inside.mean(axis=1)


def in_cof_domain(center_lat, y_offsets_deg, lat_res, lat_limit=LAT_LIMIT, min_coverage=0.7):
    """Boolean array (n,): True where at least min_coverage of the box lies within |latitude| <= lat_limit. min_coverage <= 0 keeps every point."""
    center_lat = np.asarray(center_lat, dtype='float64')
    if min_coverage <= 0:
        return np.ones(center_lat.shape, dtype=bool)
    return lat_box_coverage(center_lat, y_offsets_deg, lat_res, lat_limit) >= min_coverage - 1e-12


def store_in_cof_domain(ds, lat_limit=LAT_LIMIT, min_coverage=0.7):
    """
    in_cof_domain for an ETC 2D dataset: the centre is cof_lat (storm_lat if there is none), the rows are y * lat_res of the dataset.
    Returns a boolean numpy array over the time (point) dimension.
    """
    if min_coverage <= 0:
        return np.ones(ds.sizes['time'], dtype=bool)
    _, lat_res = get_grid_resolution(ds)
    name = 'cof_lat' if 'cof_lat' in ds else ('storm_lat' if 'storm_lat' in ds else None)
    if name is None:
        raise ValueError("the dataset has neither cof_lat nor storm_lat: the latitude of the ETC centre is needed for the latitude rule")
    return in_cof_domain(ds[name].values, ds['y'].values.astype('float64') * lat_res, lat_res, lat_limit, min_coverage)
