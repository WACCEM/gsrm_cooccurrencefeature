# Combined Tracking Mask Creation Procedure

**Reference script:** `scripts/combine_tracking_masks.py`  
**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This procedure combines four independently produced feature tracking datasets—MCS (mesoscale convective systems), AR (atmospheric rivers), TC (tropical cyclones), and ETC (extratropical cyclones)—into a single spatially and temporally consistent zarr store on a HEALPix grid. The output provides a unified multi-feature mask dataset suitable for co-occurrence analysis.

Each feature type is tracked by a separate algorithm and stored in different formats (zarr for MCS; NetCDF files for AR, TC, and ETC). This step harmonizes their time coordinates, spatial dimensions, and variable names before merging and writing to a common output.

---

## Summary Flowchart

```
Input: MCS zarr + AR NetCDF files + TC NetCDF files + ETC NetCDF files
         |
         v
[Step 1] Load Datasets
  └─ MCS: open zarr store directly
  └─ AR, TC, ETC: open individual NetCDF files, concatenate, sort, deduplicate
         |
         v
[Step 2] Calendar Harmonization
  └─ Check calendar type of AR dataset (reference calendar)
  └─ Convert MCS time coordinate to match AR calendar if they differ
         |
         v
[Step 3] Common Time Selection
  └─ Intersect time arrays of all four datasets
  └─ Subset each dataset to the common time period
         |
         v
[Step 4] Dimension and Coordinate Standardization
  └─ Rename 'ncol' → 'cell' if present
  └─ Drop 'lat' and 'lon' coordinate variables from all datasets
         |
         v
[Step 5] Dataset Merging and Variable Cleanup
  └─ Merge all four datasets along the cell/time dimensions
  └─ Rename variables to standard names (e.g., AR_count_index → ar_mask)
  └─ Drop unwanted variables (pr, ETC_binary_tag, sfcWind)
         |
         v
[Step 6] Rechunking and zarr Output
  └─ Compute HEALPix zoom level from nside attribute
  └─ Rechunk to (time=28, cell=12×4^zoom)
  └─ Write to zarr store
         |
         v
Output: zarr file with mcs_mask, ar_mask, tc_mask, etc_mask on common time grid
```

---

## Key Steps at a Glance

- **Step 1 — Load Datasets:** The MCS dataset is read directly from a zarr store, while AR, TC, and ETC datasets are loaded by concatenating multiple NetCDF files along the time dimension. Each multi-file dataset is sorted by time and checked for duplicate time steps, which are removed before further processing.

- **Step 2 — Calendar Harmonization:** Different GSRM outputs may use different calendar conventions (e.g., proleptic gregorian vs. standard). The AR dataset serves as the reference, and the MCS time coordinate is converted to match if necessary.

- **Step 3 — Common Time Selection:** The four datasets are intersected to retain only time steps present in all of them, ensuring that the merged output has no missing values arising from misaligned time coverage.

- **Step 4 — Dimension and Coordinate Standardization:** The spatial dimension is standardized to `cell` (renaming `ncol` if needed), and redundant lat/lon coordinate variables are dropped to avoid conflicts during merging.

- **Step 5 — Variable Merging and Renaming:** The four datasets are merged into one, variables are renamed to short standard names (`ar_mask`, `tc_mask`, `etc_mask`), and ancillary variables not needed for co-occurrence analysis are dropped.

- **Step 6 — Rechunking and Output:** The merged dataset is rechunked to match the HEALPix grid structure — one full HEALPix face per spatial chunk and 28 time steps (one week of 6-hourly data) per time chunk — then written to a zarr store for downstream analysis.

---

## Step 1 — Load Datasets

The four feature datasets are loaded separately before any harmonization:

**MCS dataset** is stored as a zarr store and opened directly:
- Contains `mcs_mask` and cloud type variables from `make_mcs_swath_masks.py`
- Opened with `mask_and_scale=False` to preserve integer track numbers

**AR, TC, and ETC datasets** are stored as collections of NetCDF files (one per month or year) and are loaded by:
1. Opening each file individually with lazy loading (`chunks={}`)
2. Concatenating all files along the `time` dimension
3. Sorting the concatenated dataset by time to ensure monotonic order
4. Detecting and removing duplicate time steps (with a warning if any are found)

---

## Step 2 — Calendar Harmonization

GSRM outputs may use non-standard calendars (e.g., `360_day`, `noleap`, `julian`). Because the AR dataset is produced alongside the ETC and TC datasets, it is used as the **reference calendar**.

- If the MCS time coordinate uses a different calendar (standard/proleptic gregorian), its time values are converted to match the AR calendar using `convert_to_matching_calendar`.
- If both already use a standard calendar, no conversion is performed.

This step is necessary because `xarray` cannot merge datasets with heterogeneous calendar types.

---

## Step 3 — Common Time Selection

A four-way intersection of time arrays is computed:

$$T_{\text{common}} = T_{\text{MCS}} \cap T_{\text{AR}} \cap T_{\text{TC}} \cap T_{\text{ETC}}$$

Each dataset is then subsetted to $T_{\text{common}}$ with `.sel(time=common_times)`. If the intersection is empty, processing stops with a warning. This step ensures the merged output has no gaps or time mismatches between feature types.

---

## Step 4 — Dimension and Coordinate Standardization

Before merging, spatial dimensions and coordinates are standardized across all four datasets:

- **Dimension renaming:** If a dataset uses `ncol` as the spatial dimension name (a convention from some GSRM outputs), it is renamed to `cell` to match the HEALPix convention used by the other datasets.
- **Coordinate dropping:** The `lat` and `lon` coordinate variables are dropped from all datasets. These will be reattached from the grid description (e.g., via `egh.attach_coords`) in downstream processing rather than embedded in the merged file, which avoids conflicts between datasets that may encode lat/lon slightly differently.

---

## Step 5 — Dataset Merging and Variable Cleanup

The four standardized datasets are merged with `xr.merge` using `combine_attrs='drop_conflicts'` and `compat='override'` to resolve any minor attribute inconsistencies.

After merging, variables are renamed and cleaned up:

| Original name | Standard name | Note |
|--------------|--------------|------|
| `AR_count_index` | `ar_mask` | AR track identification number |
| `TC_count_index` | `tc_mask` | TC track identification number |
| `ETC_int_tag` | `etc_mask` | ETC track identification number |
| `pr` | *(dropped)* | Precipitation not needed in combined mask |
| `ETC_binary_tag` | *(dropped)* | Superseded by `etc_mask` |
| `sfcWind` | *(dropped)* | Not needed for co-occurrence analysis |

The `mcs_mask` variable from the MCS dataset is retained unchanged.

---

## Step 6 — Rechunking and zarr Output

The merged dataset is rechunked to a layout optimized for the HEALPix grid and downstream time-series access:

- **Spatial chunk size:** $12 \times 4^z$ cells, where $z$ is the HEALPix zoom level derived from the `healpix_nside` attribute. This corresponds to one complete HEALPix base pixel at the given resolution, which is the natural spatial unit for HEALPix data.
- **Time chunk size:** 28 time steps, corresponding to one week of 6-hourly output ($7 \text{ days} \times 4 \text{ steps/day}$).

The rechunked dataset is then written to a consolidated zarr store using Dask for parallel I/O.

---

## Output Variables

All variables are on the HEALPix `cell` dimension with a common 6-hourly time coordinate:

| Variable | Source | Description |
|----------|--------|-------------|
| `mcs_mask` | MCS swath zarr | MCS swath mask; pixel value = MCS track number (0 = no MCS) |
| `ar_mask` | AR NetCDF | AR mask; pixel value = AR track/count index (0 = no AR) |
| `tc_mask` | TC NetCDF | TC mask; pixel value = TC track number (0 = no TC) |
| `etc_mask` | ETC NetCDF | ETC mask; pixel value = ETC integer tag (0 = no ETC) |

This combined dataset is the direct input to the co-occurrence feature (COF) identification procedure described in [cof_identification.md](cof_identification.md).
