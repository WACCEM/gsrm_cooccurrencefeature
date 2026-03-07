# Monthly Precipitation Statistics by COF Type Procedure

**Reference script:** `scripts/calc_monthly_rainmap_by_cof.py`  
**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This procedure computes monthly precipitation statistics for each co-occurrence feature (COF) category on a HEALPix grid. Starting from the COF mask dataset (output of `make_cooccurrence_masks.py`) and a precipitation source (catalog or local zarr), each grid cell is attributed to one or more feature categories and its precipitation is accumulated over the month. The output is a single NetCDF file containing, for every month and every HEALPix cell, the total precipitation, feature occurrence counts, and precipitating-hour counts for all categories.

The COF mask dataset is the final product of a four-script pipeline:
1. `make_mcs_swath_masks.py` — creates aggregated MCS swath masks and cloud type precipitation on HEALPix
2. `combine_tracking_masks.py` — merges MCS swath masks with AR, TC, and ETC NetCDF tracking files into a single `allmasks` zarr
3. `make_cooccurrence_masks.py` — reads the `allmasks` zarr and produces the `cofmasks` zarr with isolated and co-occurrence masks for all feature pairs
4. **This script** — reads the `cofmasks` zarr and computes monthly precipitation statistics

The feature categories covered are:

- **All instances** of each original feature (MCS, AR, ETC, TC), regardless of co-occurrence state
- **Isolated** features (MCS-only, AR-only, ETC-only) that do not overlap with any other feature
- **Two-way COFs:** MCS-AR, MCS-ETC, AR-ETC
- **Three-way COF:** MCS-AR-ETC
- **Non-MCS cloud types** outside all feature footprints: deep convective, stratiform, non-deep convective, drizzle

---

## Summary Flowchart

```
Input: COF mask zarr  +  precipitation source (catalog or local zarr)
         |
         v
[Step 1] Data Loading and Calendar Harmonization
  └─ Open COF mask zarr (all COF category masks + cloud type precipitation)
  └─ Load precipitation from catalog or pre-processed zarr
  └─ Convert time coordinates to matching calendar if they differ
         |
         v
[Step 2] Common Time Selection and Time Range Subsetting
  └─ Intersect time arrays of COF mask and precipitation datasets
  └─ Subset to configured start/end date range
  └─ Attach precipitation variable to COF mask dataset
         |
         v
[Step 3] Monthly Grouping
  └─ Resample dataset by calendar month (month start frequency)
  └─ Dispatch each month to process_month_chunked() (serial or parallel)
         |
         v
[Step 4] Per-Month Chunked Processing
  └─ Detect time interval (hours) from first two timestamps
  └─ Loop over day-chunks (default: 6 days) within the month
         |
         v
[Step 5] Union Mask Construction (per chunk)
  └─ Combine perspective masks for each COF pair into union 2-way mask
  └─ Combine three perspective masks into union 3-way mask
         |
         v
[Step 6] Precipitation and Count Accumulation (per chunk)
  └─ For each category: sum(pr × Δt where mask > 0) → total precipitation (mm)
  └─ For each category: sum(mask > 0) → occurrence count (hours)
  └─ For each category: sum(pr > threshold & mask > 0) → precipitating hours
  └─ Cloud types: apply all-feature exclusion mask; use pre-computed dc/st/nd/dz_pr
         |
         v
[Step 7] Write Monthly Output to NetCDF
  └─ Stack monthly results into (time, cell) arrays
  └─ Attach lat/lon coordinates and variable attributes
  └─ Write compressed NetCDF4 file
         |
         v
Output: NetCDF file with monthly precipitation, count, and precipitating-hour maps
        for all COF categories and cloud types
```

---

## Key Steps at a Glance

- **Step 1 — Data Loading:** The COF mask zarr (`{source_name}_cofmasks_hp8_v1.zarr`, output of `make_cooccurrence_masks.py`) containing all isolated and co-occurrence masks plus cloud type precipitation variables is opened alongside a matched precipitation dataset from either an intake catalog or a pre-regridded local zarr. Time coordinates are converted to a common calendar when the two sources differ.

- **Step 2 — Time Harmonization:** A common time intersection ensures that only time steps present in both the mask and precipitation datasets are used. An optional start/end datetime subset further restricts the analysis to a configured period, with robust handling of both standard and cftime calendar types.

- **Step 3 — Monthly Grouping:** The merged dataset is resampled at monthly frequency. Each month is processed independently, enabling straightforward parallelization via Dask where available.

- **Step 4 — Chunked Processing:** To limit peak memory usage, each month is processed in sub-chunks of configurable length (default: 6 days). Statistics are accumulated incrementally across chunks using running sums.

- **Step 5 — Union Mask Construction:** Because the COF identification procedure stores separate perspective masks for each feature within a pair (e.g., `mcs_ar_overlap_mask` and `ar_mcs_overlap_mask`), these are combined into a single spatial union mask before computing areal precipitation statistics. The same is done for all three perspectives of the three-way COF.

- **Step 6 — Statistics Accumulation:** Three statistics are computed for each category at each grid cell: total precipitation (mm), occurrence count (hours present), and precipitating-hour count (hours with precipitation above threshold). Cloud type statistics additionally exclude all feature footprints and use the pre-computed frequency-weighted precipitation variables (`dc_pr`, `st_pr`, `nd_pr`, `dz_pr`).

- **Step 7 — Output:** Monthly results are stacked into `(time, cell)` arrays and written to a compressed NetCDF4 file with full variable attributes, enabling direct use in downstream climatological analysis.

---

## Step 1 — Data Loading and Calendar Harmonization

Two datasets are required:

**COF mask dataset** (zarr): `{source_name}_cofmasks_hp8_v1.zarr` stored in `/pscratch/sd/w/wcmca1/hackathon/cof_masks/`. This is the output of `make_cooccurrence_masks.py`, which itself reads the combined tracking mask zarr (`{source_name}_allmasks_hp8_v1.zarr`) produced by `combine_tracking_masks.py`. The upstream pipeline is:

```
make_mcs_swath_masks.py
  └─ /hackathon/mcs_masks/{source}_mcs_masks_hp8.zarr
       └─ (input to combine_tracking_masks.py along with AR/TC/ETC NetCDF files)

combine_tracking_masks.py
  └─ /hackathon/all_masks/{source}_allmasks_hp8_v1.zarr
       └─ (input to make_cooccurrence_masks.py)

make_cooccurrence_masks.py
  └─ /hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr  ← THIS SCRIPT'S INPUT
```

The `cofmasks` zarr contains:
- Original feature masks: `mcs_mask`, `ar_mask`, `etc_mask`, `tc_mask`
- Isolated masks: `mcs_isolated_mask`, `ar_isolated_mask`, `etc_isolated_mask`
- Two-way overlap masks (both perspectives each): `mcs_ar_overlap_mask`, `ar_mcs_overlap_mask`, etc.
- Three-way overlap masks (all three perspectives): `mcs_ar_etc_overlap_mask`, `ar_mcs_etc_overlap_mask`, `etc_mcs_ar_overlap_mask`
- Cloud type masks and frequency-weighted precipitation (originating from `make_mcs_swath_masks.py`): `cloud_types`, `dc_pr`, `st_pr`, `nd_pr`, `dz_pr`

**Precipitation dataset**: loaded from an intake catalog (most GSRM sources) or from a pre-regridded local zarr (for sources not yet in the catalog, such as IMERG, SCREAM, NICAM, UM, CASESM2). Liquid and ice precipitation components are summed and converted to mm h⁻¹ using a source-specific conversion factor from `config_sources.yaml`:

| Source | Liquid var | Ice var | Convert factor | Notes |
|--------|-----------|---------|---------------|-------|
| IMERG | `precipitation` | — | 1.0 | Already mm h⁻¹ |
| SCREAM | `pr` | `prs` | 3 600 000 | m s⁻¹ → mm h⁻¹ |
| ICON | `pr` | `prs` | 3 600 | kg m⁻² s⁻¹ → mm h⁻¹ |
| NICAM | `pr` | — | 3 600 | kg m⁻² s⁻¹ → mm h⁻¹ |
| UM | `pr` | — | 3 600 | kg m⁻² s⁻¹ → mm h⁻¹ |
| CASESM2 | `pr` | `prs` | 3 600 | kg m⁻² s⁻¹ → mm h⁻¹ |
| IFS | `pr` | — | 3 600 | kg m⁻² s⁻¹ → mm h⁻¹ |

If the COF mask dataset and precipitation dataset use different calendar conventions (e.g., `datetime64` vs. `DatetimeNoLeap`), the mask dataset's time coordinate is converted to match the precipitation dataset's calendar before any time operations.

---

## Step 2 — Common Time Selection and Time Range Subsetting

A two-way intersection of time arrays selects only steps common to both datasets:

$$T_{\text{common}} = T_{\text{COF}} \cap T_{\text{pr}}$$

The precipitation variable is then added to the COF mask dataset as `ds["pr"]`. If a `start_datetime` and `end_datetime` are specified in the configuration, the merged dataset is further subsetted to that range using a robust method that handles both `numpy.datetime64` and `cftime` objects.

---

## Step 3 — Monthly Grouping and Parallel Dispatch

The merged dataset is grouped by calendar month using `xarray`'s `.resample(time='1MS')`. Each monthly group is dispatched to `process_month_chunked()` either:

- **Serially:** one month at a time in a loop, or
- **In parallel:** all months submitted simultaneously to a Dask distributed cluster, with progress tracked via `dask.distributed.progress`.

Results from all months are gathered into a list before writing output.

---

## Step 4 — Per-Month Chunked Processing

Within each month, processing is split into time sub-chunks to control memory usage:

$$\text{chunk size (steps)} = \lfloor \text{chunk\_days} \times \frac{24}{\Delta t} \rfloor$$

where $\Delta t$ is the detected time interval in hours (inferred from the first two timestamps) and `chunk_days` defaults to 6. Running sums are initialized on the first chunk and accumulated over subsequent chunks.

---

## Step 5 — Union Mask Construction

The COF identification procedure (see [cof_identification_procedure.md](cof_identification_procedure.md)) stores each feature's perspective separately within a co-occurrence pair. For statistics, the spatial union across perspectives is needed:

**Two-way union masks:**
$$M_{\text{MCS-AR}} = (M_{\text{mcs\_ar}} > 0) \cup (M_{\text{ar\_mcs}} > 0)$$
$$M_{\text{MCS-ETC}} = (M_{\text{mcs\_etc}} > 0) \cup (M_{\text{etc\_mcs}} > 0)$$
$$M_{\text{AR-ETC}} = (M_{\text{ar\_etc}} > 0) \cup (M_{\text{etc\_ar}} > 0)$$

**Three-way union mask:**
$$M_{\text{3-way}} = (M_{\text{mcs\_ar\_etc}} > 0) \cup (M_{\text{ar\_mcs\_etc}} > 0) \cup (M_{\text{etc\_mcs\_ar}} > 0)$$

The union operation is implemented by summing binary versions of each perspective mask and testing whether the sum exceeds zero. The result is a binary mask for each category covering all pixels belonging to that co-occurrence type, regardless of which feature "owns" the pixel.

---

## Step 6 — Precipitation and Count Accumulation

For each category $c$ and each time chunk with $N_t$ time steps and time interval $\Delta t$ (hours), three statistics are accumulated:

**Total precipitation (mm):**
$$P_c = \sum_{t=1}^{N_t} \text{pr}(t) \cdot \mathbf{1}[M_c(t) > 0] \cdot \Delta t$$

**Occurrence count (hours):**
$$K_c = \sum_{t=1}^{N_t} \mathbf{1}[M_c(t) > 0] \cdot \Delta t$$

**Precipitating-hour count (hours above threshold $\theta$):**
$$H_c = \sum_{t=1}^{N_t} \mathbf{1}[M_c(t) > 0 \;\text{and}\; \text{pr}(t) > \theta] \cdot \Delta t$$

where $\theta = 0.1\;\text{mm\,h}^{-1}$ by default.

### Cloud Type Statistics

For cloud types, an additional all-feature exclusion mask is applied to ensure that cloud type precipitation is counted only at cells not covered by any feature mask:

$$M_{\text{all features}} = \bigcup \{M_{\text{MCS}}, M_{\text{AR}}, M_{\text{ETC}}, M_{\text{TC}}, M_{\text{iso}}, M_{\text{2-way}}, M_{\text{3-way}}\}$$

Cloud type precipitation uses the pre-computed frequency-weighted variables (`dc_pr`, `st_pr`, `nd_pr`, `dz_pr`) from the swath mask pipeline rather than raw per-timestep precipitation, since those variables already encode the temporal aggregation structure of the 6-hourly swath windows:

$$P_{\text{dc}} = \sum_t \text{dc\_pr}(t) \cdot \mathbf{1}[\neg M_{\text{all features}}(t)] \cdot \Delta t$$

and analogously for `st_pr`, `nd_pr`, `dz_pr`. Occurrence counts for cloud types are computed as the number of time steps at which the respective `x_pr` variable is positive and no feature mask is active.

---

## Step 7 — Output Variables

All output variables have shape `(time, cell)` with a monthly time dimension and written to a compressed NetCDF4 file:

| Variable | Units | Description |
|----------|-------|-------------|
| `precipitation` | mm | Total precipitation (all sources) |
| `ntimes` | count | Number of sub-daily time steps in the month |
| `mcs_precipitation` | mm | Precipitation under any MCS footprint |
| `ar_precipitation` | mm | Precipitation under any AR footprint |
| `etc_precipitation` | mm | Precipitation under any ETC footprint |
| `tc_precipitation` | mm | Precipitation under any TC footprint |
| `mcs_iso_precipitation` | mm | Precipitation under isolated MCS |
| `ar_iso_precipitation` | mm | Precipitation under isolated AR |
| `etc_iso_precipitation` | mm | Precipitation under isolated ETC |
| `mcs_ar_precipitation` | mm | Precipitation under MCS-AR 2-way COF |
| `mcs_etc_precipitation` | mm | Precipitation under MCS-ETC 2-way COF |
| `ar_etc_precipitation` | mm | Precipitation under AR-ETC 2-way COF |
| `mcs_ar_etc_precipitation` | mm | Precipitation under MCS-AR-ETC 3-way COF |
| `dc_precipitation` | mm | Non-feature deep convective precipitation |
| `st_precipitation` | mm | Non-feature stratiform precipitation |
| `nd_precipitation` | mm | Non-feature non-deep convective precipitation |
| `dz_precipitation` | mm | Non-feature drizzle precipitation |
| `*_count` | hour | Hours of feature presence (for each category above) |
| `*_precipitation_count` | hour | Hours with precipitation above threshold |

---

## Parameter Summary

| Parameter | Default | Description |
|-----------|---------|-------------|
| Aggregation | Monthly | Output time resolution |
| `chunk_days` | 6 days | Sub-chunk size for within-month processing |
| `pcp_thresh` | 0.1 mm h⁻¹ | Precipitation threshold for precipitating-hour count |
| HEALPix zoom | 8 | Spatial resolution of input and output grids |
| Time interval $\Delta t$ | Auto-detected | Hours between consecutive time steps |
