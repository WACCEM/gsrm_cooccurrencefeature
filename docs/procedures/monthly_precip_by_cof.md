# Monthly Precipitation Statistics by COF Type Procedure

**Reference script:** `scripts/calc_monthly_rainmap_by_cof.py`  
**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This procedure computes monthly precipitation statistics for each co-occurrence feature (COF) category on a HEALPix grid. Starting from the COF mask dataset (output of `make_cooccurrence_masks.py`), which also carries the total precipitation `tot_pr` that the cloud types were classified with, each grid cell and time step is attributed to exactly one feature category (by priority, see Step 6) and its precipitation is accumulated over the month. The output is a single NetCDF file containing, for every month and every HEALPix cell, the total precipitation, feature occurrence time-step counts, and precipitating time-step counts for all categories.

The COF mask dataset is the final product of a four-script pipeline:
1. `make_mcs_swath_masks.py` — creates aggregated MCS swath masks, cloud type precipitation and the window-mean total precipitation (`tot_pr`) on HEALPix
2. `combine_tracking_masks.py` — merges MCS swath masks with AR, TC, and ETC NetCDF tracking files into a single `allmasks` zarr (for IMERG, whose AR/TC/ETC masks come from ERA5, this step is `combine_era5_imerg_tracking_masks.py`)
3. `make_cooccurrence_masks.py` — reads the `allmasks` zarr and produces the `cofmasks` zarr with isolated and co-occurrence masks for all feature pairs
4. **This script** — reads the `cofmasks` zarr and computes monthly precipitation statistics

The feature categories covered are:

- **All instances** of each original feature (MCS, AR, ETC), regardless of co-occurrence state (these overlap the categories below)
- **Isolated** features (MCS-only, AR-only, ETC-only) that do not overlap with any other feature
- **Two-way COFs:** MCS-AR, MCS-ETC, AR-ETC
- **Three-way COF:** MCS-AR-ETC
- **TC** (tropical cyclone) not already assigned to one of the categories above
- **Non-MCS cloud types** outside all feature footprints: deep convective, stratiform, non-deep convective, drizzle

The twelve categories that are not "all instances" (isolated ×3, two-way ×3, three-way, TC, cloud types ×4) are mutually exclusive at every grid cell and time step and add up to the total precipitation.

---

## Summary Flowchart

```
Input: COF mask zarr (COF category masks, cloud type precipitation, total precipitation tot_pr)
         |
         v
[Step 1] Data Loading
  └─ Open COF mask zarr; stop if it has no tot_pr (rerun Steps 1-3)
         |
         v
[Step 2] Time Range Subsetting
  └─ Subset to configured start/end date range
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
  └─ Assign each grid cell and time step to at most one feature category, in priority order
  └─ For each category: sum(tot_pr × Δt where assigned) → total precipitation (mm)
  └─ For each category: sum(assigned) → occurrence count (time steps)
  └─ For each category: sum(tot_pr > threshold & assigned) → precipitating count (time steps)
  └─ Cloud types: apply all-feature exclusion mask; use pre-computed dc/st/nd/dz_pr
         |
         v
[Step 7] Budget Check and Write Monthly Output to NetCDF
  └─ Check that the twelve exclusive categories add up to the total (log per month)
  └─ Stack monthly results into (time, cell) arrays
  └─ Attach lat/lon coordinates and variable attributes
  └─ Write compressed NetCDF4 file
         |
         v
Output: NetCDF file with monthly precipitation, count, and precipitating-count maps
        for all COF categories and cloud types
```

---

## Key Steps at a Glance

- **Step 1 — Data Loading:** The COF mask zarr (`{source_name}_cofmasks_hp8_v1.zarr`, output of `make_cooccurrence_masks.py`) containing all isolated and co-occurrence masks, the cloud type precipitation variables and the total precipitation `tot_pr` is opened. No separate precipitation product is read, so the masks, the cloud types and the total share one time axis and one precipitation definition.

- **Step 2 — Time Range Subsetting:** An optional start/end datetime subset restricts the analysis to a configured period, with robust handling of both standard and cftime calendar types.

- **Step 3 — Monthly Grouping:** The merged dataset is resampled at monthly frequency. Each month is processed independently, enabling straightforward parallelization via Dask where available.

- **Step 4 — Chunked Processing:** To limit peak memory usage, each month is processed in sub-chunks of configurable length (default: 6 days). Statistics are accumulated incrementally across chunks using running sums.

- **Step 5 — Union Mask Construction:** Because the COF identification procedure stores separate perspective masks for each feature within a pair (e.g., `mcs_ar_overlap_mask` and `ar_mcs_overlap_mask`), these are combined into a single spatial union mask before computing areal precipitation statistics. The same is done for all three perspectives of the three-way COF.

- **Step 6 — Statistics Accumulation:** Each grid cell and time step is first assigned to at most one feature category by priority. Three statistics are then computed for each category: total precipitation (mm), occurrence count (number of sub-daily time steps present), and precipitating count (number of sub-daily time steps with precipitation above threshold). Cloud type statistics additionally exclude all feature footprints and use the pre-computed frequency-weighted precipitation variables (`dc_pr`, `st_pr`, `nd_pr`, `dz_pr`).

- **Step 7 — Budget Check and Output:** For each month the script reports how far the sum of the twelve exclusive categories is from the total precipitation, then stacks the monthly results into `(time, cell)` arrays and writes a compressed NetCDF4 file with full variable attributes, enabling direct use in downstream climatological analysis.

---

## Step 1 — Data Loading

One dataset is required:

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
- Total precipitation (originating from `make_mcs_swath_masks.py`): `tot_pr`

**Total precipitation (`tot_pr`)**: `make_mcs_swath_masks.py` writes `tot_pr` (mm h⁻¹), the mean over each 6-hourly aggregation window of the hourly precipitation `pr` that the cloud types are classified with. It is not masked by the MCS swath, and it uses the same denominator as `dc_pr`, `st_pr`, `nd_pr` and `dz_pr` (the number of hourly steps in the window, missing steps counting as zero). Outside MCS swaths those four therefore add up to `tot_pr`, and precipitation that has no cloud type shows up as an unattributed remainder instead of being hidden by a mismatch between two precipitation products. `combine_tracking_masks.py` and `make_cooccurrence_masks.py` carry the variable through unchanged; it is not called `pr` because `combine_tracking_masks.py` drops any variable with that name. A COF store made before this change has no `tot_pr` and is rejected: rerun Steps 1-3.

Earlier versions of this script read a separate 6-hourly precipitation product (from an intake catalog or a pre-regridded zarr) and added a frozen-precipitation field where one was configured. That product differed from the hourly precipitation Step 1 used (ICON: snow counted twice and the time label at the end of the window; SCREAM: liquid only and about 3 h offset; UM: no snow term while Step 1 added one), which left a large positive or negative remainder in the budget. It is no longer read.

Each model's own precipitation field is used as it is, and no separate frozen (snow) field is added for any model. What the field contains therefore differs between sources:

| Source | Field used in Step 1 | Contains |
|--------|---------------------|----------|
| IMERG | `precipitation` | The IMERG precipitation estimate; rain at pixels with missing Tb is removed in Step 1 (Tb has gaps, about 0.3% of the rain) |
| SCREAM | hourly `pr` (`precip_total_surf_mass_flux`) | Total including snow |
| ICON | `pr` | Total including snow (`prs` is its snow part and is not added) |
| NICAM | `pr` | Single precipitation field |
| UM | `pr` | Rain only (stratiform rainfall flux) |
| CASESM2 | `pr` | Single precipitation field |

The conversion to mm h⁻¹ uses `pcp_convert_factor` from `config_mcs_tbpf_<source>.yml`.

---

## Step 2 — Time Range Subsetting

If a `start_datetime` and `end_datetime` are specified in the configuration, the dataset is subsetted to that range using a robust method that handles both `numpy.datetime64` and `cftime` objects. No intersection with a separate precipitation product and no calendar conversion are needed, because the total precipitation is stored in the same dataset as the masks.

---

## Step 3 — Monthly Grouping and Parallel Dispatch

The dataset is grouped by calendar month using `xarray`'s `.resample(time='1MS')`. Each monthly group is dispatched to `process_month_chunked()` either:

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

The COF identification procedure (see [cof_identification.md](cof_identification.md)) stores each feature's perspective separately within a co-occurrence pair. For statistics, the spatial union across perspectives is needed:

**Two-way union masks:**
$$M_{\text{MCS-AR}} = (M_{\text{mcs\_ar}} > 0) \cup (M_{\text{ar\_mcs}} > 0)$$
$$M_{\text{MCS-ETC}} = (M_{\text{mcs\_etc}} > 0) \cup (M_{\text{etc\_mcs}} > 0)$$
$$M_{\text{AR-ETC}} = (M_{\text{ar\_etc}} > 0) \cup (M_{\text{etc\_ar}} > 0)$$

**Three-way union mask:**
$$M_{\text{3-way}} = (M_{\text{mcs\_ar\_etc}} > 0) \cup (M_{\text{ar\_mcs\_etc}} > 0) \cup (M_{\text{etc\_mcs\_ar}} > 0)$$

The union operation is implemented by summing binary versions of each perspective mask and testing whether the sum exceeds zero. The result is a binary mask for each category covering all pixels belonging to that co-occurrence type, regardless of which feature "owns" the pixel.

---

## Step 6 — Precipitation and Count Accumulation

### Category assignment

The category masks written by `make_cooccurrence_masks.py` are not disjoint: a track can appear in two pair lists, and an "isolated" MCS, AR or ETC can lie inside the whole-object footprint of a three-way category. Adding the categories as they are would count part of the precipitation twice (2-3.5% at 60°S-60°N and 4-9% at mid-latitudes in the September 2026 reprocessing). Each grid cell and time step is therefore assigned to at most one of the eight feature categories, in the same priority order that `calc_stormtype_extreme_precip_spatial.py` uses:

MCS-AR-ETC > MCS-AR > MCS-ETC > AR-ETC > isolated MCS > isolated AR > isolated ETC > TC

The two-way and three-way masks are the union masks of Step 5; the isolated masks and the TC mask are used as stored. A cell that qualifies for several categories is counted only in the first of them, so `tc_*` means "TC not already assigned to a higher-priority category". The all-instances `mcs_*`, `ar_*` and `etc_*` variables are not exclusive and are computed from the original masks. (That the Step 3 masks overlap is a known issue of `make_cooccurrence_masks.py` to be followed up; the priority order removes its effect on the totals here but not necessarily on other uses of the masks, such as per-category counts or track statistics.)

### Accumulation

For each category $c$, with $A_c(t)$ true where $M_c(t) > 0$ and no higher-priority category holds, and each time chunk with $N_t$ time steps and time interval $\Delta t$ (hours), three statistics are accumulated:

**Total precipitation (mm):**
$$P_c = \sum_{t=1}^{N_t} \text{tot\_pr}(t) \cdot \mathbf{1}[A_c(t)] \cdot \Delta t$$

**Occurrence count (time steps):**
$$K_c = \sum_{t=1}^{N_t} \mathbf{1}[A_c(t)]$$

**Precipitating count (time steps above threshold $\theta$):**
$$H_c = \sum_{t=1}^{N_t} \mathbf{1}[A_c(t) \;\text{and}\; \text{tot\_pr}(t) > \theta]$$

where $\theta = 1\;\text{mm\,h}^{-1}$ by default (`--pcp_thresh`). If elapsed hours are needed, multiply these count variables by the detected time interval $\Delta t$ stored in the output global attributes.

### Cloud Type Statistics

For cloud types, an additional all-feature exclusion mask is applied to ensure that cloud type precipitation is counted only at cells not covered by any feature mask:

$$M_{\text{all features}} = \bigcup \{M_{\text{MCS}}, M_{\text{AR}}, M_{\text{ETC}}, M_{\text{TC}}, M_{\text{iso}}, M_{\text{2-way}}, M_{\text{3-way}}\}$$

Cloud type precipitation uses the pre-computed frequency-weighted variables (`dc_pr`, `st_pr`, `nd_pr`, `dz_pr`) from the swath mask pipeline rather than raw per-timestep precipitation, since those variables already encode the temporal aggregation structure of the 6-hourly swath windows:

$$P_{\text{dc}} = \sum_t \text{dc\_pr}(t) \cdot \mathbf{1}[\neg M_{\text{all features}}(t)] \cdot \Delta t$$

and analogously for `st_pr`, `nd_pr`, `dz_pr`. Occurrence counts for cloud types are computed as the number of time steps at which the respective `x_pr` variable is positive and no feature mask is active.

### Budget check

`tot_pr` is the precipitation the cloud types were classified with and the eight feature categories are exclusive, so the twelve exclusive categories (isolated MCS, AR and ETC; the three two-way COFs; the three-way COF; TC; and the four cloud types) add up to `precipitation` wherever every window has a cloud type. For each month the script logs how far the worst wet cell (more than 1 mm) and the domain total are from that sum, and warns when a cell differs by more than 0.1%. A remainder means precipitation in windows or places without a cloud-type classification, typically because the brightness temperature is missing. In the test months the remainder is zero to numerical precision for SCREAM, ICON, UM, NICAM and CASESM2. For IMERG (February 2020) it is 0.24% of the 60°S-60°N precipitation, from windows without IR brightness temperature (2020-02-26 06 and 12 UTC, and parts of 2020-02-04 and 2020-02-26 18 UTC).

---

## Step 7 — Output Variables

All output variables have shape `(time, cell)` with a monthly time dimension and written to a compressed NetCDF4 file:

| Variable | Units | Description |
|----------|-------|-------------|
| `precipitation` | mm | Total precipitation (`tot_pr` integrated over the month) |
| `ntimes` | count | Number of sub-daily time steps in the month |
| `mcs_precipitation` | mm | Precipitation under any MCS footprint (all instances, not exclusive) |
| `ar_precipitation` | mm | Precipitation under any AR footprint (all instances, not exclusive) |
| `etc_precipitation` | mm | Precipitation under any ETC footprint (all instances, not exclusive) |
| `tc_precipitation` | mm | Precipitation under TC not already assigned to a higher-priority category |
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
| `*_count` | count | Number of sub-daily time steps with feature presence (for each category above) |
| `*_precipitation_count` | count | Number of sub-daily time steps with precipitation above threshold |

The twelve exclusive categories (`mcs_iso`, `ar_iso`, `etc_iso`, `tc`, `mcs_ar`, `mcs_etc`, `ar_etc`, `mcs_ar_etc`, `dc`, `nd`, `st`, `dz`) add up to `precipitation`; their `*_count` and `*_precipitation_count` variables are exclusive in the same way. The output attributes `precipitation_source` and `category_assignment` record this, and the eight feature categories carry `assignment = priority-exclusive`.

---

## Parameter Summary

| Parameter | Default | Description |
|-----------|---------|-------------|
| Aggregation | Monthly | Output time resolution |
| `chunk_days` | 6 days | Sub-chunk size for within-month processing |
| `pcp_thresh` | 1 mm h⁻¹ | Precipitation threshold for precipitating-hour count |
| Category assignment | Priority-exclusive | MCS-AR-ETC > MCS-AR > MCS-ETC > AR-ETC > isolated MCS > AR > ETC > TC (Step 6) |
| HEALPix zoom | 8 | Spatial resolution of input and output grids |
| Time interval $\Delta t$ | Auto-detected | Hours between consecutive time steps |
