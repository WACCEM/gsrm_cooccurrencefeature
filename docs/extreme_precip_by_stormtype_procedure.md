# Extreme Precipitation by Storm Type Procedure

**Reference scripts:**
- `scripts/calc_extreme_precip_thresholds.py` — compute per-cell precipitation percentile thresholds
- `scripts/calc_stormtype_extreme_precip_spatial.py` — attribute extreme precipitation to storm types

**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

These two scripts together produce spatial maps of extreme precipitation attributed to each storm type (COF category). The analysis is split into two stages:

1. **Threshold computation** (`calc_extreme_precip_thresholds.py`): for each HEALPix cell, compute the $N$th-percentile precipitation over the full analysis period. This yields a spatially varying threshold $\tau_p(x)$ that defines what counts as "extreme" locally.

2. **Attribution** (`calc_stormtype_extreme_precip_spatial.py`): at each time step and each cell, test whether precipitation exceeds $\tau_p(x)$, then assign the extreme event to exactly one storm type following a strict priority hierarchy. Counts and precipitation amounts are accumulated over all time steps and written to a NetCDF file.

The analysis uses three independent input datasets:
- **COF mask zarr** (`{source_name}_cofmasks_hp8_v1.zarr`): output of `make_cooccurrence_masks.py`
- **Percentile threshold file** (`{source_name}_precip_percentiles_6h_hp8_v1.nc`): output of `calc_extreme_precip_thresholds.py`
- **Precipitation dataset**: intake catalog or pre-regridded local zarr (same sources as the monthly precipitation script)

---

## Summary Flowchart

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 STAGE 1: calc_extreme_precip_thresholds.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Input: precipitation (catalog or local zarr)
         |
         v
[Step 1] Load and convert precipitation to mm/h
[Step 2] Apply minimum precipitation threshold (default: 0.01 mm/h)
         └─ Values below threshold set to NaN (excluded from percentile)
[Step 3] Optionally resample to target time duration (default: 6h)
[Step 4] Compute quantile across all time steps at each cell
         └─ e.g., 90th, 95th percentile
         |
         v
Output: {source}_precip_percentiles_6h_hp8_v1.nc
        └─ Variables: pr_p90, pr_p95, ... (shape: cell)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 STAGE 2: calc_stormtype_extreme_precip_spatial.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Input: COF mask zarr  +  percentile threshold file  +  precipitation
         |
         v
[Step 1] Load COF masks, percentile thresholds, precipitation
[Step 2] Align calendars; select common time steps; optional date range subset
         |
         v
[Step 3] For each percentile (P90, P95, ...):
         |
         v
[Step 4] Process all time steps in parallel (Dask)
  └─ Per time step:
       ├─ Build union masks for all 2-way and 3-way COF pairs
       ├─ Identify extreme cells: pr > τ_p(x)
       └─ Assign each extreme cell to exactly one storm type (priority order):
            3-way COF > MCS-AR 2-way > MCS-ETC 2-way > AR-ETC 2-way
            > isolated MCS > isolated AR > isolated ETC > TC
            > cloud types (DC > ND > ST > DZ) > unassigned
         |
         v
[Step 5] Accumulate counts and precipitation amounts over all time steps
[Step 6] Compute precipitation fractions per storm type
         |
         v
Output: {source}_stormtype_extreme_precip_{Pxx}_hp8_v1.nc
        └─ Per-cell counts, precipitation sums, and fractions for each storm type
```

---

## Key Steps at a Glance

- **Stage 1, Step 1–2 — Threshold Conditioning:** Precipitation values below a minimum threshold (default 0.01 mm/h) are excluded from the percentile calculation by setting them to NaN. This prevents near-zero drizzle values from depressing the computed extreme threshold in regions with frequent light precipitation.

- **Stage 1, Steps 3–4 — Percentile Computation:** Precipitation is optionally resampled to a target time duration (default: retain 6-hourly resolution), then the $N$th percentile is computed independently at each HEALPix cell across all time steps using `xarray.quantile` with NaN-skipping. The result is a spatially varying threshold map.

- **Stage 2, Steps 1–2 — Data Loading and Alignment:** The COF mask zarr, threshold file, and precipitation dataset are loaded and aligned to a common time coordinate, handling calendar differences between model outputs and applying an optional analysis date range.

- **Stage 2, Steps 3–4 — Time-Step Attribution:** For each time step, cells where precipitation exceeds the local percentile threshold are identified, then each extreme cell is assigned to exactly one storm type using a strict priority hierarchy (highest-order COF first). This ensures that every extreme event is counted once and assigned to its most complex storm type.

- **Stage 2, Steps 5–6 — Spatial Accumulation and Fraction Calculation:** Extreme event counts and precipitation amounts are summed over all time steps at each cell. Precipitation fractions are then computed as the ratio of each storm type's extreme precipitation to the total extreme precipitation, providing a normalized measure of each type's relative contribution.

---

## Stage 1 — Precipitation Percentile Threshold Computation

### Input

Precipitation loaded identically to other scripts: from an intake catalog or a pre-regridded local zarr, with liquid and ice components summed and converted to mm h⁻¹.

### Minimum Precipitation Filter

Before computing quantiles, values below a minimum threshold $\epsilon$ (default: $0.01\;\text{mm\,h}^{-1}$) are replaced with NaN:

$$\tilde{P}(t, x) = \begin{cases} P(t, x) & \text{if } P(t, x) \geq \epsilon \\ \text{NaN} & \text{otherwise} \end{cases}$$

This ensures that the percentile is computed only over meaningful precipitation events and is not dragged down by trace amounts. NaN values are skipped during quantile computation (`skipna=True`).

### Temporal Resampling (Optional)

If a time duration other than 6h is requested (e.g., daily), precipitation is resampled by taking the mean over each duration window before computing percentiles.

### Percentile Computation

For each cell $x$, the $p$th-percentile threshold is:

$$\tau_p(x) = Q_p\bigl(\{\tilde{P}(t, x) : \tilde{P}(t, x) \text{ is not NaN}\}\bigr)$$

where $Q_p$ is the empirical quantile function. Multiple percentiles (e.g., 90th and 95th) are computed in a single pass.

### Output

NetCDF file with shape `(cell,)`, one variable per percentile: `pr_p90`, `pr_p95`, etc.

---

## Stage 2 — Storm Type Attribution of Extreme Precipitation

### Step 1 — Data Loading

Three datasets are loaded:
- **COF masks**: `{source_name}_cofmasks_hp8_v1.zarr` from `/pscratch/sd/w/wcmca1/hackathon/cof_masks/`
- **Percentile thresholds**: `{source_name}_precip_percentiles_6h_hp8_v1.nc` from `/pscratch/sd/w/wcmca1/hackathon/extreme_precip/`
- **Precipitation**: same catalog/local zarr sources as `calc_monthly_rainmap_by_cof.py`

### Step 2 — Calendar and Time Alignment

The mask dataset calendar is converted to match the precipitation dataset if they differ, then a common time intersection is taken. An optional analysis date range is applied, consistent with the `start_datetime`/`end_datetime` values in `config_sources.yaml`.

### Step 3 — Extreme Cell Identification

At each time step $t$ and cell $x$, a cell is flagged as having extreme precipitation if:

$$E(t, x) = \mathbf{1}\bigl[P(t, x) > \tau_p(x)\bigr]$$

### Step 4 — Priority-Based Storm Type Assignment

Each extreme cell is assigned to exactly one storm type following a **strict descending priority order**. An `assigned_mask` accumulates all already-classified cells so that no cell is counted twice:

| Priority | Category | Condition |
|----------|----------|-----------|
| 1 | MCS-AR-ETC (3-way) | Cell in 3-way COF union mask |
| 2 | MCS-AR (2-way) | Cell in MCS-AR union mask, not yet assigned |
| 3 | MCS-ETC (2-way) | Cell in MCS-ETC union mask, not yet assigned |
| 4 | AR-ETC (2-way) | Cell in AR-ETC union mask, not yet assigned |
| 5 | Isolated MCS | Cell in isolated MCS mask, not yet assigned |
| 6 | Isolated AR | Cell in isolated AR mask, not yet assigned |
| 7 | Isolated ETC | Cell in isolated ETC mask, not yet assigned |
| 8 | TC | Cell in TC mask, not yet assigned |
| 9 | Deep convective (DC) | `cloud_types == 1`, no feature mask, not yet assigned |
| 10 | Non-deep convective (ND) | `cloud_types == 2`, no feature mask, not yet assigned |
| 11 | Stratiform (ST) | `cloud_types == 3`, no feature mask, not yet assigned |
| 12 | Drizzle (DZ) | `cloud_types == 4`, no feature mask, not yet assigned |
| 13 | Unassigned | Extreme but not classified by any above |

For cloud types, tracked storm pixels are removed from the cloud type mask before classification to enforce mutual exclusivity. The union masks for 2-way and 3-way COFs are constructed as in `calc_monthly_rainmap_by_cof.py`.

All time steps are processed in parallel using `dask.delayed`, with one delayed task per time step dispatched to a Dask worker pool.

### Step 5 — Spatial Accumulation

For each storm type $c$, counts and precipitation amounts are summed over all $N_t$ time steps:

$$\text{count}_c(x) = \sum_{t=1}^{N_t} \mathbf{1}[E(t,x) \text{ and assigned to } c]$$
$$\text{precip}_c(x) = \sum_{t=1}^{N_t} P(t,x) \cdot \mathbf{1}[E(t,x) \text{ and assigned to } c]$$

### Step 6 — Precipitation Fraction

The fraction of extreme precipitation attributed to storm type $c$ at cell $x$ is:

$$f_c(x) = \frac{\text{precip}_c(x)}{\text{precip}_{\text{total}}(x)}$$

where $\text{precip}_{\text{total}}(x) = \sum_c \text{precip}_c(x)$ is the total extreme precipitation at that cell. Cells with zero total extreme precipitation are set to NaN.

---

## Output Variables

The attribution output has shape `(cell,)` with one file per percentile threshold, saved to `/pscratch/sd/w/wcmca1/hackathon/extreme_precip/`:

| Variable | Description |
|----------|-------------|
| `total_extreme_count` | Total count of extreme precipitation occurrences (all storm types) |
| `total_extreme_precip` | Sum of all extreme precipitation amounts (mm h⁻¹) |
| `{type}_count` | Extreme event count attributed to each storm type |
| `{type}_precip` | Extreme precipitation sum for each storm type (mm h⁻¹) |
| `{type}_fraction` | Fraction of extreme precipitation for each storm type |

where `{type}` is one of: `mcs_isolated`, `ar_isolated`, `etc_isolated`, `tc`, `mcs_ar_2way`, `mcs_etc_2way`, `ar_etc_2way`, `mcs_ar_etc_3way`, `dc`, `nd`, `st`, `dz`, `unassigned`.

---

## Parameter Summary

| Parameter | Default | Description |
|-----------|---------|-------------|
| HEALPix zoom | 8 | Spatial resolution |
| Percentiles | P90 | One output file per percentile |
| Time duration | 6h | Time resolution for threshold calculation |
| Min precipitation filter | 0.01 mm h⁻¹ | Exclude near-zero values from threshold computation |
| Quantile method | `linear` | Interpolation method for `xarray.quantile` |
| Dask workers | 8 | Workers for parallel time-step processing |
