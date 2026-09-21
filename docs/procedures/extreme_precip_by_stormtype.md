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

2. **Attribution** (`calc_stormtype_extreme_precip_spatial.py`): at each time step and each cell, test whether precipitation exceeds $\tau_p(x)$, then assign the extreme event to the storm type category it belongs to. The underlying COF masks are constructed so that each grid cell belongs to at most one category in the large majority of cases (see the note under Step 4), so this assignment is usually unambiguous; a priority order is applied as a deterministic tie-break for the residual cases where a cell could otherwise match more than one category. Counts and precipitation amounts are accumulated over all time steps and written to a NetCDF file. The precipitation used in this stage is `tot_pr`, stored in the COF mask zarr next to the masks and the cloud types (see [monthly_precip_by_cof.md](monthly_precip_by_cof.md)); `calc_monthly_rainmap_by_cof.py` applies the same priority order, so its categories add up to the total precipitation.

The attribution uses two input datasets:
- **COF mask zarr** (`{source_name}_cofmasks_hp8_v1.zarr`): output of `make_cooccurrence_masks.py`, which also carries the total precipitation `tot_pr` (mm h⁻¹, written by `make_mcs_swath_masks.py`)
- **Percentile threshold file** (`{source_name}_precip_percentiles_6h_hp8_v1.nc`): output of `calc_extreme_precip_thresholds.py`, computed from a 6-hourly precipitation product (intake catalog or pre-regridded local zarr)

**Thresholds from `tot_pr` (adopted 2026-09-20).** The attribution takes the precipitation from the COF store (`tot_pr`), so the thresholds are computed from the same field: Step 1's window-mean `tot_pr`, through `calc_extreme_precip_thresholds.py --input_zarr <data root>/mcs_masks/{source}_mcs_masks_hp8.zarr --input_var tot_pr`. IMERG is the exception: its thresholds come from the non-IR 6-hourly IMERG store (`IMERG_V7_6H_zoom8_20190101_20211231.zarr`, `--input_var precipitation`), which is the same field as `tot_pr` within 60S-60N (window labels identical, correlation 0.9998 at the same label, sum ratio 1.0000; the thresholds equal the earlier production file bit for bit) and has data at all latitudes, whereas the IR-based `tot_pr` is 0 poleward of 60 (the IR input store covers 59.87S-59.87N; poleward of that `Tb` and `precipitation` are NaN). Full-record comparison of the `tot_pr`-based thresholds with the thresholds of each source's own 6-hourly product (median ratio P90 / P95; share of cells differing by more than 10% at P90): SCREAM 1.19 / 1.24 (74%), ICON 1.01 / 1.01 (45%), NICAM 1.00 / 1.00 (19%), UM 1.00 / 1.00 (2.5%), CASESM2 1.00 / 1.00 (0%), IMERG 1.00 / 1.00 (0.3%). SCREAM differs because its 6-hourly product is built from a 3-hourly data stream whose spatial processing differs from Step 1's hourly one (see `docs/AUDIT_FINDINGS.md`), is liquid only, and its windows were offset by about 3 h from Step 1's (corrected in the new file `scream_pr6h_z8_aligned.zarr`, which the thresholds no longer use); ICON's product counted snow twice. The pipeline runner (`docs/procedures/run_cof_pipeline.md`) passes these inputs.

---

## Summary Flowchart

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 STAGE 1: calc_extreme_precip_thresholds.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Input: --input_zarr (Step 1's tot_pr; IMERG: the non-IR 6-hourly store)
       or, without it, the source's own 6-hourly precipitation (catalog or local zarr)
         |
         v
[Step 1] Load and convert precipitation to mm/h
[Step 2] Apply minimum precipitation threshold (default 0.1 mm/h; the pipeline runner passes it explicitly)
         └─ Values below threshold set to NaN (excluded from percentile)
[Step 3] Optionally resample to target time duration (default: 6h)
[Step 4] Compute quantile across all time steps at each cell, in blocks of cells
         └─ e.g., 90th, 95th percentile
         |
         v
Output: {source}_precip_percentiles_6h_hp8_v1.nc
        └─ Variables: pr_p90, pr_p95, ... (shape: cell)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 STAGE 2: calc_stormtype_extreme_precip_spatial.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Input: COF mask zarr (masks, cloud types, tot_pr)  +  percentile threshold file
         |
         v
[Step 1] Load COF masks with tot_pr and percentile thresholds (stop if tot_pr is missing)
[Step 2] Optional date range subset
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
Output: {source_name}_stormtype_spatial_{pxx}{date_suffix}.nc
        └─ Per-cell total extreme precipitation plus counts and fractions for each storm type
        └─ {pxx} is the lowercase percentile name (e.g. p90); {date_suffix} is empty unless
           --start_date/--end_date were given, in which case it is _{start}_{end}
```

---

## Key Steps at a Glance

- **Stage 1, Step 1–2 — Threshold Conditioning:** Precipitation values below a minimum threshold (default 0.1 mm/h) are excluded from the percentile calculation by setting them to NaN. This prevents near-zero drizzle values from depressing the computed extreme threshold in regions with frequent light precipitation.

- **Stage 1, Steps 3–4 — Percentile Computation:** Precipitation is optionally resampled to a target time duration (default: retain 6-hourly resolution), then the $N$th percentile is computed independently at each HEALPix cell across all time steps using `xarray.quantile` with NaN-skipping. The result is a spatially varying threshold map.

- **Stage 2, Steps 1–2 — Data Loading:** The COF mask zarr (with `tot_pr`) and the threshold file are loaded and an optional analysis date range is applied. Masks and precipitation share one time axis, so no calendar conversion or time intersection with a separate precipitation product is needed.

- **Stage 2, Steps 3–4 — Time-Step Attribution:** For each time step, cells where precipitation exceeds the local percentile threshold are identified, then each extreme cell is assigned to the storm type category it belongs to. Because the COF masks are constructed to be mutually exclusive at each grid cell in the large majority of cases, this assignment is usually determined directly by which category's mask the cell falls in, not by the order categories are checked; a priority order (highest-order COF first) is applied only to deterministically resolve the residual cases where a cell could match more than one category, ensuring every extreme event is still counted exactly once.

- **Stage 2, Steps 5–6 — Spatial Accumulation and Fraction Calculation:** Extreme event counts and precipitation amounts are summed over all time steps at each cell. Precipitation fractions are then computed as the ratio of each storm type's extreme precipitation to the total extreme precipitation, providing a normalized measure of each type's relative contribution.

---

## Stage 1 — Precipitation Percentile Threshold Computation

### Input

Two ways to give the precipitation:

- **`--input_zarr STORE`** (used by the pipeline): any Zarr store on the HEALPix grid, for any source; its zoom must match `--zoom` (checked). `--input_var` names the variable
  (default: the config's `varname_precip_liq`); when it is given the field is used as it is, with no frozen precipitation added and `--input_factor` defaulting to 1, because a named
  variable such as `tot_pr` is already in mm h⁻¹. A guard rejects an input whose domain-mean precipitation is outside 0.001-20 mm h⁻¹ (a wrong variable or factor). All time steps of the store are used unless
  `--start_time` / `--end_time` are given (the config's dates apply only to the normal input). The variable, the factor, the period and the number of frames are written into the file
  (`input_zarr`, `input_var`, `input_factor`, `frames_in_input`, `precipitation_source`).
- **Without it**: each source's own 6-hourly precipitation, from an intake catalog or a pre-regridded local zarr, with liquid and ice components summed and converted to mm h⁻¹ (the earlier way; see the status note above for why the pipeline no longer uses it).

The quantiles are taken in blocks of cells (`--cell_chunk_size`, or sized from the available memory and the store's chunk width), which bounds the memory when several sources run on one node; the values equal
xarray's `quantile` up to float32 round-off. The helpers come from `calc_extreme_precip_thresholds_1h.py`, whose way of making 6-hourly values from hourly ones (threshold the hourly values, then
average the wet hours) is not used here: it is a different quantity from `tot_pr`.

### Minimum Precipitation Filter

Before computing quantiles, values below a minimum threshold $\epsilon$ (default: $0.1\;\text{mm\,h}^{-1}$) are replaced with NaN:

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

Two datasets are loaded:
- **COF masks and precipitation**: `{source_name}_cofmasks_hp8_v1.zarr` from `/pscratch/sd/w/wcmca1/hackathon/cof_masks/`. The precipitation is the variable `tot_pr`; the script stops with a message if the store has none (rerun Steps 1-3 with the current scripts).
- **Percentile thresholds**: `{source_name}_precip_percentiles_6h_hp8_v1.nc` from `/pscratch/sd/w/wcmca1/hackathon/extreme_precip/`

### Step 2 — Time Range

An optional analysis date range is applied, consistent with the `start_datetime`/`end_datetime` values in `config_sources.yaml`. The masks and `tot_pr` are already on the same time axis.

### Step 3 — Extreme Cell Identification

At each time step $t$ and cell $x$, a cell is flagged as having extreme precipitation if:

$$E(t, x) = \mathbf{1}\bigl[P(t, x) > \tau_p(x)\bigr]$$

### Step 4 — Priority-Based Storm Type Assignment

Each extreme cell is assigned to the storm type category it belongs to. An `assigned_mask` accumulates all already-classified cells so that no cell is counted twice.

**Mutual exclusivity of storm type categories.** The COF mask categories are constructed so that, for the large majority of grid cells and time steps, a cell belongs to at most one category: each tracked feature (MCS, AR, ETC) is partitioned into isolated vs. two-way/three-way co-occurrence at the track level (see [cof_identification.md](cof_identification.md)), and non-tracked cloud types (DC/ND/ST/DZ) are explicitly re-excluded from every one of the eight tracked-feature categories immediately before classification, regardless of what the upstream masks contain. In this large majority of cases, the assignment below is determined entirely by which category's mask a cell falls in, and the order in which categories are checked has no effect on the result.

This mutual exclusivity is not, however, a mathematically guaranteed invariant of the mask-construction pipeline for the eight tracked-feature categories (priorities 1–8). Three specific, structural situations can leave a cell matching more than one candidate category:

- A feature track can independently satisfy co-occurrence criteria with two different partner feature types without those two partners themselves being linked (three-way promotion currently only triggers when an ETC track bridges an AR pairing and an MCS pairing; the same check is not made for an MCS or AR track bridging two other pairings).
- Two tracks of different feature types can physically share pixels while their *mutual* overlap fraction stays below the threshold needed to register a co-occurrence pair, leaving both classified as isolated despite the spatial overlap.
- An MCS track with tropical cyclone (TC) overlap below the MCS-TC removal threshold retains all of its pixels, including any that coincide with TC pixels; the MCS-TC test is applied once, in Step 1, and Step 3 only logs it (AR and ETC do not have this gap, since they are excluded from TC pixels unconditionally at the pixel level, not by an overlap-fraction threshold).

For the residual cells affected by these situations, the priority order below acts as a deterministic tie-break — it does not represent a scientific ranking of storm-type importance, only a fixed convention (favoring three-way, then two-way, then isolated co-occurrence structure) that ensures every extreme event is still counted exactly once. Cloud types (priorities 9–12) are unaffected by any of this: they are unconditionally excluded from all eight tracked-feature categories, so the priority order between a tracked feature and a cloud type is always inert. As an optional diagnostic, the fraction of extreme cells affected by the three situations above (cells where more than one candidate category's mask is true before priority resolution) can be quantified directly from the intermediate masks in `process_single_timestep()`, to report a concrete rate rather than a qualitative "rare."

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
| 10 | Non-deep convective (ND) | `cloud_types == 3`, no feature mask, not yet assigned |
| 11 | Stratiform (ST) | `cloud_types == 2`, no feature mask, not yet assigned |
| 12 | Drizzle (DZ) | `cloud_types == 4`, no feature mask, not yet assigned |
| 13 | Unassigned | Extreme but not classified by any above |

For cloud types, tracked storm pixels are removed from the cloud type mask before classification to enforce mutual exclusivity. The union masks for 2-way and 3-way COFs are constructed as in `calc_monthly_rainmap_by_cof.py`.

All time steps are processed in parallel using `dask.delayed`, with one delayed task per time step dispatched to a Dask worker pool.

### Step 5 — Spatial Accumulation

For each storm type $c$, counts and precipitation amounts are summed over all $N_t$ time steps:

$$\text{count}_c(x) = \sum_{t=1}^{N_t} \mathbf{1}[E(t,x) \text{ and assigned to } c]$$
$$\text{precip}_c(x) = \sum_{t=1}^{N_t} P(t,x) \cdot \mathbf{1}[E(t,x) \text{ and assigned to } c]$$

where $P(t,x)$ is `tot_pr`.

### Step 6 — Precipitation Fraction

The fraction of extreme precipitation attributed to storm type $c$ at cell $x$ is:

$$f_c(x) = \frac{\text{precip}_c(x)}{\text{precip}_{\text{total}}(x)}$$

where $\text{precip}_{\text{total}}(x) = \sum_c \text{precip}_c(x)$ is the total extreme precipitation at that cell. Cells with zero total extreme precipitation are assigned a fraction value of 0.0 in the output — including `unassigned_frac`. Downstream consumers that need a true residual/unassigned category should read `unassigned_frac` directly (rather than re-deriving it as `1 - sum(other fractions)`) and mask cells where `total_extreme_count == 0`, since at those cells every `{type}_frac` is 0.0, not undefined.

Because `tot_pr` is the precipitation the cloud types were classified with, `unassigned` is zero: 0.000% of the 60°S-60°N extreme precipitation at P90 and P95 for all six sources over the full records (IMERG included, since Step 1 removes the rain at pixels with missing Tb; before that rule it was 0.016%). In the earlier production files, which took the precipitation from a separate product, it was 0.10-0.32%.

---

## Output Variables

The attribution output has shape `(cell,)` with one file per percentile threshold, saved to `/pscratch/sd/w/wcmca1/hackathon/extreme_precip/`:

| Variable | Description |
|----------|-------------|
| `total_extreme_count` | Total count of extreme precipitation occurrences (all storm types) |
| `total_extreme_precip` | Sum of all extreme precipitation values over extreme time steps (summed over samples, not time-integrated by the sampling interval; not a physical accumulated depth) |
| `{type}_count` | Extreme event count attributed to each storm type |
| `{type}_frac` | Fraction of extreme precipitation for each storm type |

where `{type}` is one of: `mcs_isolated`, `ar_isolated`, `etc_isolated`, `tc`, `mcs_ar_2way`, `mcs_etc_2way`, `ar_etc_2way`, `mcs_ar_etc_3way`, `dc`, `nd`, `st`, `dz`, `unassigned`.

Per-type precipitation sums are accumulated internally to compute `{type}_frac`, but they are not written as separate `{type}_precip` variables by the current script.

---

## Parameter Summary

| Parameter | Default | Description |
|-----------|---------|-------------|
| HEALPix zoom | 8 | Spatial resolution |
| Percentiles | P90 | One output file per percentile |
| Time duration | 6h | Time resolution for threshold calculation |
| Min precipitation filter | 0.1 mm h⁻¹ | Exclude near-zero values from threshold computation (`--min_precip_threshold`). The value the existing threshold files record; `run_all_extreme_precip_thresholds.sh` and the pipeline runner pass it explicitly; the 1-hourly script has the same default |
| Quantile method | `linear` | Interpolation method for `xarray.quantile` |
| Dask workers | 8 | Workers for parallel time-step processing |
