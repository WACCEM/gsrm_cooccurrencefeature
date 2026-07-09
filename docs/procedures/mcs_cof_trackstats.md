# MCS Co-occurring Feature (COF) Track Statistics Procedure

**Reference scripts:**
- `scripts/extract_mcs_cof_tracks.py` — extract per-track COF overlap flags and AR/ETC partner track IDs
- `scripts/combine_mcs_cof_trackstats_multisource.py` — combine MCS track statistics with COF overlap data across all sources

**Reference notebook:** `notebooks/plot_mcs_cof_trackstats_multisource.ipynb` — COF-type classification, box-plot, and significance-test analysis

**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This procedure asks whether an MCS's own lifecycle statistics (lifetime, size, intensity, rain
production) differ depending on the type of co-occurring feature (COF) it is associated with. It
operates on individual MCS **tracks** — i.e., on lifetime-aggregated statistics per track — rather
than on gridded, per-time-step masks as in the total- and extreme-precipitation pipelines
([monthly_precip_by_cof.md](monthly_precip_by_cof.md), [extreme_precip_by_stormtype.md](extreme_precip_by_stormtype.md)).

This is Analysis 4 of the four analyses described in [pipelines/overview.md](../pipelines/overview.md)
and reuses the COF mask zarr (Step 3 of the shared upstream pipeline;
see [cof_identification.md](cof_identification.md)) together with PyFLEXTRKR MCS track-statistics
NetCDF files, which are natively hourly rather than 6-hourly.

The pipeline has three stages:

1. **`extract_mcs_cof_tracks.py`** (per source) — scans the COF mask zarr and, for each MCS track,
   records at which 6-hourly time windows it co-occurs with AR and/or ETC. These 6-hourly COF flags
   are then mapped onto the MCS track's native 1-hourly time steps, and the AR/ETC partner track
   IDs at each co-occurring cell are extracted.
2. **`combine_mcs_cof_trackstats_multisource.py`** (all sources) — merges each source's MCS track
   statistics with its COF flags, optionally removes tropical MCS tracks, converts the result to a
   tidy per-time-step DataFrame, and concatenates all six sources into one analysis-ready parquet
   file.
3. **`plot_mcs_cof_trackstats_multisource.ipynb`** — classifies each MCS track into one of four COF
   types based on which co-occurrence category dominates its lifetime, computes per-track lifetime
   statistics, splits tracks into ocean and land subsets, and compares statistics across COF types
   using box plots and a non-parametric significance test.

Six sources are analyzed: observations (IMERGv7, displayed as **OBS**) and five GSRMs (SCREAM, ICON,
UM, NICAM, CASESM2).

---

## Summary Flowchart

```
Input: COF masks zarr (Step 3, shared pipeline) + MCS track statistics NetCDF (PyFLEXTRKR, 1-hourly)
         |
         v
[Step 1] extract_mcs_cof_tracks.py   (per source)
  └─ Scan COF zarr: record (MCS track ID, 6-h timestamp) → overlap category
  └─ Map 6-hourly COF flags onto the MCS track's 1-hourly time grid
  └─ Extract AR/ETC partner track IDs at 2-way and 3-way overlap cells
         |
         v
Output: {source}_mcs_cof_tracks_2d.nc, {source}_mcs_cof_flags.parquet,
        {source}_mcs_trackstats_cof.parquet
         |
         v
[Step 2] combine_mcs_cof_trackstats_multisource.py   (all sources)
  └─ Load & merge MCS track stats + COF 2D data per source
  └─ Optional tropical filter (remove tracks with lifetime-median |meanlat| < 20°)
  └─ Convert to a tidy per-time-step DataFrame; reduce merger/PF dimensions
  └─ Concatenate all sources
         |
         v
Output: mcs_cof_trackstats_allsources.parquet
         |
         v
[Step 3] plot_mcs_cof_trackstats_multisource.ipynb
  └─ Classify each track into 1 of 4 COF types (dominant lifetime overlap fraction)
  └─ Compute per-track lifetime statistics (max / min / sum / median / mean)
  └─ Split into ocean and land MCS subsets by lifetime land fraction
  └─ Box plots + Mann-Whitney U significance test vs. the Isolated baseline
         |
         v
Output: Box-plot figures (PDF) + significance-test tables (HTML)
```

---

## Key Steps at a Glance

- **Step 1 — COF Flag Extraction:** For each MCS track, the COF mask zarr is scanned to determine,
  at every 6-hourly time step, whether the track's footprint falls in the isolated, MCS-AR, MCS-ETC,
  or MCS-AR-ETC category. Because MCS track statistics are natively hourly, each 6-hourly COF flag is
  broadcast to the (up to six) hourly steps it temporally covers.

- **Step 1 — Partner Track ID Extraction:** At every co-occurring cell, the specific AR and/or ETC
  track ID involved is also recorded (minimum ID retained if more than one partner is present),
  enabling later linkage between an MCS and the synoptic-scale system(s) it co-occurs with.

- **Step 2 — Multi-Source Combination:** Each source's MCS track statistics and COF flags are merged
  on their shared (tracks, times) grid, converted to a tidy long-format table (one row per valid
  hourly track time step), and concatenated into a single multi-source parquet file.

- **Step 2 — Tropical Filtering:** Tracks whose lifetime-median absolute latitude falls below a
  configurable threshold (default 20°) are optionally excluded, restricting the analysis to
  extratropical MCS — the population most relevant to AR/ETC co-occurrence.

- **Step 3 — COF Type Classification:** Each track is assigned to one of four COF types (Isolated,
  MCS+AR, MCS+ETC, MCS+AR+ETC) based on which overlap category occupied the largest fraction of its
  lifetime, not merely whether it ever overlapped.

- **Step 3 — Statistical Comparison:** Per-track lifetime statistics (e.g., duration, area, minimum
  brightness temperature, rain rate) are compared across the four COF types, separately for ocean and
  land MCS subsets, using a Mann-Whitney U test against the Isolated baseline.

---

## Step 1 — Extract MCS COF Flags and Partner Track IDs

### Track ID Convention

COF mask values are 1-based (0 = no feature present); MCS track statistics use a 0-based `tracks`
dimension index. The mapping is: stats index $i$ $\leftrightarrow$ COF mask value $i+1$.

### COF Time to Hourly Mapping

`make_mcs_swath_masks.py` floors each hourly input time to the nearest $N$-hour boundary (default
$N=6$). A COF timestamp $T$ therefore represents all hourly steps in $[T,\, T+N\,\text{h})$. Each MCS
track's hourly `base_time` steps are floored the same way and matched against the set of COF
timestamps at which that track's ID appears in each overlap category, so that every hourly step
inherits the COF classification of the 6-hour window it falls within.

### Flag and Partner Variables

For each MCS track and hourly time step, four binary overlap flags are recorded from the
corresponding COF mask (`mcs_ar_overlap_mask`, `mcs_etc_overlap_mask`, `mcs_ar_etc_overlap_mask`,
`mcs_isolated_mask`):

| Flag | Meaning |
|------|---------|
| `cof_flag_mcs_ar` | MCS-AR 2-way co-occurrence present |
| `cof_flag_mcs_etc` | MCS-ETC 2-way co-occurrence present |
| `cof_flag_mcs_ar_etc` | MCS-AR-ETC 3-way co-occurrence present |
| `cof_flag_isolated` | MCS isolated (no co-occurrence) |

Flag values: `1` = overlap present, `0` = valid step with no overlap, `-1` = padding (invalid
`base_time`, i.e. beyond the track's actual lifetime).

Partner track IDs are extracted from the COF mask's paired perspective variables (e.g.
`mcs_ar_overlap_mask` for the MCS-side footprint and `ar_mcs_overlap_mask` for the AR-side footprint
of the same MCS-AR pair; see [cof_identification.md](cof_identification.md) Step 6):

| Variable | Meaning |
|----------|---------|
| `ar_tracknum` | AR track ID co-occurring with the MCS (2-way + 3-way combined) |
| `etc_tracknum` | ETC track ID co-occurring with the MCS (2-way + 3-way combined) |
| `ar_tracknum_3way` | AR track ID, restricted to 3-way MCS-AR-ETC cells only |
| `etc_tracknum_3way` | ETC track ID, restricted to 3-way MCS-AR-ETC cells only |

Partner ID values: `>0` = partner track ID (minimum ID stored if multiple partners overlap one
cell), `0` = no overlap at this valid step, `-1` = padding.

### Output

Three files are written per source to `/pscratch/sd/w/wcmca1/hackathon/cof_masks/stats/`:

| File | Content |
|------|---------|
| `{source}_mcs_cof_tracks_2d.nc` | 2D `(tracks, times)` arrays: the four flag variables, the four partner-ID variables, and `base_time` (for reference) |
| `{source}_mcs_cof_flags.parquet` | One row per MCS track: boolean `flag_*` (any overlap over the track's lifetime) and `frac_*` (fraction of valid lifetime steps with that overlap) for each of the four categories |
| `{source}_mcs_trackstats_cof.parquet` | The full MCS track statistics table (one row per track), with the scalar flag/fraction columns appended; only written when `--trackstats` is provided |

If `--trackstats` is omitted, only a simplified per-track boolean flags parquet is written (no 2D
arrays, no fractions).

---

## Step 2 — Combine Multi-Source Track Statistics and COF Overlap Data

### Loading and Merging

For each source, the MCS track-statistics NetCDF and the Step 1 `{source}_mcs_cof_tracks_2d.nc` file
are opened (with `decode_times=False`, `mask_and_scale=False` so `base_time` stays as raw float64
epoch seconds) and merged on their shared `(tracks, times)` grid. `base_time` is dropped from the COF
file before merging, since it is already present in the track-statistics file.

### Tropical Track Filtering

A track is considered tropical, and removed, if the lifetime median of its `meanlat` (ignoring
padding) has an absolute value below a configurable threshold (default **20°**):

$$\text{keep track if } \left|\text{median}_t(\text{meanlat}(t))\right| \geq 20°$$

This restricts the analysis to extratropical MCS, the population where AR/ETC co-occurrence is
physically relevant. Setting the threshold to 0 disables the filter.

### Dimension Reduction

MCS track-statistics variables have several distinct native shapes, each handled differently when
flattening to a tidy per-time-step table:

| Native shape | Handling |
|--------------|----------|
| `(tracks,)` | Broadcast to every valid time step of that track |
| `(tracks, times)` | Indexed directly at each valid (track, time) pair |
| `(tracks, times, mergers)` — `merge_cloudnumber`, `split_cloudnumber` | Discarded (never read) |
| `(tracks, times, mergers)` — `merge_ccs_area`, `split_ccs_area` | Summed over the `mergers` dimension |
| `(tracks, times, nmaxpf)` | Reduced to index 0 (the largest precipitation feature, PF) |
| Any other shape | Skipped |

Only valid time steps (as determined by each track's `track_duration`) are materialized — padded
steps are never read into memory. Four columns are added: `relative_step` (0-based step index from
track initiation), `relative_time_h`, `track_duration_h`, and `dataset` (the source's display name).

### Output

A single snappy-compressed parquet file:
`/pscratch/sd/w/wcmca1/hackathon/cof_masks/stats/mcs_cof_trackstats_allsources.parquet` — one row
per valid hourly time step per track, across all six sources, with MCS track statistics and COF
overlap flags/partner IDs combined.

---

## Step 3 — COF-Type Classification and Statistical Comparison

### Derived Variables

Before classification, three derived quantities are computed: `ccs_area_all` (cold cloud shield area
including merged and split systems, `ccs_area + merge_ccs_area + split_ccs_area`), `core_area_ratio`
(`core_area / ccs_area`), and `heavy_rain_ratio` (`total_heavyrain / total_rain`); division-by-zero
results are set to NaN.

### COF-Type Classification

For each track, the fraction of its valid lifetime hourly steps spent in each 2-way/3-way category is
computed (`frac_ar`, `frac_etc`, `frac_ar_etc` — each a count of `cof_flag_* == 1` steps divided by
`track_duration`). Each track is then assigned to exactly one of four COF types by the fraction that
dominates its lifetime, applied in this order:

| Priority | COF type | Condition |
|----------|----------|-----------|
| 1 | Isolated | All three fractions equal 0 |
| 2 | MCS+AR+ETC | `frac_ar_etc` is at least as large as both `frac_ar` and `frac_etc` (3-way wins ties with either 2-way type) |
| 3 | MCS+AR | Not 3-way, and `frac_ar > frac_etc` |
| 4 | MCS+ETC | Not 3-way, and `frac_ar ≤ frac_etc` (also the tie-break when `frac_ar == frac_etc > 0`) |

This lifetime-*dominance* classification is a deliberate design choice: it assigns each track to the
co-occurrence type it spends most of its life in, rather than to every category it ever briefly
touches. Note that "Isolated" here is re-derived from the three fractions being zero, not read
directly from the `cof_flag_isolated` column extracted in Step 1 — the two should agree by
construction, since the COF mask categories are mutually exclusive at the grid-cell level (see
[cof_identification.md](cof_identification.md)), but the notebook does not depend on that column for
this classification.

### Per-Track Lifetime Statistics

Each MCS track's 1-hourly time series is reduced to a single set of lifetime statistics, using the
aggregation most physically meaningful for each variable:

| Aggregation | Variables | Output suffix |
|-------------|-----------|----------------|
| Max | `core_area`, `cold_area`, `ccs_area`, `ccs_area_all`, `pf_area`, `pf_majoraxis`, `pf_aspectratio`, `pf_maxrainrate`, `pf_rainrate`, `pf_landfrac`, `rainrate_heavyrain` | `_lt_max` |
| Min | `corecold_mintb`, `core_meantb`, `meanlat` | `_lt_min` |
| Sum | `total_rain`, `total_heavyrain` | `_lt_sum` |
| Median | `movement_speed` | `_lt_median` |
| Mean | `core_area_ratio`, `heavy_rain_ratio`, `pf_landfrac` | `_lt_mean` |

`pf_landfrac` is aggregated both ways (`pf_landfrac_lt_max` and `pf_landfrac_lt_mean`), since both the
peak and the lifetime-average land fraction are used downstream. The resulting one-row-per-track table
(`df_lt`) is merged with `track_duration_h` and the Step 3 COF type/fraction columns.

### Ocean vs. Land Subsets

Tracks are split into two (non-exhaustive, non-exclusive) subsets using lifetime land-fraction
criteria:

| Subset | Criterion | Rationale |
|--------|-----------|-----------|
| Ocean | `pf_landfrac_lt_max` < 0.05 | Track's precipitation feature never had more than 5% of its area over land |
| Land | `pf_landfrac_lt_mean` > 0.80 | Track spent more than 80% of its lifetime predominantly over land |

Both thresholds are adjustable. Box plots and significance tests (below) are produced separately for
each subset.

### Box Plots

For each subset, a 3×3 grid of box-and-whisker plots compares nine lifetime statistics (lifetime,
maximum cold-cloud-shield area, minimum brightness temperature, median movement speed, maximum PF
major-axis length, maximum PF aspect ratio, maximum rain rate, total rain, and heavy-rain ratio)
across the four COF types, with one box per data source per type. Boxes show the median and
interquartile range with whiskers at the 5th/95th percentiles; outliers beyond the whiskers are not
plotted.

### Significance Testing

For each source and each lifetime statistic, the mean percentage difference of each non-isolated COF
type (MCS+AR, MCS+ETC, MCS+AR+ETC) relative to the Isolated baseline is computed, together with a
significance test:

| Aspect | Choice |
|--------|--------|
| Comparison | Each non-isolated COF type vs. Isolated, within the same source and land/ocean subset |
| Percent difference | $(\overline{x}_\text{type} - \overline{x}_\text{isolated}) / \|\overline{x}_\text{isolated}\| \times 100$ |
| Significance test | Mann-Whitney U (two-sided), chosen because several variables (e.g. total rain, PF area) are heavily right-skewed and not well-suited to a parametric t-test |
| Significance threshold | $\alpha = 0.05$ |
| Minimum sample size | 5 valid tracks per group; comparisons with fewer are reported as not applicable |

One exception to the percent-difference formula: for minimum brightness temperature
(`corecold_mintb_lt_min`), the difference is instead normalized by the Isolated group's 10th-90th
percentile range rather than by its mean, since Tb values are all close to 200-260 K, where a
raw percentage of the mean would understate a physically meaningful shift.

### Output

Figures and tables are written to `/global/cfs/cdirs/m1867/zfeng/hk25/figures_mcs/`:

| File pattern | Content |
|--------------|---------|
| `Boxplot_MCS_by4COFtypes_3x3_ocean.pdf` | 3×3 box-plot grid, ocean MCS subset |
| `Boxplot_MCS_by4COFtypes_3x3_land.pdf` | 3×3 box-plot grid, land MCS subset |
| `Table_MCS_sig_ocean.html` | Significance-test summary table, ocean MCS subset |
| `Table_MCS_sig_land.html` | Significance-test summary table, land MCS subset |

The HTML tables are formatted for direct opening in Microsoft Word (File → Open) to preserve table
formatting.

---

## Parameter Summary

| Parameter | Default | Description |
|-----------|---------|--------------|
| `--cof-window` | 6 h | COF aggregation window; must match `--aggregation-window` in `make_mcs_swath_masks.py` |
| `--chunk-size` | 100 | COF zarr time steps loaded per processing chunk (Step 1) |
| `--tropics-lat-threshold` | 20° | Absolute lifetime-median latitude below which a track is excluded as tropical (Step 2); 0 disables the filter |
| `--time-res-h` | 1.0 h | MCS track time resolution used to derive `relative_time_h` / `track_duration_h` (Step 2) |
| `ocean_thresh` | 0.05 | Max lifetime `pf_landfrac` defining the ocean MCS subset (Step 3, notebook) |
| `land_thresh` | 0.80 | Mean lifetime `pf_landfrac` defining the land MCS subset (Step 3, notebook) |
| Significance level $\alpha$ | 0.05 | Mann-Whitney U threshold (Step 3, notebook) |
| Minimum group size | 5 tracks | Minimum sample size for a valid significance comparison (Step 3, notebook) |

## Sources

| Source key | Display name | MCS track statistics file pattern |
|------------|---------------|-----------------------------------|
| `IMERGv7` | OBS | `mcs_tracks_final_20190101.0000_20220101.0100.nc` |
| `scream` | SCREAM | `mcs_tracks_final_20190801.0000_20200901.0000.nc` |
| `icon_d3hp003` | ICON | `mcs_tracks_final_20200102.0000_20201231.2330.nc` |
| `um_glm_n2560_RAL3p3` | UM | `mcs_tracks_final_20200201.0000_20210301.0000.nc` |
| `nicam_gl11` | NICAM | `mcs_tracks_final_20200301.0000_20210301.0000.nc` |
| `casesm2_10km_nocumulus` | CASESM2 | `mcs_tracks_final_20200301.0000_20210301.0000.nc` |
