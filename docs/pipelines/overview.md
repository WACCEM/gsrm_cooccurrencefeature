# Analysis Pipeline Overview

**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This document describes the end-to-end processing pipelines for four analyses:

1. **Total precipitation by COF type** — monthly mean precipitation maps attributed to each co-occurrence feature category
2. **Extreme precipitation by COF type** — spatial distribution of extreme precipitation events attributed to each co-occurrence feature category
3. **ETC composite analysis** — 2D composite fields and spatial statistics of atmospheric variables centered at extratropical cyclone (ETC) track points, stratified by COF overlap category
4. **MCS COF track statistics** — per-track statistics of mesoscale convective system (MCS) lifecycles stratified by COF overlap type, combining data from all sources

Analyses 1, 2, and 4 share a common upstream pipeline (Steps 1–3) that produces the COF mask dataset. Analysis 3 is an independent pipeline that operates directly on ETC track files and HEALPix model output.

---

## Shared Upstream Pipeline (Steps 1–3)

Steps 1–3 are identical for both analyses and need only be run once per model source.

```
[Step 1] MCS Swath Mask Creation
  └─ Script: scripts/make_mcs_swath_masks.py
  └─ Input:  hourly MCS pixel masks + Tb + Precipitation (HEALPix zarr, catalog)
  └─ Output: /hackathon/mcs_masks/{source}_mcs_masks_hp8.zarr
             (mcs_mask, cloud_types, dc_pr, st_pr, nd_pr, dz_pr)
         |
         v
[Step 2] Combined Tracking Mask Creation
  └─ Script: scripts/combine_tracking_masks.py
  └─ Input:  Step 1 zarr  +  AR/TC/ETC NetCDF tracking files
  └─ Output: /hackathon/all_masks/{source}_allmasks_hp8_v1.zarr
             (mcs_mask, ar_mask, tc_mask, etc_mask)
         |
         v
[Step 3] COF Mask Creation
  └─ Script: scripts/make_cooccurrence_masks.py
  └─ Input:  Step 2 zarr
  └─ Output: /hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr
             (isolated masks, 2-way and 3-way overlap masks, cloud type precipitation)
```

---

## Analysis 1 — Total Precipitation by COF Type

```
[Step 3 output] /hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr
         |
         v
[Step 4] Monthly Precipitation Statistics
  └─ Script: scripts/calc_monthly_rainmap_by_cof.py
  └─ Input:  Step 3 zarr  +  precipitation (catalog or local zarr)
  └─ Output: /hackathon/cof_masks/stats/monthly/
             {source}_monthly_rainmap_cof_hp8_v1.nc
             (monthly precipitation, count, precipitating-hour maps
              for all COF categories and cloud types)
         |
         v
[Step 5] Visualization and Analysis
  └─ Notebook: notebooks/plot_cof_raintype_rank_map.ipynb
  └─ Input:  Step 4 NetCDF
  └─ Output: Figures showing precipitation maps and rankings by COF type
```

---

## Analysis 2 — Extreme Precipitation by COF Type

```
[Step 3 output] /hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr
         |
         v
[Step 4a] Precipitation Percentile Threshold Computation
  └─ Script: scripts/calc_extreme_precip_thresholds.py
  └─ Input:  precipitation (catalog or local zarr)
  └─ Output: /hackathon/extreme_precip/
             {source}_precip_percentiles_6h_hp8_v1.nc
             (per-cell P90, P95, ... thresholds; shape: cell)

[Step 4b] Extreme Precipitation Attribution
  └─ Script: scripts/calc_stormtype_extreme_precip_spatial.py
  └─ Input:  Step 3 zarr  +  Step 4a NetCDF  +  precipitation
  └─ Output: /hackathon/extreme_precip/
             {source}_stormtype_extreme_precip_{Pxx}_hp8_v1.nc
             (per-cell counts, precipitation sums, and fractions
              for each storm type at each percentile level)
         |
         v
[Step 5] Visualization and Analysis
  └─ Notebook: notebooks/plot_cof_extreme_raintype_rank_map.ipynb
  └─ Input:  Step 4b NetCDF
  └─ Output: Figures showing extreme precipitation maps and rankings by COF type
```

> **Note:** Steps 4a and 4b are independent and can be run in either order or in parallel, but both must complete before Step 5.

---

## Analysis 3 — ETC Composite Analysis

This pipeline is independent of the COF mask pipeline (Analyses 1 and 2). It operates directly on ETC track files and HEALPix model output. Data sources include ERA5+IMERG (observations) and the same six GSRM models.

```
ETC track files  COF overlap parquet    ETC track files  HEALPix catalog
       │                │                      │                │
       └────────┬───────┘                      └────────┬───────┘
                │    (Steps 1 & 2 are independent)      │
                ▼                                       ▼
[Step 1] ETC + COF Data Preparation
  └─ Script:  scripts/combine_etc_cof_data.py
  └─ Input:   ETC track text files (/hackathon/etc_tracks/)
              COF overlap tracking parquet (/hackathon/cof_masks/stats/)
  └─ Output:  /hackathon/etc_tracks/{source}_etc_cof_data.parquet
              (ETC track points merged with COF overlap flags,
               AR/MCS partner track IDs, and ETC center coordinates)
         |
         v  (Step 2 can run in parallel with Step 1)
[Step 2] 2D Variable Extraction
  └─ Script:  extract_environments/extract_etc_2d_vars.py
  └─ Method:  Batched time-slice loading; remap HEALPix → 81×81 lat-lon grid
              centered at each ETC track point (0.25° resolution, ±10° radius)
  └─ Submit:  extract_environments/submit_etc_extraction_jobs.py
              (Python script automating Slurm job array submission per model;
               one task per variable; serial shared-queue jobs, ~15 min each)
  └─ Input:   HEALPix catalog variables + ETC track file
  └─ Output:  /hackathon/etc_data/{source}/single_vars/
              etc_2d_{variable}_{suffix}.zarr  (one zarr per variable)
              Variables: pr, huss, tas, uas, vas, psl, ua850, va850,
                         rh850, uivt, vivt, zg500, and optional 3D fields
         |
         v  (requires Steps 1 & 2)
[Step 3] Combine Individual Zarr Files
  └─ Script:  extract_environments/combine_etc_2d_vars.py
  └─ Input:   Step 2 single-variable zarr files
              Step 1 ETC COF parquet (/hackathon/etc_tracks/)
  └─ Output:  /hackathon/etc_data/{source}/
              etc_2d_combined_{suffix}.zarr
              (multi-variable zarr with unit standardization applied
               and COF overlap_flag / storm metadata added)
         |
         v
[Step 4] ETC Spatial Statistics
  └─ Script:  scripts/calc_etc_spatial_stats.py
  └─ Input:   /hackathon/etc_data/{source}/etc_2d_combined_all_all.zarr
  └─ Output:  /hackathon/etc_data/stats/
              etc_spatial_stats_{source}.nc
              (per-ETC-point spatial statistics: mean/min/max within 10°
               circular radius; COF mask fractional area at 10° and 15°;
               feature-specific precipitation statistics;
               overlap_flag coordinate for COF stratification)
  └─ Batch:   scripts/run_calc_etc_spatial_stats_all.sh
              (handles Steps 1, 3 & 4 for all 6 sources; Step 2 must be
               submitted separately via submit_etc_extraction_jobs.py:
               era5, scream, icon_d3hp003, um_glm_n2560_RAL3p3,
               nicam_gl11, casesm2_10km_nocumulus)
         |
         v
[Step 5] Visualization and Analysis
  └─ Notebooks (composite figures — input: Step 3 zarr):
     notebooks/plot_etc_composites.ipynb
       → 2D ETC composite fields for one data source
     notebooks/plot_etc_composites_diff.ipynb
       → Composite differences between a model and ERA5+IMERG
     notebooks/plot_etc_composites_multipanel.ipynb
       → Multi-panel composites for all data sources (paper figure)
  └─ Notebooks (spatial statistics — input: Step 4 NetCDF):
     notebooks/plot_etc_spatialmean_stats_multisource.ipynb
       → ETC spatial mean statistics for all sources (paper figure)
     notebooks/plot_etc_spatialmean_stats_1source.ipynb
       → ETC spatial mean statistics for one source (prototype)
```

> **Performance note:** The batched time-slice extraction in `extract_etc_2d_vars.py` achieves a ~16× speedup over naïve per-storm data loading by grouping storms by timestamp and loading only the union of required HEALPix cells per time step. See [README_BATCHED_EXTRACTION.md](../../extract_environments/README_BATCHED_EXTRACTION.md) for full details.

---

## Analysis 4 — MCS COF Track Statistics

This pipeline uses the COF masks zarr (Step 3 of the shared upstream pipeline) together with MCS track statistics NetCDF files produced by PyFLEXTRKR. It extracts per-track COF overlap flags and partner feature IDs, then merges all sources into a single analysis-ready parquet file.

```
[Step 3 output] /hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr
  +  MCS track statistics NetCDF  (/hackathon/mcs/{source}/stats/)
         |
         v
[Step 1] Extract MCS COF Flags & Partner Track Numbers
  └─ Script:  scripts/extract_mcs_cof_tracks.py  (full version)
              scripts/extract_mcs_cof_flags.py   (prototype: flags only, no partner IDs)
  └─ Method:  Scans each 6-hourly COF zarr time step; records which
              (MCS track ID, timestamp) pairs appear in each overlap
              category; maps 6-hourly COF flags back onto the 1-hourly
              MCS track-stats time dimension
  └─ Input:   /hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr
              MCS track stats NetCDF  (--trackstats; optional)
  └─ Output:  /hackathon/cof_masks/stats/
              {source}_mcs_cof_tracks_2d.nc
                (2D arrays shaped (tracks, times) — cof_flag_mcs_ar,
                 cof_flag_mcs_etc, cof_flag_mcs_ar_etc,
                 cof_flag_isolated; partner IDs ar_tracknum,
                 etc_tracknum, ar_tracknum_3way, etc_tracknum_3way)
              {source}_mcs_cof_flags.parquet
                (long-format per-track COF overlap summary)
              {source}_mcs_trackstats_cof.parquet
                (flat tidy DataFrame of MCS track stats with COF flags;
                 only written when --trackstats is provided)
         |
         v
[Step 2] Combine Multi-Source Track Statistics + COF Overlap Data
  └─ Script:  scripts/combine_mcs_cof_trackstats_multisource.py
  └─ Input:   MCS track stats NetCDF  (/hackathon/mcs/{source}/stats/)
              COF 2D tracks NetCDF    (/hackathon/cof_masks/stats/)
              Sources: IMERGv7, scream, icon_d3hp003,
                       um_glm_n2560_RAL3p3, nicam_gl11,
                       casesm2_10km_nocumulus
  └─ Steps:   (a) Load & merge track-stats and COF 2D files per source
                  on the shared (tracks, times) grid
              (b) Optional tropical filter: remove tracks whose
                  lifetime-median |meanlat| < 20° (configurable)
              (c) Convert to tidy Pandas DataFrame
                  (one row per valid MCS time step)
              (d) Concatenate all sources
  └─ Output:  /hackathon/cof_masks/stats/
              mcs_cof_trackstats_allsources.parquet
              (snappy-compressed; all sources; per-timestep MCS stats
               with COF overlap flags)
         |
         v
[Step 3] Visualization and Analysis
  └─ Notebook: notebooks/plot_mcs_cof_trackstats_multisource.ipynb
  └─ Input:  Step 2 parquet
  └─ Output: Figures comparing MCS track statistics by COF type
             across all data sources (paper figure)
```

---

## Full Pipeline Diagram

```
── COF ANALYSES (1 & 2) ──────────────────────────────────────────────────────

Hourly MCS pixel masks
AR / TC / ETC NetCDF files          Precipitation (catalog or zarr)
          │                                       │
          ▼                                       │
   [Step 1] make_mcs_swath_masks.py               │
          │                                       │
          ▼                                       │
   [Step 2] combine_tracking_masks.py             │
          │                                       │
          ▼                                       │
   [Step 3] make_cooccurrence_masks.py            │
      (cofmasks zarr)                             │
          │                    ┌──────────────────┤
          │                    │                  │
          ▼                    ▼                  ▼
   ┌──────────────────┐  ┌─────────────────┐ ┌───────────────────────┐
   │    ANALYSIS 1    │  │   ANALYSIS 2a   │ │     ANALYSIS 2b       │
   │                  │  │                 │ │                       │
   │ calc_monthly_    │  │ calc_extreme_   │ │ calc_stormtype_       │
   │ rainmap_by_cof   │  │ precip_thresh.  │ │ extreme_precip_       │
   │                  │  │                 │ │ spatial               │
   │ monthly_rainmap  │  │ precip_percen-  │ │                       │
   │ _cof_hp8_v1.nc   │  │ tiles_6h_hp8    │ │ stormtype_extreme_    │
   └────────┬─────────┘  │ _v1.nc          │ │ precip_{Pxx}_hp8_v1   │
            │            └────────────────-┘ └──────────┬────────────┘
            │                    └──────────────────────┘
            ▼                                │
   plot_cof_raintype_                        ▼
   rank_map.ipynb              plot_cof_extreme_raintype_
                               rank_map.ipynb

── ETC COMPOSITE ANALYSIS (3) ────────────────────────────────────────────────

ETC tracks  COF parquet    ETC tracks    HEALPix catalog
    │            │               │             │
    └─────┬──────┘               └──────┬──────┘
          │   (Steps 1 & 2 independent) │
          ▼                             ▼
[Step 1] combine_etc_cof_data.py  [Step 2] extract_etc_2d_vars.py
 → {src}_etc_cof_data.parquet      (Slurm array: one job per variable)
          │                             │
          └──────────────┬──────────────┘
                         │
                         ▼
              [Step 3] combine_etc_2d_vars.py
                         │
               ┌─────────┴────────┐
               ▼                  ▼
     ┌──────────────────┐  ┌──────────────────┐
     │ etc_2d_combined  │  │  [Step 4]        │
     │ _{suffix}.zarr   │  │  calc_etc_       │
     │  (composites)    │  │  spatial_stats   │
     └────────┬─────────┘  │  etc_spatial_    │
              │            │  stats_{src}.nc  │
              ▼            └────────┬─────────┘
     plot_etc_composites            │
     plot_etc_composites_diff       ▼
     plot_etc_composites_   plot_etc_spatialmean_
     multipanel.ipynb       stats_*.ipynb

── MCS COF TRACK STATISTICS (4) ──────────────────────────────────────────────

[Step 3 COF masks zarr]  +  MCS track stats NetCDF
          │                           │
          └─────────────┬─────────────┘
                        │
                        ▼
           [Step 1] extract_mcs_cof_tracks.py
                        │
                        ▼
           [Step 2] combine_mcs_cof_trackstats_multisource.py
                (all 6 sources → single parquet)
                        │
                        ▼
           plot_mcs_cof_trackstats_multisource.ipynb
```

For script-level documentation, wrapper scripts, and input/output tables, see the [documentation index](../index.md).
