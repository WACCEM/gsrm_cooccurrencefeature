# Analysis Pipeline Overview

**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This document describes the end-to-end processing pipelines for two analyses:

1. **Total precipitation by COF type** — monthly mean precipitation maps attributed to each co-occurrence feature category
2. **Extreme precipitation by COF type** — spatial distribution of extreme precipitation events attributed to each co-occurrence feature category

Both analyses share a common upstream pipeline (Steps 1–3) that produces the COF mask dataset. They then diverge into separate downstream processing and visualization steps.

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

## Full Pipeline Diagram

```
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
   │ _cof_hp8_v1.nc   │  │ tiles_6h_hp8   │ │ stormtype_extreme_    │
   └────────┬─────────┘  │ _v1.nc          │ │ precip_{Pxx}_hp8_v1  │
            │            └────────────────-┘ └──────────┬────────────┘
            │                    └──────────────────────┘
            ▼                                │
   plot_cof_raintype_                        ▼
   rank_map.ipynb              plot_cof_extreme_raintype_
                               rank_map.ipynb
```

---

## Documentation Index

| Script / Notebook | Documentation |
|-------------------|--------------|
| `make_mcs_swath_masks.py` | [mcs_swath_cloud_type_procedure.md](mcs_swath_cloud_type_procedure.md) |
| `combine_tracking_masks.py` | [combine_tracking_masks_procedure.md](combine_tracking_masks_procedure.md) |
| `make_cooccurrence_masks.py` | [cof_identification_procedure.md](cof_identification_procedure.md) |
| `calc_monthly_rainmap_by_cof.py` | [calc_monthly_rainmap_by_cof_procedure.md](calc_monthly_rainmap_by_cof_procedure.md) |
| `calc_extreme_precip_thresholds.py` | [extreme_precip_by_stormtype_procedure.md](extreme_precip_by_stormtype_procedure.md) |
| `calc_stormtype_extreme_precip_spatial.py` | [extreme_precip_by_stormtype_procedure.md](extreme_precip_by_stormtype_procedure.md) |

---

## Input and Output File Summary

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `make_mcs_swath_masks.py` | Hourly MCS zarr + catalog pr/Tb | `/hackathon/mcs_masks/{source}_mcs_masks_hp8.zarr` |
| 2 | `combine_tracking_masks.py` | Step 1 zarr + AR/TC/ETC NetCDF | `/hackathon/all_masks/{source}_allmasks_hp8_v1.zarr` |
| 3 | `make_cooccurrence_masks.py` | Step 2 zarr | `/hackathon/cof_masks/{source}_cofmasks_hp8_v1.zarr` |
| 4 (A1) | `calc_monthly_rainmap_by_cof.py` | Step 3 zarr + catalog pr | `/hackathon/cof_masks/stats/monthly/{source}_monthly_rainmap_cof_hp8_v1.nc` |
| 4a (A2) | `calc_extreme_precip_thresholds.py` | Catalog pr | `/hackathon/extreme_precip/{source}_precip_percentiles_6h_hp8_v1.nc` |
| 4b (A2) | `calc_stormtype_extreme_precip_spatial.py` | Step 3 zarr + Step 4a nc + catalog pr | `/hackathon/extreme_precip/{source}_stormtype_extreme_precip_{Pxx}_hp8_v1.nc` |
| 5 (A1) | `plot_cof_raintype_rank_map.ipynb` | Step 4 (A1) nc | Figures |
| 5 (A2) | `plot_cof_extreme_raintype_rank_map.ipynb` | Step 4b (A2) nc | Figures |
