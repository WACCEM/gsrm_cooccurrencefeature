# Documentation Index

Start with the [README](../README.md) for the project overview, then use the [analysis pipeline overview](pipelines/overview.md) for the end-to-end workflow, and finally open the procedure docs for implementation details of individual processing steps.

## Pipeline and Procedure Docs

| Topic | Documentation |
|-------|---------------|
| End-to-end analysis pipelines | [pipelines/overview.md](pipelines/overview.md) |
| MCS swath masks and cloud type classification | [procedures/mcs_swath_cloud_type.md](procedures/mcs_swath_cloud_type.md) |
| Combined tracking mask creation | [procedures/combine_tracking_masks.md](procedures/combine_tracking_masks.md) |
| Co-occurrence feature identification | [procedures/cof_identification.md](procedures/cof_identification.md) |
| Monthly precipitation by COF type | [procedures/monthly_precip_by_cof.md](procedures/monthly_precip_by_cof.md) |
| Extreme precipitation by storm type | [procedures/extreme_precip_by_stormtype.md](procedures/extreme_precip_by_stormtype.md) |
| Batched ETC environment extraction | [../extract_environments/README_BATCHED_EXTRACTION.md](../extract_environments/README_BATCHED_EXTRACTION.md) |
| ETC variable scaling | [../extract_environments/README_VARIABLE_SCALING.md](../extract_environments/README_VARIABLE_SCALING.md) |
| ETC extraction optimization notes | [../extract_environments/OPTIMIZATION_NOTES.md](../extract_environments/OPTIMIZATION_NOTES.md) |

## Script and Notebook Documentation

| Script / Notebook | Documentation |
|-------------------|---------------|
| `make_mcs_swath_masks.py` | [procedures/mcs_swath_cloud_type.md](procedures/mcs_swath_cloud_type.md) |
| `combine_tracking_masks.py` | [procedures/combine_tracking_masks.md](procedures/combine_tracking_masks.md) |
| `combine_era5_imerg_tracking_masks.py` | Needs standalone documentation |
| `make_cooccurrence_masks.py` | [procedures/cof_identification.md](procedures/cof_identification.md) |
| `calc_monthly_rainmap_by_cof.py` | [procedures/monthly_precip_by_cof.md](procedures/monthly_precip_by_cof.md) |
| `calc_extreme_precip_thresholds.py` | [procedures/extreme_precip_by_stormtype.md](procedures/extreme_precip_by_stormtype.md) |
| `calc_stormtype_extreme_precip_spatial.py` | [procedures/extreme_precip_by_stormtype.md](procedures/extreme_precip_by_stormtype.md) |
| `extract_etc_2d_vars.py` | [../extract_environments/README_BATCHED_EXTRACTION.md](../extract_environments/README_BATCHED_EXTRACTION.md) |
| `combine_etc_cof_data.py` | Documented in [pipeline overview](pipelines/overview.md) |
| `combine_etc_2d_vars.py` | Documented in [pipeline overview](pipelines/overview.md) |
| `calc_etc_spatial_stats.py` | Documented in [pipeline overview](pipelines/overview.md) |
| `extract_mcs_cof_tracks.py` | Documented in [pipeline overview](pipelines/overview.md) |
| `extract_mcs_cof_flags.py` | Prototype documented in [pipeline overview](pipelines/overview.md) |
| `combine_mcs_cof_trackstats_multisource.py` | Documented in [pipeline overview](pipelines/overview.md) |

## Wrapper Shell Scripts

These scripts automate running one or more Python processing scripts across model sources in a single invocation. They are the primary entry points for batch execution.

| Script | Analysis | Python Script(s) Called | Sources | Notes |
|--------|----------|-------------------------|---------|-------|
| `slurm/slurm_make_mcs_swath_masks.sh` | Shared Step 1 | `make_mcs_swath_masks.py` | 6 sources via task file | SLURM job array (`--array=1-6`); reads commands from `tasks_make_mcs_swath_masks.txt` |
| `run_combine_tracking_masks_all.sh` | Shared Step 2 | `combine_tracking_masks.py` | scream, icon, nicam, um, casesm2 | Optionally pass a single source as argument |
| `slurm/slurm_make_cooccurrence_masks.sh` | Shared Step 3 | `make_cooccurrence_masks.py` | 6 sources via task file | SLURM job array (`--array=1-6`); reads commands from `tasks_make_cooccurrence_masks_all.txt` |
| `slurm/slurm_calc_monthly_rainmap_by_cof.sh` | A1 Step 4 | `calc_monthly_rainmap_by_cof.py` | 6 sources via task file | SLURM job array (`--array=1-6`); reads commands from `tasks_calc_monthly_rainmap_by_cof.txt` |
| `run_all_extreme_precip_thresholds.sh` | A2 Step 4a | `calc_extreme_precip_thresholds.py` | all sources from config | Uses `config_sources.yaml` |
| `run_all_extreme_precip_thresholds_1h.sh` | A2 Step 4a | `calc_extreme_precip_thresholds_1h.py` | all sources from config | For 1-hourly input data; uses `config_sources_1h.yaml` |
| `run_stormtype_extreme_precip.sh` | A2 Step 4b | `calc_stormtype_extreme_precip_spatial.py` | scream, icon, nicam, um, casesm2, IR_IMERG | Optionally pass a single source as argument |
| `extract_environments/submit_etc_extraction_jobs.py` | A3 Step 2 | `extract_etc_2d_vars.py` | all sources | SLURM job array for each group of variables |
| `run_calc_etc_spatial_stats_all.sh` | A3 Steps 1, 3, 4 | `combine_etc_cof_data.py` -> `combine_etc_2d_vars.py` -> `create_etc_composites.py` -> `calc_etc_spatial_stats.py` | era5, scream, icon, um, nicam, optionally casesm2 | Step 2 extraction must run separately via `submit_etc_extraction_jobs.py` before this script |
| `run_extract_mcs_cof_tracks_all.sh` | A4 Step 1 | `extract_mcs_cof_tracks.py` | all sources | Produces `{source}_mcs_cof_tracks_2d.nc` |
| `combine_mcs_cof_trackstats_multisource.py` | A4 Step 2 | `combine_mcs_cof_trackstats_multisource.py` | all sources | Produces `mcs_cof_trackstats_allsources.parquet` |

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
| 1 (A3) | `combine_etc_cof_data.py` | ETC track text files + COF overlap tracking parquet | `/hackathon/etc_tracks/{source}_etc_cof_data.parquet` |
| 2 (A3) | `extract_etc_2d_vars.py` | HEALPix catalog data + ETC track file | `/hackathon/etc_data/{source}/single_vars/etc_2d_{var}_{suffix}.zarr` |
| 3 (A3) | `combine_etc_2d_vars.py` | Step 2 (A3) zarr + Step 1 (A3) parquet | `/hackathon/etc_data/{source}/etc_2d_combined_{suffix}.zarr` |
| 4 (A3) | `calc_etc_spatial_stats.py` | Step 3 (A3) zarr | `/hackathon/etc_data/stats/etc_spatial_stats_{source}.nc` |
| 5 (A3) | Composite notebooks | Step 3 (A3) zarr | Figures (2D composite maps) |
| 5 (A3) | Spatial stats notebooks | Step 4 (A3) nc | Figures (ETC spatial mean statistics) |
| 1 (A4) | `extract_mcs_cof_tracks.py` | Step 3 cofmasks zarr + MCS track stats nc | `/hackathon/cof_masks/stats/{source}_mcs_cof_tracks_2d.nc`, `{source}_mcs_trackstats_cof.parquet` |
| 2 (A4) | `combine_mcs_cof_trackstats_multisource.py` | Step 1 (A4) nc + MCS track stats nc for all sources | `/hackathon/cof_masks/stats/mcs_cof_trackstats_allsources.parquet` |
| 3 (A4) | `plot_mcs_cof_trackstats_multisource.ipynb` | Step 2 (A4) parquet | Figures (MCS track statistics by COF type) |

## Assets

| Asset | Description |
|-------|-------------|
| [assets/COF_schematic.pdf](assets/COF_schematic.pdf) | COF flowchart schematic |
