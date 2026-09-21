# Running Analysis 3 (ETC Composites) in One Go (`run_etc_pipeline.py`)

**Reference scripts:**
- `scripts/run_etc_pipeline.py` — runs every step of Analysis 3 for any set of sources, in dependency order, on one node (reuses the engine of `run_cof_pipeline.py`)
- `config/config_etc_pipeline.yaml` — the sources, CPU slots and estimated durations
- `extract_environments/extract_etc_2d_vars.py` — extracts 2D fields around every ETC point (`--pr_source`, `--max_missing_fraction`)
- `extract_environments/link_etc_env_stores.py` — links the unchanged environment stores into a run
- `extract_environments/combine_etc_2d_vars.py`, `scripts/combine_etc_cof_data.py`, `scripts/create_etc_composites.py`, `scripts/calc_etc_spatial_stats.py` — the other steps
- `slurm/run_interactive_etc_pipeline.sh` — launcher for an interactive node

**Author:** Zhe Feng | zhe.feng@pnnl.gov

---

## What Analysis 3 does, and what the runner replaces

Analysis 3 composites 2D fields (161 x 161 points at 0.25°, ±20° around the storm) around every point of every extratropical-cyclone (ETC) track and stratifies them by the COF overlap flag of the ETC
(0 isolated, 1 MCS only, 2 AR only, 3 MCS and AR). Before the runner, Step 1 (`combine_etc_cof_data.py`) and Steps 3-4 were run by hand or by `run_calc_etc_spatial_stats_all.sh` (a source list inside the file), and Step 2
(the extraction, 27-33 variables per source) was submitted separately, one Slurm array per variable group, by `submit_etc_extraction_jobs.py` (which still works and shares its settings with the runner).

## The dependency graph

For each source (era5 = the observations, scream, icon_d3hp003, nicam_gl11, um_glm_n2560_RAL3p3, casesm2_10km_nocumulus):

```
etc_cof ---------------------------------\
pr ---------------------------------------\
mask_<variable> x 7 (the COF masks) --------+--> combine --> composites
link_env (unchanged environment stores) ---/             \-> stats
```

| Step | Script | Output (under `--data-root`) |
|------|--------|-------------------------------|
| `etc_cof` | `combine_etc_cof_data.py` | `etc_tracks/{src}_etc_cof_data.parquet` (track points with the Step 3 overlap flag and the partner track IDs) |
| `pr` | `extract_etc_2d_vars.py` | `etc_data/{src}/single_vars/etc_2d_pr_all_all.zarr` |
| `mask_<variable>` | `extract_etc_2d_vars.py --cof_mask` | `.../etc_2d_<variable>_all_all.zarr`, one task for each of `mcs_ar_etc_overlap_mask`, `ar_mcs_etc_overlap_mask`, `etc_mcs_ar_overlap_mask`, `etc_ar_overlap_mask`, `etc_mcs_overlap_mask`, `mcs_etc_overlap_mask`, `ar_etc_overlap_mask` |
| `link_env` | `link_etc_env_stores.py` | symbolic links to the environment stores (`huss`, `tas`, `psl`, winds, humidity, geopotential, ...) |
| `combine` | `combine_etc_2d_vars.py` | `etc_data/{src}/etc_2d_combined_all_all.zarr` |
| `composites` | `create_etc_composites.py` | `etc_data/stats/{src}/etc_2d_composite_{nh,sh}_{all,isolated,mcs_only,ar_only,3way}.nc` |
| `stats` | `calc_etc_spatial_stats.py` | `etc_data/stats/etc_spatial_stats_{src}.nc` |

A step starts when the steps it needs have exited with status 0; a failed step skips the later steps of its own source only. The scheduler, the memory gate, the markers, `--resume` and the protection of
outputs that the runner did not create are those of `run_cof_pipeline.py` (see [run_cof_pipeline.md](run_cof_pipeline.md)); the COF products and the track files are inputs and are not under `--data-root`.

## Running it

```bash
cd /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature
bash slurm/run_interactive_etc_pipeline.sh --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/etc_round1 --dry-run   # graph, commands, estimate
bash slurm/run_interactive_etc_pipeline.sh --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/etc_round1             # the run, in the background
```

| Option | Effect |
|--------|--------|
| `--data-root DIR` | required: receives `etc_tracks/`, `etc_data/`, the markers and logs. Use a test area first; pointing it at the production tree stops at once unless `--force` (the existing outputs have no markers) |
| `--cof-root DIR` | where the Step 3 products (`cof_masks/`) are read from (default: the production tree) |
| `--env-from DIR` | folder with `<source>/single_vars/` environment stores to link (default `/pscratch/sd/w/wcmca1/hackathon/etc_data/`) |
| `--sources um scream` | only these sources (names or aliases such as `imerg`, `obs`, `icon`) |
| `--steps combine composites stats`, `--from combine` | only these steps (`mask` stands for the seven mask tasks) |
| `--resume`, `--force`, `--dry-run`, `--preflight-only`, `--skip-preflight`, `--max-slots`, `--stagger-sec`, `--min-free-gb` | as in `run_cof_pipeline.py` |
| `--step-args STEP "ARGS"` | extra arguments for one step, e.g. `--step-args pr "--start_date 2020-03-01 --end_date 2020-03-08"` for a short test |

## Precipitation and time matching

`pr` of the models is Step 1's `tot_pr`, read from `cof_masks/{source}_cofmasks_hp8_v1.zarr`: the mean of the hourly precipitation over the 6-hour window that starts at the label time T, the window
that the COF masks, the COF flags and Analyses 1, 2 and 4 use. The observations keep IMERG 6-hourly (`tot_pr` of the IMERG source is 0 poleward of 60°, and the ETC points reach 70°). The separate 6-hourly files and the
catalog 6-hour mean that Analysis 3 used before are still available with `--pr_source legacy`; they differ from `tot_pr` as follows (correlation with `tot_pr` at label offsets of -6 h, 0 and +6 h, three 4-day blocks per store, |latitude| < 55°):

| Source | `pr` before | -6 h | 0 h | +6 h | Reading |
|--------|-------------|------|-----|------|---------|
| UM | local 6-hourly file | 0.28-0.34 | **1.000** | 0.28-0.34 | same window |
| CASESM2 | local 6-hourly file | 0.38-0.43 | **1.000** | 0.38-0.42 | same window |
| IMERG (ERA5 source) | `IMERG_V7_6H_zoom8` | 0.38-0.49 | **0.998-1.000** | 0.39-0.48 | same window (kept) |
| NICAM | local 6-hourly file (shifted) | 0.21-0.30 | **0.94-0.95** | 0.39-0.47 | window matches; 3-hourly means against hourly values |
| SCREAM | `scream_pr6h_z8` (liquid only) | 0.61-0.63 | 0.72-0.74 | 0.25-0.27 | window centred on T: straddles two windows |
| ICON | catalog `PT6H` mean | **0.967-0.969** | 0.40-0.48 | 0.17-0.28 | end-stamped: the value at T is the mean of [T-6 h, T), a full window early |

The track files are exactly 6-hourly (00, 06, 12, 18 UTC), so the matching itself is exact. **A track time that has no frame in the source is NaN, never the nearest frame.** `sel(time, method="nearest")` never fails: it returns the last (first)
frame for a time after (before) the record, which the old extraction would do against the current COF stores, whose periods are shorter than the track periods. On the extraction's point list (|latitude| ≤ 70°) that concerns
ICON 3299 of 22594 points (14.6%: 2020-01-01 and 2021-01/02; the COF period is one year), UM 749 of 23287 (3.2%: 20-31 January 2020), SCREAM 121 of 22301 (0.5%: seven consecutive windows, 20 April 06 UTC to 21 April 18 UTC 2020, are missing from the COF store); NICAM, CASESM2 and the observations none.
In the March data these points have `overlap_flag` NaN, `pr` from the catalog, and masks that are NaN for ICON and partly finite for UM (36%) and SCREAM (12%); now their `pr` and masks are NaN and the `overlap_flag` stays NaN, so they stay identifiable and
drop out of the composites of a single flag and out of the statistics that are filtered by `overlap_flag`. They are still in the `all` composites (see "Latitude rule and radii" below) and in the environment variables, which do not depend on the COF period. Each extraction records `n_points`, `n_points_time_missing` and `n_points_failed` in the store attributes and fails before writing when more points are
missing than `--max_missing_fraction` (the registry passes the expected share per source: SCREAM 0.008, ICON 0.15, UM 0.035, others 0) or when a slice could not be loaded.

`combine_etc_2d_vars.py` applies the source-specific factor to a native flux (SCREAM ×3.6·10⁶, the others ×3600) but no longer scales a `pr` that is already in mm h⁻¹; it fails when the stores do not share one storm-point list (it used to print a warning and
merge by position), when a variable cannot be loaded (`--allow-missing-variables` to skip it), and when the mean of `pr` is outside 0.01-1 mm/h (about 0.1 mm/h over an ETC neighbourhood).

## The environment variables are not extracted again

The environment variables do not depend on the COF products or on the precipitation. Their stores from the March extraction (`etc_data/{source}/single_vars/`, restored from CFS on 2026-09-20 after the scratch purge emptied them: 25/19/19/19/19/20 environment
stores for SCREAM/ICON/NICAM/UM/CASESM2/ERA5) are linked into the run by `link_env` after each store has been checked (complete chunk by chunk, with metadata, and with the number of storm points of the current track file). UM and CASESM2, whose catalogs are online only, therefore
never touch the network in this pipeline. To re-extract them, use `submit_etc_extraction_jobs.py` as before.

## Preflight

Before anything starts the runner checks (directory listings only): zarr 2 in the environment; for each source the track file and the Step 3 overlap parquet, the COF store (chunk by chunk; the IMERG 6-hourly store for the observations' `pr`), and the environment stores
(complete, metadata present, the right number of storm points). `scripts/check_zarr_store.py` reports a store without array metadata as not complete, which is what a purged store looks like.

## Resources and time

Measured in the full run of 2026-09-20 (all six sources, one node with 256 CPUs and 503 GB; minutes and peak resident memory of one step):

| Step | Models | ERA5 |
|------|--------|------|
| `pr` or one of the seven `mask_*` extractions | 4.4-6.5 min, 2.4-3.1 GB | 15.8-16.7 min, 7.7-10.3 GB |
| `etc_cof`, `link_env` | seconds | seconds |
| `combine` | 1.9-2.7 min, 50-75 GB | 5.9 min, 195 GB |
| `composites` | 0.8-1.3 min, 22-47 GB | 2.2 min, 61 GB |
| `stats` | 0.3-0.5 min, at most 25 GB | 0.8 min, 19 GB |

The extraction reads whole time chunks of the source once, not once per storm time. The 60 extraction steps of the six sources ran side by side in 19.3 min (at most 44 CPU slots in use); combine, composites and stats took another 9.9 min, so the whole core tier takes about 30 min
(`--dry-run` estimates 27 min). The `combine` step holds every variable in memory (about 2.3 GB per variable for 22,000 storm points, ERA5 three times that), so it takes many CPU slots in the registry and only two run at a time.

## Latitude rule and radii of the composites and statistics

- **Polar ETCs are left out.** The ETC tracks are global (the extraction keeps |latitude| <= 70°), but the COF products exist only equatorward of 60°: the MCS masks, and with them the overlap flags and the MCS/AR masks, are not available poleward of that, so an ETC there is
  "isolated" or "AR only" by construction (before 2026-09-21 nothing removed them: 20-26% of the northern and 35-44% of the southern points the composites select lie poleward of 60°, and 42-59% of the "isolated" ETC points of the statistics were polar). An ETC is used only if at least
  `--min-lat-coverage` of its box (161 rows, ±20°) lies within `|latitude| <= --lat-limit` (60°): `create_etc_composites.py` default 0.8 (centre within 48°, keeps 32-43% of the northern and 26-33% of the southern points), `calc_etc_spatial_stats.py` default 0.5 (centre within 60°, the centroid
  rule; keeps 67-72% of the points; the excluded points are dropped from the output, so a track is truncated where it leaves the domain). 1.0 means the centre within 40°, 0 switches the rule off. The same value in both scripts gives the same sample (`--step-args composites "--min-lat-coverage 0.5"`
  changes one run of the runner). The output files record the rule and the number of points (`lat_limit`, `min_lat_coverage`, `n_points`, `n_points_hemisphere_before_lat_rule` or `n_points_before_lat_rule`); look at the number of points of the small categories (northern `ar_only` composites have 28-100 points for
  SCREAM, ICON, ERA5 and NICAM at 0.8, 5-13 tracks of the type ETC+AR in the northern statistics for ICON, SCREAM and NICAM).
- **The radii of the statistics are in degrees** (fixed 2026-09-21): `x` and `y` are offsets in grid points (0.25° apart) and are converted with the global attributes `lon_res` and `lat_res` before they are compared with `--basic-radius`, `--pr-radius` and `--mask-radii` (10, 10 and 10 15 by default,
  40 and 60 grid points). The circle is a circle in degrees (like the ones drawn in the composite notebooks), so it is narrower in km towards the poles. Statistics made before that date used the radii as grid points (2.5° and 3.75°); the attribute `radius_units = degrees` marks the new files.
- **The `all` composites include the points that have no frame of the source** (with the shares of the composites' own hemisphere selection, before the latitude rule: ICON 16.7% of the northern and 12.6% of the southern points, UM 3.7% and 2.8%, SCREAM 0.6% and 0.5%, the others none): their `pr` is NaN and skipped, but
  the mask frequencies (`ar_freq`, `mcs_freq`, `etc_freq`) and the exclusive precipitation (`pr_ar`, `pr_mcs`, `pr_etc`) count them as zeros, so these variables, and the fractions `prfrac_mcs` and `prfrac_ar` made from them, are low by that share in the `all` files (ICON NH: MCS frequency maximum 0.126 as stored against 0.151
  from the frame points only). The files of a single flag and the statistics filtered by `overlap_flag` are not affected. Not changed; see the audit findings.
- The masks of the combined store are 0 where there is no object (the March stores had NaN); the composites and statistics use `mask > 0`.

Details and the comparison with the March data are in `docs/AUDIT_FINDINGS.md` (Analysis 3 section).

## Troubleshooting

- A step that failed: read `<data-root>/pipeline_logs/<run id>/<source>_<step>.log`, fix the cause, run again with `--resume`.
- `... have no frame in the source ... more than --max_missing_fraction`: the COF store (or the source) is shorter than the track file by more than expected; check the periods, or pass a larger `--step-args pr "--max_missing_fraction X"`.
- `stores are not aligned`: a single-variable store was extracted from another track file or with other filters than the rest; extract it again.
- `no array metadata found`: the store is empty (scratch purge); restore it from CFS.
- To run the old wrapper for Steps 1, 3 and 4 of one source: `scripts/run_calc_etc_spatial_stats_all.sh SOURCE` (paths hard-coded to production).
