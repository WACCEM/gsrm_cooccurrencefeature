# Re-running the analyses after a tracking update (TC, ETC or AR masks)

When the TempestExtremes products change upstream (2026-09-25: quasi-stationary storms removed from the ETC and TC tracks, then a corrected TC filter), not everything has to be redone. This page says what depends on what, which steps to run for which change, how to make the ETC
environment stores, and what to check. Commands use a new data root (never the production tree until the results are reviewed).

## What reads what

| Product | Reads TC masks | Reads ETC masks / tracks | Reads AR masks | Notes |
|---|---|---|---|---|
| Step 1 `make_mcs_swath_masks.py` (mcs_masks) | **yes**: MCS-TC exclusion and cloud types | no | no | `tot_pr` does not depend on any tracking |
| Step 2 `combine_tracking_masks.py` (all_masks) | yes (`tc_mask`) | yes (`etc_mask`) | yes (`ar_mask`) | reads the netCDF of `dir_te` |
| Step 3 `make_cooccurrence_masks.py` (cof_masks, `*_etc_overlap_tracking.parquet`) | yes: also zeroes cloud types, AR and ETC pixels under TC | yes | yes | |
| monthly (Analysis 1), attribution (Analysis 2) | via Step 3 | via Step 3 | via Step 3 | |
| thresholds (Analysis 2) | no | no | no | from `tot_pr`; identical when Step 1's `tot_pr` is |
| Analysis 3 `etc_cof` | via the Step 3 parquet | the ETC stitched-node file | | |
| Analysis 3 `pr` | no | the track file | no | `tot_pr` at the ETC points |
| Analysis 3 environment stores (`env_<store>`, `link_env`) | **no** | **the ETC stitched-node file only** | no | values depend on the point's time and position only |
| Analysis 3 seven `mask_*` stores, `combine`, `composites`, `stats` | via Step 3 | via Step 3 | via Step 3 | |
| Analysis 4 (`extract_mcs_cof_tracks.py`, multi-source combine) | via Step 3 | via Step 3 | via Step 3 | |

Bryce's `filter_TCs_and_ETCs.py` removes ETC storms that coincide with the **raw** TC nodes (`tc_stitched_nodes.txt`) before the quasi-stationary filter, so a change of the qs-filtered TC files does not change the ETC stitched-node file (check it: below).

## Which steps to redo

| What changed | Step 1 | Steps 2, 3, monthly, attribution | thresholds | Analysis 3 | Analysis 4 |
|---|---|---|---|---|---|
| TC masks only (ETC track file identical) | **yes** | yes | keep (or 1.4 min to redo) | `etc_cof`, `mask_*`, combine, composites, stats; `link_env` (environment stores unchanged); `pr` unchanged (6 min to redo) | yes |
| ETC and/or AR only (TC identical) | **no** (link `mcs_masks/`) | yes | keep | as above, and the environment route below if the ETC file changed | yes |
| TC and ETC | yes | yes | keep | environment route below | yes |

Measured 2026-09-25 (node with 256 CPUs; ICON): Step 1 15 min (13-15 for the models with 24 workers each, all side by side; IMERG 55 min), Step 2 0.5, Step 3 6, monthly 1.2, thresholds 1.4, attribution 6.4 (about 30 min per source); Analysis 3 without environment extraction about 12 min;
Analysis 4 under 1 min per source.

## Procedure

0. **Check the inputs** for every source (`scripts/check_tracking_inputs.py`, see its help): new ETC file vs the old one and vs the previous round's, ETC masks vs the file. It prints the environment route. With `--tc-track-file` it also checks the TC side: how each raw TC storm was treated (input | unchanged | trimmed | split | removed | output, the terms of Bryce's logs `<src>_qs_log_TCs_<criteria>.csv`) and the TC masks against the TC stitched-node file, in the first, the last and the first changed month. Do not go on with a `MISMATCH`.
   **Verify the content, not the file name.** Several versions of a track file coexist (`tc_stitched_nodes.qs_filter_r30_d48_wrong.txt`, `r30_d48`, `r15_d48`, `r15_d96`) and the name of the file in a list or a script is not proof of which one the masks were built from. On 2026-09-25 the list of input files named the `r30_d48` TC files, while the TC masks in use were built from `r15_d96` (Bryce's final TC criteria: 1.5 GCD and 96 h, plus the TC minima); the check told them apart at once (pointed at `r30_d48` it reports `MISMATCH`, exit status 1).
1. **New data root**, e.g. `/pscratch/sd/w/wcmca1/hackathon/tmp/<name>/`. For an ETC/AR-only update link the unchanged `mcs_masks/` into it.
2. **COF pipeline:** `bash slurm/run_interactive_cof_pipeline.sh --data-root ROOT --sources ... --analysis both` (all steps; steps that are present are not run again with `--resume`).
3. **Analysis 3** per route (the environment stores depend on the ETC file only):

   | Route | When | Command |
   |---|---|---|
   | `link` | new ETC file identical to the previous round's | `bash slurm/run_interactive_etc_pipeline.sh --data-root ROOT --cof-root ROOT --tracks-dir PREV/etc_tracks --env-from PREV/etc_data/ --sources ...` |
   | `reuse` | every point of the new file is in the old file (removal, trimming, splitting only) | `... --tracks-dir ROOT/etc_tracks --reuse-env --sources ...` (copy the new ETC file into `ROOT/etc_tracks/` under the registry's file name) |
   | `extract` | new points exist | `... --tracks-dir ROOT/etc_tracks --extract-env --sources ...` (online catalogs: 4 readers at most, hours) |

4. **Analysis 4:** `extract_mcs_cof_tracks.py --cof-dir ROOT/cof_masks --output-dir ROOT/cof_masks/stats ...` per source, then `combine_mcs_cof_trackstats_multisource.py --cof-dir ROOT/cof_masks/stats/ --output-dir ROOT/cof_masks/stats/`.
5. **Checks:** exact comparison of the three mask stores and the monthly file against the previous round (`tot_pr` and `ar_mask` must be identical when only TC changed; ETC and TC footprints change), thresholds identical, the number of storm points in the A3 stores equal to the track file's, then the notebooks.

## Two rules that keep the environment stores right

- The rows of an environment store are matched to a new track file by (time, lon, lat), never by storm ID (IDs are running counters: removing or splitting one storm renumbers all later ones). `--reuse-env` checks that the old store is the extraction of the old file row by row
  and that every new point is found, and fails otherwise. Checked 2026-09-25: for all six sources every point of the new file is in the old file, in the same order (ICON: 19 of 19 rebuilt stores identical to a fresh extraction; SCREAM 16 of 16).
- Points without a frame in the catalog are NaN and recorded, never the nearest frame, and never a crash below the sanity limit ([run_etc_pipeline.md](run_etc_pipeline.md)).

## ERA5 masks and the IMERG (observation) chain

The observations take their AR, TC and ETC masks from ERA5 (`ERA5_AR_TC_ETC_hp8_v1.zarr`, made by `scripts/remap_era5_masks_healpix.py` from Bryce's `hackathon_pre/era5_tracking_etc_nocoldcoreonly/`), and only 2019-2021 is used (the `*_full_*` tracking of 1979-2021 is not).

- The remap reads `AR_tracks_era5_*.nc`, `TC_test_tracks_era5_*.nc` and `ETC_test_tracks_era5_*.nc` from **one folder** and keeps only the times that all three have; IMERG's Step 2 intersects times in the same way. A folder with 20 months of ETC masks makes a 20-month store without a message.
  Before Step 2 check that the three inputs cover the same 36 months (`check_tracking_inputs.py`) and that the new store has 4384 six-hourly steps (2019-01-01 to 2021-12-31T18).
- Verify the content here too. On 2026-09-26 the folder held 20 monthly ETC masks (201908-202103) and 1096 daily TC masks built from the *short* TC file (2019-08-01 to 2021-03-26, 143 storms), which are empty outside that window (2019-02-01: 0 TC cells, 10,626 in the March masks).
  r15_d96 changes none of the 244 TC storms of 2019-2021 (the March masks `TC_test_tracks_old_SN_*` equal them), so the ERA5 TC masks of 2019-2021 are the March ones.
- Only the ETC masks of ERA5 change with the quasi-stationary filter. IMERG Step 1 reads the TC footprint only, so with an unchanged footprint (compare the daily masks with `old_SN` on every day) it is not re-run: link the production `mcs_masks/IMERGv7_mcs_masks_hp8.zarr`
  (for NICAM and UM, whose TC masks did not change either, Step 1 came out bit-identical). The thresholds (`v1_2014_2024`) do not read the masks and are linked. Steps 2, 3, monthly, attribution, Analysis 3 (route `reuse`) and Analysis 4 follow the procedure above.
