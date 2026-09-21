# COF residual audit findings (2026-09-19, node nid004162)

All diagnostics were read-only. The scripts and raw outputs are not in the repository; they are in
`/global/homes/f/feng045/.claude_tmp_verify/audit/` (`audit_*.py`, `budget_*.json`, `orphans_*.csv/json`,
`tcfrac_*.csv/json`, `log_*.txt`), and the one-month test tooling is in `.../fixtest/` and `.../thresholds/`.

## Status (updated 2026-09-19, evening)

| Finding | Resolution |
|---|---|
| 1. Precipitation-product mismatch | Step 1 writes `tot_pr`, the window mean of each model's own hourly `pr` (no frozen term); Steps 2-3, the monthly script and the extreme attribution read it (cddd292). One-month tests: residual 0.0000% in every box for SCREAM, ICON, UM, NICAM and CASESM2 |
| 2. Category double counting | The monthly script assigns each pixel-frame to one category by priority (d2a77e8). The overlap in Step 3's pair lists is not fixed (follow-up 1) |
| 3. Step 3 MCS-TC re-filter | Log-only; Step 1's pooled hourly test is the only exclusion (044d0e7) |
| 4. Zero-byte chunks | Restored from CFS and verified (section 4) |
| Found while testing: Step 1 lost frames silently | Partial windows are written, retries actually re-run, failures are reported and the exit status is 1 (60c18c1) |

Test months: February 2020 (SCREAM, ICON, UM, IMERG), March 2020 (CASESM2), January 2021 (NICAM). IMERG keeps a residual of +0.24% at
60S-60N because some windows have no IR brightness temperature (2020-02-26 06 and 12 UTC, parts of 2020-02-04 and 2020-02-26 18 UTC);
that is missing input, not a pipeline error. The old IMERG residual of -1.7% hid it under the category overlap.

A full-scale re-process of SCREAM (last section) gives a residual of 0.0000% in all 13 months and in every region of the notebook's residual metric.

**Decisions of 2026-09-20 (after the full re-process of all six sources; sections at the end):**
- The extreme-precipitation thresholds are computed from Step 1's `tot_pr` (`calc_extreme_precip_thresholds.py --input_zarr ... --input_var tot_pr`); IMERG's come from the non-IR 6-hourly store, which has data at all latitudes (commit 2d3907e).
- Step 1 removes precipitation at pixels with missing Tb, since they cannot be classified (d2bcaab). This removes the IMERG residual.
- The whole of Analyses 1 and 2 runs from one dependency-aware runner, `scripts/run_cof_pipeline.py` (378ecf4), and every script takes its data root from `COF_DATA_ROOT` (dfa60ed).
- The test-area outputs are not promoted to production yet; another round of re-processing may be needed.
- Held for after this phase (see "Follow-ups after this phase" at the end): which SCREAM 6-hourly file `extract_etc_2d_vars.py` reads, and the Step 3 pair-list overlap (option C stays).

## Budget identity (validated to ~1e-5 mm per cell against calc_monthly_rainmap_by_cof.py's own function)

Residual R = sum_{m==0}(P - C) - sum_{m>=2}(m-1)P, with m = number of the 8 plotted category flags at a
pixel-frame, P = the monthly script's 6-hourly precipitation, C = cloud-type precipitation (dc+st+nd+dz).

## 1. The large Residual / over-count is a precipitation-PRODUCT mismatch, not missing classification

Step 1 classifies from hourly precipitation; the monthly script takes the TOTAL from a separate 6-hourly product.
Rebuilding the total from the same hourly data Step 1 used makes `classified_short`, `classified_over`,
`raw_only` and `unclassified` exactly 0 in all four sources and both seasons.

| source | Step 1 hourly pr | monthly total P | effect |
|---|---|---|---|
| ICON | PT1H inst, `pr` only (total precip, no `prs` in this product) | PT6H mean: `pr` (total) + `prs` (snow part) | snow counted twice (P ~ 2x where it snows: 0.0496+0.0462 vs 0.0506 mm/h in 40-60N,60-110E, Feb); PT6H mean is end-labelled (first label 06:00), so P(T) covers (T-6h,T]: best-match offset -5 h, r=0.966 |
| SCREAM | `scream2D_hrly` pr = precip_total_surf_mass_flux | `scream_pr6h_z8.zarr` pr = precip_liq_surf_mass_flux (liquid only; `prs` build line commented out in avg_pr_healpix.py) | snow missing (0.0105 vs 0.0000 mm/h); also built from end-labelled 3-hourly means: best offset -2..-3 h (r=0.88) |
| UM | pr (rain) + prs (snow) | 6-hourly `pr` only; `varname_precip_ice` absent from config_sources.yaml; the file's `prs` is all NaN | snow missing (B/C 0.26-0.39 in snowy regions); alignment fine (r=0.997) |
| NICAM | pr | pr | consistent (r~0.93-1.0, ratio ~1.00) |

What-if (60S-60N / W.China-Russia winter box), Residual %:

| source | as run | totals rebuilt from Step 1 hourly | ICON fix (PT6H pr only, taken at T+6h) |
|---|---|---|---|
| ICON | +3.1 / +28.2 | -3.4 / -0.5 | -3.0 / -0.6 |
| SCREAM | -3.5 / -1992 | -2.8 / 0.0 | |
| UM | -2.6 / -248 | -0.2 / 0.0 | |
| NICAM | -2.1 / +0.6 | -2.5 / 0.0 | |

## 2. True category double counting (E) is real but small: 2-3.5% at 60S-60N, 4-9% in mid-latitude bands

Dominated (76-96%) by different objects overlapping but assigned to different categories:
mcs_ar+mcs_etc (one MCS in two pair lists), mcs_ar+ar_etc (one AR in two pair lists),
mcs_iso+mcs_ar_etc (an "isolated" MCS inside a 3-way category's whole AR/ETC footprint).

## 3. Step 3 MCS-TC re-filter (hypothesis confirmed)

- Removed (frame, track) pairs: ICON 265, SCREAM 189 (1 frame unreadable), UM 219 (4 frames unreadable), NICAM 243.
- Every removed track passes Step 1's pooled test (f_pool < 0.10 in 100%); Step 1's TC(T) == Step 2's TC(T) in all frames.
- Swath-basis fraction (Step 3) median 12.5-12.9%: borderline. Union footprint >= 0.10 in 72-81%; the rest from coverage-priority pixel loss.
- A per-hour test would have removed 86-90% of them (Step 1 as coded pools the window, not per hour).
- Orphaned rain: 0.09-0.24% of total precipitation at 60S-60N, ~80% of removals equatorward of 30 deg. Not the mid-latitude cause.

## 4. Data integrity: 19 zero-byte files in the hourly MCS-mask stores (Step 1 inputs) - restored

scream 5, um 4, nicam 3, casesm2 2, IMERGv7 5; all mtime 2026-09-18 19:16-19:17; ICON untouched. Sixteen are chunk files inside the stores'
chunk grids and were unreadable (blosc decompression error): scream 5, um 4, nicam 3, casesm2 2, IMERGv7 2. The other three (IMERGv7) are
stale leftovers outside the store's current 24 x 262144 chunk grid (the store was rewritten on 2026-09-17), which no read ever touched. No
store had a missing chunk file, so nothing was silently read as zeros.

Restored on 2026-09-19 from the CFS backup (`/global/cfs/cdirs/wcm_shr/hk25/mcs/`) with `cp -p --remove-destination`, after checking that the
other 43,152 files of the five stores are md5-identical to CFS. The 19 files now match CFS in md5, size and mtime, are owned by feng045, and
a read of every chunk of all stores reports 0 unreadable. Fourteen of them were owned by wcmca1 with mode 644, so a process running as
feng045 could not have truncated them; that fits a filesystem-side event, cause unknown. The Step 1 outputs of 2026-09-18 (08:24-09:11) are
unaffected because they ran before the files were emptied. See `zero_byte_chunks.txt` and `restore_zero_byte_chunks.sh` in the audit
directory.

## Open follow-ups

1. **Step 3 category masks are not disjoint.** A track can be in two pair lists, and an "isolated" object can lie inside a 3-way category's
   whole footprint. Priority hides this in the monthly and extreme totals only; per-category counts and the A4 track statistics may still see it.
   Measured on the full SCREAM record (1577 frames, `overlap_full.py`, audit tooling): the precipitation that would be counted twice if the categories were
   simply added (E) is 2.18% of the total at 60S-60N (1.8-2.9% per month), 5.8% at 30-60N (up to 8.7% in northern winter), 4.95% at 30-60S and 0.04% at 20S-20N.
   93% of E is different objects overlapping, 7% is one object in two categories. By combination of categories (the ten largest combinations, 81% of E, so these are lower bounds): an isolated MCS inside a 3-way footprint
   (`mcs_iso` + `mcs_ar_etc`) is at least 45% of E; an AR in both the MCS-AR and AR-ETC lists at least 15%; an MCS in both MCS-AR and MCS-ETC 8%; MCS-ETC plus 3-way 9%.
   At track level (per 6-h frame): an AR is in both MCS-AR and AR-ETC lists at 0.26 per frame (4.4% of AR track-frames, 367 of 1577 frames); an MCS in both
   MCS-AR and MCS-ETC at 0.18 per frame (0.13% of MCS track-frames); an ETC bridging MCS-ETC and AR-ETC never (the existing promotion to 3-way works);
   an MCS in a 2-way and a 3-way list never. Options (they change definitions or per-category counts, so they need a decision):
   (A) also promote to 3-way when an AR or an MCS, not only an ETC, bridges two pairings, which removes the AR and MCS double-membership cases (about a quarter of E);
   (B) make "isolated" exclusive of 3-way footprints in Step 3, so that the per-category counts agree with what the priority order already does to the totals
   (targets the largest piece, 45% of E, and is not ambiguous); (C) leave Step 3 and keep the priority order in the consumers, as now.
   **Decision (2026-09-19): option C.** Step 3 is left unchanged for this round and the priority order stays in the consumers. To be followed up after this round, with A and B still open.
2. **Extreme thresholds.** `calc_extreme_precip_thresholds.py` still uses each source's 6-hourly product while the attribution uses `tot_pr`.
   One-month comparison (thresholds from `tot_pr` / from the product): UM, IMERG and CASESM2 identical (median ratio 1.000); NICAM same mean, per-cell
   correlation 0.94; ICON 0.56-0.88 at high and northern mid-latitudes (snow counted twice in the product); SCREAM 1.21-1.25 at the median.
   **Decision (2026-09-20): adopted for all sources**, implemented in the script itself as `--input_zarr` / `--input_var` (IMERG: the non-IR 6-hourly store). What follows is the comparison that led to it. As decided, the thresholds are computed from the `tot_pr` of the reprocessed Step 1 output (`thresholds_from_tot_pr.py`, audit tooling; it reuses
   `calc_precip_percentiles` and `write_netcdf` of the existing script and refuses to overwrite), saved to
   `/pscratch/sd/w/wcmca1/hackathon/extreme_precip_cof/` and compared with the existing files in `/pscratch/sd/w/wcmca1/hackathon/extreme_precip/`.
   SCREAM, full record (1578 windows): thresholds from `tot_pr` against the existing file (P90 / P95): median ratio 1.19 / 1.24 over all cells (mean 1.27 / 1.32; 10th-90th percentile 0.93-1.65 / 0.95-1.73), 74% / 78% of the cells differ by more than 10%, correlation 0.94; by band 1.28 / 1.36 at 0-30S, 1.38 / 1.44 at 0-30N, 1.10 / 1.12 at 30-60S and 1.12 / 1.15 at 30-60N. The polar bands differ most (60-90S correlation 0.03 / 0.01, where the liquid-only product holds almost no precipitation; 60-90N median ratio 0.98 / 1.02, correlation 0.70 / 0.72). The attribution with the new thresholds is compared with the existing ones in the last section. The same comparison is to be done for the other five sources once their Step 1 is rerun.
3. **SCREAM 6-hourly product** (`scripts/avg_pr_healpix.py` -> `scream_pr6h_z8.zarr`, still read by the threshold script and by
   `extract_environments/extract_etc_2d_vars.py`):
   - **Window alignment.** It is built from `scream_ne120`, 3-hourly averages that are labelled at the end of each interval (the first label is
     2019-08-01 03:00). `resample(time='6h', label='left', closed='left')` put the means labelled T and T+3 h in the window labelled T, which covers
     [T-3 h, T+3 h) instead of [T, T+6 h), the window Step 1 aggregates over. Matching the window to Step 1's raises the agreement with the 6-hourly
     `tot_pr`: over the full record (392 windows, every 4th) the correlation with `tot_pr` at the same label is 0.86 in the tropics and 0.90 at 30-60N for the aligned product, against 0.74 and 0.72 for the old one; the aligned product peaks at that label only (0.47 / 0.37 one window either side) while the old one has a second peak one window later (0.62) (an earlier 16-window measurement gave 0.74 -> 0.87 in the tropics). The rest of the disagreement is not timing (see the last bullet).
     Fixed in `avg_pr_healpix.py`: `time_label='end'` moves the stamps back by half the interval in the model's own calendar, before the conversion to standard
     dates (a plain `closed='right'` leaves two wrong windows at the noleap gap), and the result is written to a NEW file, `scream_pr6h_z8_aligned.zarr`;
     `scream_pr6h_z8.zarr` is untouched. `extract_etc_2d_vars.py` pairs the product with an instantaneous track time (`sel(time=storm_time, method='nearest')`),
     for which the window centred on T may be preferable, so it still reads the old file until that is decided.
   - **Liquid only.** The `prs` line in `vars_to_include` is commented out (the 3-hourly `pr` is `precip_liq_surf_mass_flux`).
   - **Edges.** The noleap calendar mapped to standard dates leaves four all-NaN windows on 2020-02-29 (the aligned product has 4 empty and 0 partial windows;
     the old one had a single 3-hourly mean in its first and last windows). `resample().mean()` does not check that a window is complete; the script now logs it.
   - **The two products differ in spatial structure.** The processing chains, as understood (catalog notes plus the analyst's description): the hourly stream is native ne1024 (about 3.25 km) output
     Delaunay-remapped to zoom 10 and conservatively coarse-grained to level 8; `scream_ne120` is an online SCREAM output stream at ne120 (probably a
     conservative coarsening of ne1024) remapped to level 8 with Delaunay. Compared at the SAME instant (32 instants from 2020-07-11; hourly stream at label
     L against `scream_ne120_inst` valid at L), in the tropics / 30-60N:

     | field | P99 (mm/h) | wet fraction (>= 0.1 mm/h) | mean (mm/h) |
     |---|---|---|---|
     | hourly stream at L | 19.2 / 13.9 | 11.6% / 7.4% | 0.1788 / 0.0900 |
     | 3-hourly instantaneous at L (ne120 chain) | 10.1 / 8.6 | 15.9% / 9.7% | 0.1788 / 0.0899 |
     | 3-hourly average labelled L (ne120 chain) | 8.2 / 6.6 | 18.2% / 11.2% | 0.1786 / 0.0894 |

     The pixel correlation between the two fields at the same label is 0.77 (0.83 at 30-60N) with identical means, so the spatial chains differ substantially before any
     averaging is involved; averaging inside the ne120 chain (instantaneous to 3-h average) trims P99 by about 20% more. Which step inside each chain causes
     the difference is not isolated. Window length and label alignment do not explain it. The earlier 6-hourly comparison (6-h means built from the
     3-hourly product against 6-h means of the hourly `pr`) showed the same: wet fraction 20% against 16% and P99 6.8 against 9.8 mm/h in the tropics.
4. **IMERG.** Its AR, TC and ETC masks come from ERA5, so its Step 2 is `combine_era5_imerg_tracking_masks.py` (no command-line options; it carries
   `tot_pr`); `slurm/tasks_all_IR_IMERG.txt` lists only Step 3 and the monthly script.
5. **Partial windows.** Windows with fewer than six hourly steps are written from the steps present, including single-step windows at gaps and record
   ends (SCREAM 2020-04-20T00 and 2020-09-01T00, CASESM2 2021-03-01T00). No minimum number of steps is enforced. Inventory of the Step 1 input time axes:
   SCREAM 4 such windows, NICAM 1, CASESM2 2, IMERG 6 (one with 3 steps, 2019-09-28T18), ICON 0, UM 0. Before the key-mismatch fix, the windows whose first hourly
   step was off the aligned hour were dropped: SCREAM 2019-08-01T00 and 2020-04-22T00, NICAM 2020-03-01T00, CASESM2 2020-03-01T00 and IMERG 2019-03-16T00.
   They become valid frames when those sources are rerun.

## Full-scale SCREAM re-process (2026-09-19, nid004154; test area, production untouched)

The redirected copies of Steps 1-3, the monthly script and the attribution script ran on the whole SCREAM record with the committed code
(`c512a93`, `pyflex-dev`), outputs in `/pscratch/sd/w/wcmca1/hackathon/tmp/scream_full/`.

| Check | Result |
|---|---|
| Step 1 | exit 0 in 15.5 min with 48 workers; 1578 of 1578 windows written; 0 failed, 0 empty, no retry needed; 4 partial windows listed (2019-08-01T00 and 2020-04-22T00 with 5 steps, 2020-04-20T00 and 2020-09-01T00 with 1) |
| Step 1 against production | identical in every cell of every frame for `mcs_mask`, `cloud_types` and `dc/st/nd/dz_pr`, except the two formerly dropped frames (2019-08-01T00 and 2020-04-22T00: empty in production, valid now) |
| `tot_pr` | no NaN; `dc+st+nd+dz = tot_pr` at classified non-swath cells to 8e-6 mm/h; no rain without a cloud type |
| Frame ledger | Step 1 1578 -> Step 2 1577 -> Step 3 1577 -> monthly `ntimes` 1577, identical to production; the one frame that drops (2020-09-01T00, a single-step window) drops in both |
| Steps 2, 3, monthly | exit 0 (37 s, 6.7 min, 1.4 min); `tot_pr` identical in Steps 1-3; COF `mcs_mask` equals Step 2's in every frame (Step 3 removes nothing); against production the new COF store has 0.25% more MCS pixels (370 track-frames in 11% of the frames), never fewer |
| Monthly closure | worst wet cell 0.00006% in any of the 13 months (production: domain-total residual -1.2% to -3.7% per month) |
| Boxes, whole record, residual before -> after | 60S-60N -2.30%, Sahel +1.54%, NW N. America -9.51%, W. China/Russia -20.00%, NH 30-60N -8.50%, SH 30-60S -8.59%, deep tropics +0.94% -> 0.0000% in all |
| Notebook definition (mean of cell fractions), before -> after | tropics land +0.85% / ocean +1.26%; N.H. land -13.32% / ocean -4.69%; S.H. land -7.94% / ocean -8.84% -> 0.0000% everywhere (map: `scream_residual_maps_old_new.png` in the audit tooling folder) |
| Total precipitation, new / old | 1.046 at 60S-60N (3.124 -> 3.267 mm/day); W. China/Russia 1.352, NW N. America 1.266, SH 30-60S 1.145, NH 30-60N 1.125, deep tropics 1.000: snow is now in the total |
| Extreme attribution (P90 / P95) | unassigned share 0.192% / 0.206% in production -> 0.000% with either threshold file; shares sum to 1 and the type counts add up in every cell |

The extreme attribution shows why thresholds and precipitation must come from the same field. With the existing thresholds and the new precipitation, extreme precipitation is
1.39 / 1.61 times production, the drizzle share falls from 1.6% to 0.005% (P90; consistent with the old drizzle share coming from the product mismatch) and the tracked categories
drop (MCS-AR-ETC 11.0% -> 9.5%, TC 3.9% -> 3.0%). With thresholds from `tot_pr` the total is 1.16 / 1.17 times production and the tracked categories stay close to production
(P90: isolated MCS 51.7% against 54.3%, MCS-AR 6.1% against 6.4%, MCS-AR-ETC 10.6% against 11.0%, TC 3.3% against 3.9%), while deep convective rises from 9.5% to 14.1% and
non-deep convective from 5.2% to 6.9%. The Step 3 overlap on this store is in follow-up 1.

Is `make_mcs_swath_masks.py` ready for the other sources? Yes. The SCREAM record ran end to end with exit 0 in every stage (24 min in total on the node; Step 1 15.5 min at 48 workers) and no retry was needed in 1578 windows, so retries under a real failure are covered by the injection tests only. Steps 1-3 use the same code for the other sources; all six Step 1 configs carry TC settings.

## Full-scale re-process of the other five sources (2026-09-19/20, nid004150; test area, production untouched)

Same committed code and procedure as SCREAM (`pyflex-dev`, redirected copies, Steps 1-3 and the monthly script chained; IMERG Step 2 is
`combine_era5_imerg_tracking_masks.py`). Outputs are in `/pscratch/sd/w/wcmca1/hackathon/tmp/full_reprocess/`; thresholds from Step 1's `tot_pr` are in
`/pscratch/sd/w/wcmca1/hackathon/extreme_precip_cof/`, and the production files were not touched. All stages exited 0 for every source. Only `pr` is used (frozen
precipitation is ignored), so the UM Step 1 fields differ from production where production had added `prs` (92-99% of the differing cells are poleward of 40 degrees).

| Source | Windows | Step 1 vs production | Monthly closure | Notebook residual, \|lat\|<=60, land / ocean, before -> after | Total precipitation new / old |
|---|---|---|---|---|---|
| ICON | full record | identical | worst wet cell <= 0.00006%, residual 0.0000% | +8.48% / +2.83% -> 0.0000% | 0.959 (0.737 in W. China/Russia: the snow counted twice is gone) |
| NICAM | full record | identical apart from 2020-03-01T00 (formerly dropped) | 0.0000% | -0.78% / -2.32% -> 0.0000% | 1.000 |
| CASESM2 | full record | identical apart from 2020-03-01T00 (formerly dropped) | 0.0000% | -0.13% / -0.64% -> 0.0000% | 1.000 |
| UM | 1576 frames | differs from production where production added `prs` | worst wet cell 0.00005%, 0.0000% | -35.41% / -4.72% -> 0.0000% | 1.000 |
| IMERG | 4384 (36 months) | identical apart from 2019-03-16T00 (formerly dropped) | **not zero**, see below | -0.69% / -2.54% -> +0.16% / +0.14% | 1.000 |

The NICAM and CASESM2 frame 2020-03-01T00 becomes valid in Step 1 and is dropped by Step 2 in both the old and the new run, so the frame ledger (Step 1 -> Step 2 -> Step 3 -> monthly) is unchanged
in all five sources. The IMERG frames are the same in the old and new run at every stage.

**IMERG does not close to zero.** The domain residual over the 36 months is +0.14% at 60S-60N (production -2.28%), +0.05% to +0.60% per month (largest in 2019-01 to 2019-05), and single wet cells reach
up to 100% in the worst month. The notebook fractions are +0.16% land / +0.14% ocean at 60S-60N. The cause is in Step 1, not in Steps 2-3 or the monthly script: the Step 1 store is identical to production, and
at cells that have a cloud type outside the MCS swath `dc+st+nd+dz` differs from `tot_pr` by up to 16.7 mm/h in single frames, with 472,591 pixel-frames (0.014% of all) that have rain but no cloud type.
By the classification code, an hour with rain and no cloud type outside an MCS can only be one whose Tb is missing (all four conditions compare `tb` with the threshold). Measured afterwards on the
input Step 1 actually reads, `IR_IMERG_V7_1H_zoom8_20190101_20211231.zarr` (Step 1 overrides the config's `zoom: 9` with the mask zoom, `make_mcs_swath_masks.py` line 1272): between 60S and 60N `precipitation`
is never NaN, and `Tb` is NaN in 0.30% of cell-hours; in a sample of every third day, 0.33% of the in-band rain falls on NaN-Tb cells, and the monthly share follows the monthly residual (r = 0.85; the share
is 0.06% to 1.2% per month, largest in 2019-01 to 2019-05, 2019-09, 2020-02, 2021-09 and 2021-12). The share is larger than the residual, presumably because part of that rain lies inside MCS swaths, where it is attributed
to the MCS whatever Tb is; that part was not checked. The fix is a decision about what to do with rain at missing-Tb pixels (for example put it in a separate "no Tb" category or count it as drizzle/non-deep),
not a change in the closure logic. **Resolved 2026-09-20:** Step 1 removes the rain at pixels with missing Tb (d2bcaab); results in the next-to-last section.

**IMERG input coverage.** The Step 1 input has data only within 59.87S-59.87N: poleward of 60 both `Tb` and `precipitation` are NaN in every cell and hour (not zero). Step 1 turns missing hours into zero
(`tot_pr` is the sum with NaN as 0 over the hours present), so `tot_pr` is 0 there and the `tot_pr`-based thresholds are empty poleward of 60 (no wet sample). The existing IMERG thresholds were built from
`IMERG_V7_6H_zoom8_20190101_20211231.zarr` (`calc_extreme_precip_thresholds.py`, the non-IR IMERG product, built from `IMERG_V7_1H_zoom9`), which has data at all latitudes (rain frequency 10% at 60-70N, 3-6% at
70-90N); within 60S-60N it is the same field as `tot_pr` (same window labels; at the same label r = 0.9998 and sum ratio 1.0000, at the neighbouring windows r = 0.43). No chunk is missing in either store (3288 of 3288
per variable in the IR 1H zoom-8 store, none zero-byte or truncated; 549 of 549 in the non-IR 6H store). The IR 1H zoom-9 store, which the config names, has lost `.zgroup`, `.zattrs`, `.zmetadata` and the per-array
`.zarray` files (chunk file count is complete, 13152 per variable) and cannot be opened by zarr as it is; Step 1 does not read it.

**Thresholds and attribution, `tot_pr`-based against the existing files (median ratio; unassigned share of extreme precipitation at 60S-60N, production | new code with existing thresholds | new code with `tot_pr` thresholds):**

| Source | Median ratio P90 / P95 | Cells differing > 10% (P90) | Unassigned P90 | Extreme total vs production, existing thr. / `tot_pr` thr. (P90) |
|---|---|---|---|---|
| SCREAM | 1.19 / 1.24 | 73.9% | 0.192% \| 0.000% \| 0.000% | 1.387 / 1.155 |
| ICON | 1.01 / 1.01 | 44.8% | 0.218% \| 0.000% \| 0.000% | 0.988 / 0.965 |
| NICAM | 1.00 / 1.00 | 19.1% | 0.234% \| 0.000% \| 0.000% | 1.003 / 1.000 |
| UM | 1.00 / 1.00 | 2.5% | 0.102% \| 0.000% \| 0.000% | 1.001 / 1.003 |
| CASESM2 | 1.00 / 1.00 | 0.0% | 0.281% \| 0.000% \| 0.000% | 1.000 / 1.000 |
| IMERG | 1.00 / 1.00 | 0.3% | 0.176% \| 0.016% \| 0.016% | 1.000 / 1.000 |

Only SCREAM (its thresholds came from the 6-hourly product built from another data stream) and, to a lower degree, ICON (snow in `tot_pr`) have thresholds that move with `tot_pr`. For NICAM, UM, CASESM2
and IMERG the existing thresholds are already what `tot_pr` gives. Shares sum to 1 and type counts add up to `total_extreme_count` in every source (validation: PASS x5). No decision to adopt the `tot_pr` thresholds has been made.

**Correction to an earlier note.** The first attempt to run the five chains at once ended when the Slurm job died (exit 137). I first attributed this to memory exhaustion. A measurement on the next node
(24 workers = 14 GB resident; five chains together = 131 GB of 503 GB) does not support that, so the cause of that loss is unknown. The runs afterwards were staggered and completed without incident.

## Round 2: the pipeline runner on all six sources (2026-09-20, nid004144; test area `tmp/round2`)

All six sources, Steps 1-3, the monthly map, the thresholds (from `tot_pr`; IMERG from the non-IR 6-hourly store) and the attribution, in one run of `slurm/run_interactive_cof_pipeline.sh` with the committed code
(378ecf4): 36 steps, all exit status 0, 94.6 min wall time, into `/pscratch/sd/w/wcmca1/hackathon/tmp/round2/` (production untouched at that time; the outputs were promoted to production on 2026-09-20 evening, follow-up 3). Time and memory per step are in `docs/procedures/run_cof_pipeline.md`.

| Check | Result |
|---|---|
| Five models against round 1 (manual chains), every stage | Bit-identical: Steps 1, 2 and 3 (all 7, 11 and 22 variables, every frame), the monthly file (43 variables), the thresholds (P90, P95; same NaN cells) and both attribution files. The orchestrated run equals the manual chains, and the Step 1 rule for missing Tb changes nothing for the models |
| IMERG against round 1, Steps 1-3 | Only `tot_pr` differs (8,201,903 cell-frame values, 0.24% of all, in all 4384 frames, at most 26 mm/h); `mcs_mask`, the cloud types and `dc/st/nd/dz_pr` are identical; the frames are identical in every stage (ledger: 4384, none dropped) |
| IMERG Step 1 | `max\|dc+st+nd+dz - tot_pr\|` at classified non-swath cells 3.8e-6 mm/h (16.7 mm/h in round 1); pixel-frames with rain but no cloud type 0 (472,591); identical to production apart from the formerly dropped frame 2019-03-16T00 |
| IMERG monthly closure | Worst wet cell 0.00005% in every one of the 36 months (100% in round 1); domain residual 0.0000% in every month (+0.14% at 60S-60N before); the notebook's residual fractions 0.0000 in tropics, N.H. and S.H., land and ocean (production -0.69% land / -2.54% ocean at 60S-60N) |
| IMERG total precipitation | 0.33% lower than production and round 1 (60S-60N mean 3.0026 to 2.9926 mm/day; 0.1-1.4% per month, largest in 2019-01 and 2019-02): exactly the rain that fell on pixels without Tb |
| IMERG thresholds | Bit-identical to the production file (773,094 finite cells, same NaN cells, 87.6% of the cells poleward of 60 have a threshold); 60S-60N equals the IR-based `tot_pr` field (correlation 0.9998 at the same window label) |
| Attribution, all six sources | Shares sum to 1, type counts add up, unassigned 0.000% at P90 and P95 (IMERG 0.176% / 0.183% in production, 0.016% / 0.015% in round 1). IMERG's type shares move by at most 0.12 percentage points against round 1 (isolated MCS at P90 56.20% to 56.32%); total extreme precipitation 0.996 / 0.995 of production |

The IMERG monthly counts of precipitating hours (for example `ar_etc_precipitation_count`) change at some cells, because a window whose rain was only at pixels without Tb no longer counts as precipitating.

## Follow-ups after this phase

Held on purpose while the pipeline runner and the second full round are done; none of them changes the results of Analyses 1 and 2 as they stand.

1. **`extract_etc_2d_vars.py` and the SCREAM 6-hourly file.** It still reads `scream_pr6h_z8.zarr` (windows centred on the label, as built) and pairs it with an
   instantaneous track time through `sel(time=storm_time, method='nearest')`. The aligned file `scream_pr6h_z8_aligned.zarr` (window [T, T+6 h), the window Step 1 uses) exists
   but nothing reads it. Which one the ETC composites should use is a question for that analysis (Analysis 3); the extreme thresholds no longer read either file.
2. **Step 3 pair-list overlap** (follow-up 1 above, options A and B). Option C stays: Step 3 is unchanged and the priority order stays in the consumers. The impact
   measured on SCREAM is small (2.2% of the precipitation at 60S-60N counted twice if the categories were added, 0% in the totals because of the priority order).
   **Decided 2026-09-20: option A (bridge promotion) is not adopted.** An ETC is the main dynamical driver of the AR and MCS it touches, so the existing ETC-bridge rule is physically justified; promoting
   through an AR or MCS bridge is less so, and the orphan case is deliberately left as an edge case. Option B was not part of the decision and stays open. The plan, the measured effect, the rationale and
   what a revisit would need are in [step3_bridge_promotion_future_work.md](step3_bridge_promotion_future_work.md).
3. **Promotion of the test-area outputs to production.** *Done 2026-09-20.* The 54 round-2 items (`mcs_masks`, `all_masks`, `cof_masks` with `stats/` overlap files and `stats/monthly/`, `extreme_precip`; 192.6 GB) were moved by rename from `tmp/round2`
   into the production folders of `/pscratch/sd/w/wcmca1/hackathon/`. The previous production (119.7 GB, made 2026-09-18/19; thresholds 2026-03) and the round-1 interim thresholds (`extreme_precip_cof`) are in
   `/pscratch/sd/w/wcmca1/hackathon/_prev_production_20260918/` (README, `MANIFEST_pscratch.tsv`, sizes). Verified after the move: sizes of all 108 items equal the sizes measured before, the exact chunk-key scan finds all 18 zarr stores complete
   (289 arrays), the frame counts are 1460, 1460, 4384, 1459, 1577, 1576 in Steps 2 and 3 (Step 1 holds one extra final window for CASESM2, NICAM and SCREAM, as before), every nc and parquet file opens. The ERA5 mask store, the March
   Analysis 4 files and everything else in those folders were not touched. The round-2 runner markers and logs stay in `tmp/round2`.
   The CFS backup folders were prepared the same way: the 54 older same-named items (121.1 GB, dated 2025-11 to 2026-07) moved into `/global/cfs/cdirs/wcm_shr/hk25/_prev_production_20260918/`, so that a Globus copy of the
   new pscratch items mirrors production exactly (the ERA5 mask store, which the IMERG Step 1 reads from CFS, stays). The list of the 54 new items is `NEW_ITEMS_for_globus.txt` in the pscratch archive folder.
4. **The unexplained loss of the first node** (Slurm job killed with status 137 82 s after five chains were started at once). The cause is unknown; the runner staggers starts,
   gates on available memory and records the memory of every step, which should show it if it happens again.
5. **`--min_precip_threshold`.** *Done.* The default of `calc_extreme_precip_thresholds.py` was 0.01 mm/h while the existing threshold files were made with 0.1 mm/h (through the run script); it is now 0.1
   like the 1-hourly script, and the runner and `run_all_extreme_precip_thresholds.sh` still pass it explicitly. No stored file changes; only a run without the flag would have differed.
6. **`write_zarr` for short stores.** *Done.* `src/zarr_tools.write_zarr` (Step 2) and the writer inside `combine_era5_imerg_tracking_masks.py` divided by zero when the store had fewer time steps than
   the time chunk (24 and 28); it only mattered for tiny test stores (a real record has more than 1400 frames). Both now skip the "make the chunks more even" adjustment when there is no full chunk
   (`chunks > 0`) and write one chunk. Every store that worked before is written with the same chunks and the same data (checked for 12 lengths from 1 to 1578 frames in each writer); the SCREAM test chain with 24 test steps that failed in round 2 now runs through Steps 1-3.
7. **IFS** (`ifs_tco3999_rcbmf`) is in `config_sources.yaml` and has a Step 1 config but no tracking inputs yet; it is not enabled in `config/config_pipeline.yaml`.
