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
- The test-area outputs were not promoted to production at that time; round 2 was promoted on 2026-09-20 evening (follow-up 3 at the end).
- Held for after this phase, settled on 2026-09-20 (see "Analysis 3" and "Follow-ups after this phase" at the end): Analysis 3 reads Step 1's `tot_pr` for the models instead of a 6-hourly SCREAM file; the Step 3 pair-list overlap stays as it is (option C, and the bridge promotion is not adopted).

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

## Analysis 3 (ETC composites): review, fixes and the core-tier run (2026-09-20, node nid004147; test area `tmp/etc_round1`)

Analysis 3 composites 2D fields (161 x 161 points at 0.25°) around every ETC track point and stratifies them by the Step 3 overlap flag of the ETC (0 isolated, 1 MCS only, 2 AR only, 3 MCS and AR). The review covered the
four steps (`combine_etc_cof_data.py`; the extraction `extract_etc_2d_vars.py`; `combine_etc_2d_vars.py`; `create_etc_composites.py` and `calc_etc_spatial_stats.py`), how the precipitation of each source is matched to the track times, and what was on disk.
The runner, the commands and the time matching are in [procedures/run_etc_pipeline.md](procedures/run_etc_pipeline.md). The checks below are audit tooling in `/global/homes/f/feng045/.claude_tmp_verify/round3/`
(`verify_core_tier.py`, `verify_vs_cof.py`, `mask_diff.py`, `compare_stats2.py`, `radial_pr.py`, `compare_figures.py`, `panel_metrics.py`, `prep_notebooks.py`); the outputs are in `/pscratch/sd/w/wcmca1/hackathon/tmp/etc_round1/`.

### Findings and what was done

| # | Finding | Resolution |
|---|---|---|
| 1 | **The `pr` of SCREAM, ICON and (slightly) NICAM is not the window of Steps 1-3.** Correlation with Step 1's `tot_pr` at label offsets -6 h / 0 / +6 h (three 4-day blocks per store, \|lat\| < 55): UM, CASESM2 and IMERG 1.000 at 0 h; NICAM 0.94-0.95 at 0 h; SCREAM 0.72-0.74 at 0 h (window centred on T, straddles two windows, liquid only); ICON 0.967-0.969 at -6 h (end-stamped: the value at T is the mean of [T-6 h, T), a full window early) | `pr` of the five models is Step 1's `tot_pr` (mean over [T, T+6 h), snow included), read from the COF store; ERA5 keeps IMERG 6-hourly (`tot_pr` of IMERG is 0 poleward of 60° and ETCs reach 70°). Decided 2026-09-20 (77fb876) |
| 2 | **`sel(time, method="nearest")` never fails:** a track time outside a source's record silently got the first or last frame. The COF stores of ICON (one year), UM and SCREAM are shorter than the track periods. On the extraction's point list: ICON 3299 of 22594 points (14.6%: 2020-01-01 and 2021-01/02), UM 749 of 23287 (3.2%: 20-31 January 2020), SCREAM 121 of 22301 (0.5%: seven consecutive windows, 20 April 06 UTC to 21 April 18 UTC 2020); NICAM, CASESM2 and ERA5 none. In the March data these points have `pr` from the catalog and masks that are NaN (ICON) or partly finite (UM 36%, SCREAM 12%) | exact time lookup: a track time without a frame gives NaN, is counted and recorded in the store attributes (`n_points_time_missing`), and the job fails when more is missing than `--max_missing_fraction` (registry: SCREAM 0.008, ICON 0.15, UM 0.035, others 0) (77fb876) |
| 3 | **The `etc_data` stores on pscratch were empty directory skeletons** after the scratch purge (no `.zgroup`, `.zmetadata` or `.zarray`; all `single_vars/` stores and the six combined stores), and **`scripts/check_zarr_store.py` reported such a store as complete** (it found no arrays, so nothing was missing). The March stores were restored from CFS with Globus on 2026-09-20 (store names, metadata and the `psl` point counts 22301, 22594, 21563, 23287, 17593, 66789 checked) | the checker raises when there is no array metadata or a requested array is absent (18764bb); the runner's preflight uses it and also compares the storm-point count of every environment store with the track file |
| 4 | **The combine step would scale `tot_pr` twice** (×3600, SCREAM ×3.6·10⁶, for a flux that is already in mm h⁻¹); it warned about a store with another point list and merged by position; a variable that could not be loaded was skipped | the factor is applied only to a native flux (units check, `scale_factor_applied` in the attributes); the mean of `pr` must lie in 0.01-1 mm h⁻¹; count, times, `storm_id`, lat, lon and `grid_id` of all single-variable stores must be equal; loading is strict unless `--allow-missing-variables`; exit status (772370e) |
| 5 | The COF masks were read from a hard-coded production path; the pipeline overview said 81 x 81 points (the data are 161 x 161, radius 20°); Step 3 could exit 0 with frames missing | `COF_DATA_ROOT` through `src/cof_paths.data_root()` (77fb876); documentation corrected; Step 3 fail-loud (cf695bc, follow-up 8) |
| 6 | Analysis 3 had no runner: Step 1 and Steps 3-4 by hand or by `run_calc_etc_spatial_stats_all.sh` (hard-coded source list and paths), Step 2 as one Slurm array per variable group | `scripts/run_etc_pipeline.py` with `config/config_etc_pipeline.yaml`: the same engine as `run_cof_pipeline.py` (scheduler, memory gate, markers, `--resume`, overwrite protection); 13 tasks per source (72a836e). The environment variables are not extracted again: the March stores are linked after the checks, so the online-only UM and CASESM2 catalogs are not touched |
| 7 | The notebook `plot_etc_spatialmean_stats_1source.ipynb` cannot run top to bottom: the cell that defines `track_stats` contains `... / 100track_stats.head()` (a missing line break, a SyntaxError). Its saved figures are not reproduced by the March statistics (stale) | not changed in the repository; fixed in the run copy only. Comparisons therefore use a re-run on the March data as baseline |

### Not changed on 2026-09-20 and what became of it

- **The radii of the spatial statistics were in grid points, not degrees.** `x` and `y` of the combined store are offsets in grid points (-80 to 80, 0.25° apart), and `create_circular_mask` compared them with the radius as given, so the "10°" and "15°" radii (`--basic-radius`, `--pr-radius`, `--mask-radii`,
  `*_frac_r10`, `*_frac_r15`, the attributes `basic_radius_deg`, `pr_radius_deg`, `mask_radii_deg`, the notebook's printed "10° radius") were 2.5° and 3.75° (checked: `pr_mean` of the statistics equals the mean of `pr` inside a circle of 10 grid points, correlation 1.0000, not inside 40). **Fixed 2026-09-21 (2d2b4e5);** see the second round below.
- **The "all" composites contain the points that have no frame of the source.** Their `pr` is NaN and skipped by the mean, but the mask frequencies (`xr.where(mask > 0, 1, 0)` turns NaN into 0) and the exclusive precipitation (`pr_ar`, `pr_mcs`, `pr_etc`, 0 where the mask is NaN) count them as zeros, so `*_freq` and `pr_ar/pr_mcs/pr_etc` of the `all` files, and the
  fractions `prfrac_mcs` and `prfrac_ar` made from them, are low by the share of those points in the hemisphere. Measured with the composites' own selection (`cof_lat` > 20 / < -20, before the latitude rule): ICON 16.7% NH (1840 of 11050 points) / 12.6% SH (1459 of 11544), UM 3.7% / 2.8%, SCREAM 0.6% / 0.5%, NICAM and CASESM2 none
  (the 1.5% and 4% given first for SCREAM and UM were shares of points with a NaN flag by the sign of the storm latitude, not of points without a frame). The composites of a single flag (`isolated`, `mcs_only`,
  `ar_only`, `3way`) and the statistics filtered by `overlap_flag` are not affected (verified directly 2026-09-22: a frame-less point always has `overlap_flag` NaN, and `NaN == flag` is always `False`, so the single-flag selection excludes it as a side effect with no explicit filter needed — 0 of the frame-less points in any of the six sources have a non-NaN flag); the plots that use `overlap_type = 'all'` are (frequency contours, `pr_mcs`, `pr_ar` shading and the two fractions; `pr` itself is not, nor are the environment/meteorological fields, which have no gaps). The fix is about three lines in `create_composites` (a `has_frame` selection for the `all` case); **not changed**, needs a decision.
  **Re-measured 2026-09-22 against today's production files** (`--min-lat-coverage 0.7`, i.e. after the latitude rule, not before as above — the composite files are built from the post-rule sample, so this is the share that actually determines the bias in the file): ICON 20.5% NH (n=6121) / 14.0% SH (n=3836), UM 4.6% NH (n=5680) / 2.1% SH (n=5022), SCREAM 0.5% NH (n=5777) / 0.3% SH (n=4262), ERA5/NICAM/CASESM2 0%. The relationship is exact at every pixel: `stored = frame-only-mean × (1 − frame-less share in that file)` — checked on ICON NH `mcs_freq`: stored 0.1805, frame-only-recomputed 0.2272, ratio 0.7946 = 1 − 0.2054. Each file records its own `n_points` and `n_points_hemisphere_before_lat_rule` attributes, so the current share for any file/source can be recomputed directly rather than read off this dated snapshot.
- **The masks of the combined stores have 0 where there is no object, the March stores NaN** (the current Step 3 writes 0). Composites and statistics use `mask > 0`, so this changes nothing, but a comparison must treat NaN as 0 (`mask_diff.py` does).

### Verification (all on nid004147)

| Check | Result |
|---|---|
| Unit tests | `tests/test_etc_time_matching.py` (the old code returned frames 4, 9 and 0 for a gap, a time after and a time before the record; the new code NaN, NaN, NaN), `test_etc_combine.py`, `test_check_zarr_store.py`, `test_link_etc_env_stores.py`, `test_run_etc_pipeline.py`, `test_step3_fail_loud.py`; the runner's `--self-test` and the COF runner's dry run are unchanged |
| The core-tier run, six sources | 78 steps (13 per source), 0 failures: extraction phase 60 steps, 19.3 min; combine, composites, stats 18 steps, 9.9 min (second invocation with `--resume`). Combined stores have the March point counts 66789, 22301, 22594, 21563, 23287, 17593; points without a frame 0, 121, 3299, 0, 749, 0 (as predicted); at these points `pr` and all masks are NaN and `overlap_flag` is NaN; the store attributes carry the counts |
| New stores against the round-2 COF store (independent window code, 60 points per source, 30 random and 30 with flag 3; 0.5-0.8 million mask pixels > 0 per source) | all seven masks equal pixel by pixel, and `pr` equal to `tot_pr` for the five models, at every sampled point; the sampled points without a COF frame are all NaN. The extraction reproduces the round-2 products exactly |
| `pr` against the March store at the same points (every 5th point that has a frame) | correlation / ratio of the means: UM 1.00000 / 1.0000 (0.000% of pixels differ by more than 1e-4 mm/h), CASESM2 1.00000 / 1.0000 (0.000%), ERA5 1.00000 / 1.0000 (0.000%; IMERG 6-hourly, unchanged), NICAM 0.953 / 1.001, ICON 0.459 / 0.990, SCREAM 0.712 / 1.343. These are the values that the label-offset table predicts (ICON 0.40-0.48 and SCREAM 0.72-0.74 at 0 h, NICAM 0.94-0.95) |
| Where the SCREAM and ICON `pr` differ (mean `pr` by distance from the centre, new against March) | SCREAM +16% within 5°, +37% to +42% at 5-20°; by storm latitude +82% (70S-50S), +14% (50S-30S), +3% (30S-30N), +13% (30-50N), +37% (50-70N): snow, which the old file lacked, and the aligned window. ICON -27% within 5°, +6% at 10-20°, the box mean -1.4%: the same rain, one window later, so it now lies under the masks of the same window (under the MCS mask +59% / +42% in the 3-way composites (NH / SH), under the AR and ETC masks -13% to -24%). NICAM +6% within 5° and within 1.5% elsewhere |
| Masks and flags against March (every 5th to 10th point that has a frame) | the footprints (mask > 0) are larger by 0.4-0.9% (ERA5) and by up to 7% (models: -0.3% to +7.2%, most 1-5%); 0.15-3.5% of the points have any pixel that differs; the track ID differs at at most 0.55% of the pixels where both have an object. `overlap_flag` differs at 7 (ERA5), 16 (SCREAM), 2 (ICON), 4 (NICAM), 0 (UM), 2 (CASESM2) points. These come from the round-2 Step 3 products (follow-up 3), not from Analysis 3 |
| Composites and statistics against March (3-way = the composites of the multipanel notebook; NH / SH) | **UM and CASESM2: `pr` identical** (composites and statistics). **ERA5: `pr` identical** (-0.0% in the 3-way composite, 1.000). **NICAM:** `pr` -0.1% / +0.1% (pattern correlation 0.994 / 0.993; near the centre +6% within 5°, within 1.5% elsewhere). **SCREAM:** `pr` +15% / +24% (0.95 / 0.94) in the 3-way and +23% / +47% in the `all` composite; in `all` the rain under the MCS mask rises (`pr_mcs` +57% / +75%, `pr_etc` +28% / +59%) and under the AR mask falls (`pr_ar` -12% / -7%). **ICON:** `pr` +1.0% / +0.6% (0.86 / 0.82) but `pr_mcs` +59% / +42%, `pr_ar` -18% / -13%, `pr_etc` -24% / -24% (3-way). The mask frequencies (`ar_freq`, `mcs_freq`, `etc_freq`) and the exclusive precipitation change with the masks: +0.2% to +1.2% (ERA5), +1.2% to +2.4% (NICAM), +1.8% to +7.4% (SCREAM), +2.5% to +3.0% (ICON), +2.6% to +7.1% (UM), +1.7% to +6.4% (CASESM2). The environment variables (not extracted again) change only through the few flags that changed: composite means within about 1.5% with pattern correlation of at least 0.999, more only where the mean is near zero (ω at 850 hPa in SCREAM's AR-only composite +13%, NICAM's isolated `vas`). Statistics: ERA5, UM and CASESM2 change in the mask-derived variables only (at most 3.4% of the points; UM's `pr_nonmcs_mean` through the 749 NaN points); SCREAM, ICON and NICAM in all `pr` statistics (`pr_mean` within 2.5° of the centre: SCREAM +45%, ICON -14%, NICAM +10%, mean over points where both are finite) |
| Notebooks on the new data (run copies; figures in `figures_etc/multipanel_new` and `1source_era5_new`, the repository notebooks untouched) | `plot_etc_composites_multipanel.ipynb` (NH and SH, 3-way, six sources): the March re-run reproduces the 10 existing PNGs of `/global/cfs/cdirs/m1867/zfeng/hk25/figures_etc_composites/multipanel/` pixel for pixel; new against existing, share of pixels that differ by more than 12 grey levels in each panel (mean of the 10 figures): ERA5 3.3%, CASESM2 7.2%, UM 7.8%, NICAM 11.2%, SCREAM 20.3%, ICON 21.5% (about 30% in the two rain-fraction figures for SCREAM and ICON). ERA5, UM and CASESM2 differ only by small shifts of the AR and MCS frequency contours; the precipitation and rain-fraction shading changes in SCREAM and ICON. `plot_etc_spatialmean_stats_1source.ipynb` (ERA5): all 24 code cells run without an error, 9 figures, 7 identical to the baseline and 2 differ in 0.1% and 0.7% of the pixels (shifts of a few curve segments) |

### Time and memory of a full run

| Step | Models | ERA5 |
|---|---|---|
| `pr` or one of the seven `mask_*` extractions | 4.4-6.5 min, 2.4-3.1 GB | 15.8-16.7 min, 7.7-10.3 GB |
| `etc_cof`, `link_env` | seconds | seconds |
| `combine` | 1.9-2.7 min, 50-75 GB | 5.9 min, 195 GB |
| `composites` | 0.8-1.3 min, 22-47 GB | 2.2 min, 61 GB |
| `stats` | 0.3-0.5 min, at most 25 GB | 0.8 min, 19 GB |

The extraction phase (60 steps, at most 44 CPU slots in use) took 19.3 min because the steps of all six sources ran side by side; the whole core tier takes about 30 min on one node (the dry run's estimate is 27 min).

### Second round (2026-09-21, node nid004147; test area `tmp/etc_radius_fix`): radius fixed, polar ETCs removed, why the ICON and SCREAM composites changed

**Radius fixed (2d2b4e5).** `create_circular_mask(x, y, radius_deg, lon_res, lat_res)` converts the grid-point offsets with the store's global attributes (all 12 combined stores, March and test run: 0.25, radius 20); the circle is a circle in degrees (flat, like the ones plotted in the composite notebooks). Checks:
with `--basic-radius 2.5 --mask-radii 2.5 3.75 --pr-radius 2.5` the new code reproduces the March statistics of production `etc_data/stats` (replaced in the third round) exactly (SCREAM 96 variables, ICON 81, all values and coordinates equal), so the change is only the geometry; at 60 sampled points of SCREAM the statistics equal a numpy calculation
inside the 40-point (10°, 5025 pixels) and 60-point (15°, 11289 pixels) circles to the last bit; `tests/test_etc_stats_radius.py`. Every statistics file made before this (the March production files, replaced in the third round, and the figures made from them) has radii that are 4 times too small.

**Latitude rule (234e04b).** The ETC tracks are global (the extraction keeps |latitude| <= 70°) but the COF products exist only equatorward of 60°, and nothing removed the polar ETCs: `create_etc_composites.py` selected only `cof_lat` > 20 or < -20, `calc_etc_spatial_stats.py` nothing. Of the points the composites select, **20-26% of the northern and 35-44% of the southern lie poleward of 60°**.
An ETC is now used only if at least `--min-lat-coverage` of its ±20° box (161 rows) lies within |latitude| <= `--lat-limit` (60°); the equivalent centre latitude is 40° for coverage 1.0, 48° for 0.8 (composites, default until the third round below lowered it to 0.7, centre within 52°) and 60° for 0.5 (statistics, default: your centroid rule). Share of the hemisphere points that are kept (six sources):

| minimum coverage | centre within | northern points kept | southern points kept |
|---|---|---|---|
| 0.5 | 60° | 75-80% | 57-65% |
| 0.8 (composites) | 48° | 32-43% | 26-33% |
| 1.0 | 40° | 16-23% | 13-17% |

Checked: the closed-form coverage reproduces the composites' point counts (SCREAM NH 1459 and SH 1211 for `3way`, the composite `pr` equals numpy's mean to 3e-6 mm/h) and the statistics' (15306 of 22301); `tests/test_etc_domain.py`. Sample sizes at the defaults (statistics: 67-72% of the points are kept, ERA5 47886 of 66789, SCREAM 15306 of 22301, ICON 15432 of 22594, NICAM 15305 of 21563, UM 16263 of 23287, CASESM2 11860 of 17593;
composites: SCREAM `all` NH 4513 and SH 3213 points, ICON 4745 and 3049, NICAM 3977 and 2985, UM 4378 and 4001, CASESM2 2444 and 2716, ERA5 4506 and 4172). **The polar points are mostly "isolated"** (there are no MCS masks to find an overlap with): the isolated ETC points of the statistics fall by 42-59% (ERA5 31937 to 16283, SCREAM 8071 to 3304,
ICON 6708 to 2826, NICAM 10463 to 5361, UM 13151 to 7284, CASESM2 11671 to 6756), ETC-AR by 20-47%, ETC-MCS by 6-22%, ETC-AR-MCS by 4-9%. Small samples remain: northern tracks of the type ETC-AR (dominant type of a track): ICON 5, SCREAM 6, NICAM 13, UM 48, ERA5 47, CASESM2 40; northern `ar_only` composites: SCREAM 28 points, ICON 35, ERA5 87, NICAM 100.

**Why the ICON and SCREAM composites changed** (measured with the test-run stores, the March stores and the COF store, mean over points; window = the 6 hours that `pr` averages):

| ICON, distance from the centre | old (March) | Step 1 `tot_pr` of T-6 h | new (T) | `tot_pr` of T+6 h |
|---|---|---|---|---|
| 0-2.5° | 0.552 | 0.561 | 0.318 | 0.149 |
| 0-5° | 0.387 | 0.393 | 0.280 | 0.163 |
| 5-10° | 0.175 | 0.176 | 0.173 | 0.156 |
| 10-20° | 0.099 | 0.099 | 0.107 | 0.116 |

The old ICON `pr` (catalog `PT6H` mean, end-stamped) equals the Step 1 `tot_pr` of the previous window to within 1.5% at every distance, so it is the same product one window earlier, not another product. The storms move (ICON: median six-hour displacement of the centre 2.2°, 244 km, mean eastward component 1.7°) and the rain maximum lies about 5° ahead of the centre, so the window that
starts at T puts the rain further downstream: along the track the bins at -3° and +1° fall from 0.225 and 0.209 mm/h (T-6 h) to 0.128 and 0.161 (T) and 0.096 and 0.099 (T+6 h), the maximum moves from +5° to +5°/+9°, the mean within 5° falls by 29%, that at 10-20° rises by 8%, and the box mean changes by -1.4%. The same sequence appears in UM, whose old and new windows are identical
(0-2.5°: 0.570, 0.355, 0.145 mm/h for T-6 h, T, T+6 h). So the position of the window relative to the track time matters at the level of tens of percent near the centre: every source now describes the 6 hours after the track time (about 1° downstream of the instantaneous fields on average); a window centred on T would need another product (not done).

MCS precipitation fraction (rain under the exclusive MCS mask of the 3-way composites, northern points, within 10°): ICON 16.8% (March) to 25.8% (new); with the old window (`tot_pr` of T-6 h) and the new masks it is 16.7%, so the window explains all of it; SCREAM 17.2% to 23.2%, of which aligning the window (the aligned liquid file with the new masks: 22.2%) explains 4.5 points, the more intense hourly stream
(rain inside the masks 0.605 to 0.724 mm/h, snow outside) about 1.0 and the new masks 0.5; NICAM (window unchanged) 24.3% to 23.0% and UM (`pr` identical) 15.1% to 15.6%. The MCS mask is where an MCS was during [T, T+6 h), and MCSs and their rain move 2-3° in 6 hours, about their own size, so rain of another window lies only partly under the mask: for the old ICON window the rain inside the mask was
0.424 mm/h against 0.375 outside (1.13 times), aligned 0.566 against 0.301 (1.88 times). SCREAM's mean `pr` is 35% higher than in March (+82% at 50-70°S, +3% in the tropics: snow, which the old liquid-only file lacked, and the window).

**Re-run with both fixes** (composites at 0.8, statistics with the centroid rule, radii 10° and 15°; six sources, 30 jobs, 9.3 min; 14 notebooks without an error; outputs and figures in `tmp/etc_radius_fix`). Effects, in the order radius, latitude rule, data:

| Effect | Statistics | Figures |
|---|---|---|
| Radius (production 2.5° / 3.75°, all points, to radius fixed) | `pr_mean` -29% (ERA5), -33% (CASESM2), -50% (UM, NICAM), -59% (SCREAM, ICON); `omega_850hPa_min` 59-85% and `omega_500hPa_min` 49-104% more negative; `hus_850hPa_mean` -15% to -19%; `psl_min` -0.1%; MCS fraction (10°) -8% to -33%, AR fraction +41% to +126%; 54 of 108 entries of the northern significance table change by more than 5 points or in significance | multisource figures: 5-14% of the pixels change (mean absolute difference 8-18 of 255) |
| Latitude rule (all points to centroid within 60°), sample means | `pr_mean` +11% (ICON) to +30% (SCREAM), `hus_850hPa_mean` +15% to +18%, `omega_850hPa_min` 13-15% more negative (SCREAM, ICON, NICAM, UM; ERA5 6%, CASESM2 2.5%), `psl_min` +0.5-0.8%; the composites (0.8 against unfiltered): `all` mask frequencies mostly +23% to +54% (ICON NH -0.5% to +9%), `isolated` composites change most (pattern correlation of `psl` 0.17-0.75), `3way` within a few percent for `pr` and `psl` (`pr_etc` -35% to +9%, `etc_freq` -7% to -20%) | multipanel: 41-51% of the pixels of each source's panels change (mean of the 10 figures); 1-source ERA5 (radius and rule together): all 9 figures change (4-39% of the pixels) |
| Data (March to test run, both fixed) | ERA5, UM, CASESM2: `pr_mean` identical; NICAM +2.2%; SCREAM +15.7% (`pr_max` +52.6%, `pr_mcs_mean` +28.8%, `pr_ar_mean` -21.5%); ICON -10.2% (`pr_mcs_mean` +34.9%, `pr_ar_mean` -24.0%, `pr_etc_mean` -32.9%); environment variables identical; significance table: 4 of 108 entries change in the north (Max Rain Rate of ICON and SCREAM), 6 in the south; composites (`3way` NH/SH): SCREAM `pr` +10.6%/+10.1%, `pr_mcs` +33.1%/+31.1%, `pr_ar` -18.1%/-10.1%; ICON `pr` +1.8%/+0.6%, `pr_mcs` +37.6%/+29.0%, `pr_ar` -27.2%/-21.4%, `pr_etc` -43.1%/-55.2%; mask frequencies +0.2% to +7% from the round-2 masks | multisource: at most 2.4% of the pixels change; multipanel, share of pixels changed in the panels (mean of 10 figures): ERA5 3.9%, CASESM2 5.2%, UM 7.5%, NICAM 12.1%, SCREAM 19.0%, ICON 21.3%; 1-source ERA5: 7 of 9 identical, the others 0.1% and 1.4% |

The radius and the latitude rule change the results far more than the data do: they are corrections of the analysis, the data changes (window of `pr` for SCREAM and ICON, the round-2 masks) come on top. Outputs: `tmp/etc_radius_fix/{march,new}/{composites,stats}` (six sources each), `march_radius_only/stats`, `figures/{multisource,multipanel,onesource,compare}`, executed notebooks in `notebooks/executed`.

### Third round (2026-09-21, afternoon, node nid004232; test area `tmp/allfixes_20260921`): coverage 0.7, notebook labels and counts per year, Analysis 4 regenerated, Analyses 3 and 4 promoted

**Coverage of the composites lowered to 0.7 (user decision).** The default of `create_etc_composites.py` (and of `create_composites`, `in_cof_domain`, `store_in_cof_domain`) is 0.7: the centre within 52°. The statistics keep 0.5 (centre within 60°). Sample sizes for coverage 0.50 to 0.80 in steps of 0.05, with the composites' own selection
(hemisphere = `cof_lat` beyond ±20°; ERA5 with the composite script's time filter 2019-08-01 to 2020-09-01; full tables in `lat_coverage_sample_sizes.txt`, validated against `n_points` of the composite script for ERA5 and SCREAM at 0.5, 0.7 and 0.8):

| coverage | centre within | northern points kept | southern points kept | smallest `ar_only` north / south | smallest `3way` north / south |
|---|---|---|---|---|---|
| 0.50 | 60° | 75-80% | 57-69% | 162 (ICON) / 256 (ICON) | 629 / 748 (CASESM2) |
| 0.55 | 58° | 66-74% | 48-61% | 132 / 208 (ICON) | 598 / 721 |
| 0.60 | 56° | 58-69% | 42-53% | 108 / 169 (ICON) | 565 / 704 |
| 0.65 | 54° | 51-62% | 37-47% | 82 / 127 (ICON) | 532 / 670 |
| **0.70 (composites default)** | **52°** | **44-55%** | **33-41%** | **58 / 99 (ICON)** | **489 / 647** |
| 0.75 | 50° | 38-49% | 30-37% | 45 (ICON) / 82 (SCREAM) | 446 / 615 |
| 0.80 | 48° | 32-43% | 26-33% | 28 / 60 (SCREAM) | 392 / 562 |

At 0.7 the `all` composites have 3378-6121 northern and 3529-5469 southern points per source (0.8: 2444-4745 and 2716-4172). The second-round table above counted all ERA5 points, which gives 65% of the southern points at 0.5; with the composite script's time filter for ERA5 it is 69%. The tracks of the statistics notebooks (dominant type, northern ETC-AR) hardly depend on the coverage: ICON 5-6 and SCREAM 5-6 up to 0.7, NICAM 13 at 0.5 and 8 at 0.7,
so the small ETC-AR samples of ICON, SCREAM and NICAM are a property of the data, not of the rule (UM 48 to 26, ERA5 47 to 29, CASESM2 40 to 25 from 0.5 to 0.8). The composites of the test run were made again at 0.7 (`tests/test_etc_domain.py` checks the defaults); the statistics are those of the second round (radii 10° and 15°, centroid rule).

**Notebook labels and counts per year.** In the six notebooks `plot_etc_composites`, `_diff`, `_multipanel`, `plot_etc_spatialmean_stats_1source`, `_multisource` and `plot_mcs_cof_trackstats_multisource` the COF labels are joined by "-" (ETC-MCS, ETC-AR-MCS, MCS-AR-ETC; the table labels A-M and E-A). In the three notebooks with sample sizes in their figures (1-source, multisource, MCS-COF track statistics) the legends of the box plots and lifecycle plots
show the **average count per year**: the total count divided by the years covered, `(last - first time stamp) / 365.25 days` of the classified points, which does not depend on where the record crosses a calendar year (SCREAM 2019-08-01 to 2020-08-31 is 1.09 years, not 2). Years per source in the multisource notebooks: ERA5 / OBS 3.00, SCREAM 1.09, UM 1.08, ICON 1.00, NICAM 1.00, CASESM2 1.00;
the 1-source notebook filters ERA5 to 2019-08-01 to 2020-09-01 in its own cell, so its ERA5 run covers 1.09 years. The legend title says what the numbers are ("Average tracks per year (Isolated | ETC-MCS | ETC-AR | ETC-AR-MCS)"); the printed summary tables and the significance tests use the total counts. `tests/test_notebook_labels_counts.py` (labels, the helpers `estimate_n_years` and `format_per_year`, legends on synthetic data).
The string "3=MCS+AR" remains in the description attribute of `overlap_flag` written by `combine_etc_2d_vars.py` (it shows in the dataset display of the 1-source notebook) and in some scripts and READMEs; no notebook label or logic uses it.

**Analysis 4 regenerated from the round-2 Step 3 stores.** The production files of Analysis 4 (2026-03-24) came from the older Step 3 stores. `extract_mcs_cof_tracks.py` (six sources, 121-202 s each, side by side) and `combine_mcs_cof_trackstats_multisource.py` were run again with the round-2 `cof_masks` (the MCS track statistics of SCREAM, ICON, UM and CASESM2, purged from pscratch, were copied back from CFS with rsync). The scripts are unchanged.
The number of tracks equals the MCS track-statistics files for all six sources (IMERGv7 112854, SCREAM 56487, ICON 43577, UM 39340, NICAM 44295, CASESM2 13598), and the combined parquet has the columns of the March file (77, same order) and 1,951,841 rows. Only the COF flags change (the Step 3 masks), so the type shares move a little (share of the tracks by dominant COF type, in %, new / March file):

| dataset | Isolated | MCS-AR | MCS-ETC | MCS-AR-ETC |
|---|---|---|---|---|
| OBS | 50.0 / 50.3 | 14.6 / 14.5 | 19.5 / 19.4 | 15.9 / 15.8 |
| SCREAM | 55.5 / 58.0 | 11.5 / 10.8 | 20.8 / 19.7 | 12.2 / 11.5 |
| ICON | 52.5 / 53.7 | 16.5 / 16.0 | 16.9 / 16.5 | 14.1 / 13.7 |
| UM | 56.2 / 58.2 | 16.7 / 15.9 | 15.0 / 14.4 | 12.0 / 11.5 |
| NICAM | 51.1 / 52.4 | 16.7 / 16.4 | 17.8 / 17.2 | 14.4 / 14.0 |
| CASESM2 | 47.6 / 48.9 | 20.2 / 19.8 | 18.2 / 17.7 | 14.0 / 13.5 |

The MCS-COF notebook's legend numbers of the ocean box plot were recomputed independently from the parquet and match exactly. The markdown of that notebook still says "~1.48 M rows, 74 columns" (stale before this change; the file has 1.95 M rows and 77 columns).

**Analyses 3 and 4 promoted to production (2026-09-21 15:50; user approval of the figures).** 145 items, 263.0 GB, were renamed from the test areas into production: 6 `etc_tracks/<source>_etc_cof_data.parquet`, 6 combined stores `etc_data/<source>/etc_2d_combined_all_all.zarr` (251.5 GB, from `tmp/etc_round1`), the `pr` store and the 7 overlap-mask stores of each source in `etc_data/<source>/single_vars/` (48 stores, 10.9 GB),
6 statistics files `etc_data/stats/etc_spatial_stats_<source>.nc` (from `tmp/etc_radius_fix/new/stats`), 60 composites `etc_data/stats/<source>/etc_2d_composite_*.nc` and the 19 Analysis 4 files in `cof_masks/stats/` (18 per-source files and `mcs_cof_trackstats_allsources.parquet`; from `tmp/allfixes_20260921/{composites,a4}`). The 145 previous items (266.7 GB, the March data restored from CFS on 2026-09-20) are in
`/pscratch/sd/w/wcmca1/hackathon/_prev_production_A3A4_20260921/` (README, `MANIFEST_pscratch.tsv` with the 290 renames, `SIZES_pscratch.tsv`, `NEW_ITEMS_for_globus_A3A4.txt`); the CFS counterparts of the same 145 items (the same 266.7 GB and 8547 files) moved into `/global/cfs/cdirs/wcm_shr/hk25/_prev_production_A3A4_20260921/`, so that a Globus copy of the new items mirrors production.
The environment stores (linked into the test run by `link_env`), the track text files, the other files of the folders and everything else were not touched. Verified after the move: sizes of the 145 production and the 145 archived items equal the sizes measured before, no symlink inside an item, the exact chunk-key scan finds all 54 zarr stores complete (633 arrays), the point counts of the combined and single-variable stores are 66789, 22301, 22594, 21563, 23287, 17593
(ERA5, SCREAM, ICON, NICAM, UM, CASESM2), the statistics carry `radius_units = degrees`, `min_lat_coverage` 0.5 and `lat_limit` 60, the composites 0.7 and 60, the Analysis 4 files have the expected track counts, and the 120 links of the test root to the environment stores still resolve. The six repository notebooks were then executed on the production paths (saved figures: 23 of 23, embedded figures of the 1-source notebook: 9 of 9,
and the three HTML tables identical to the ones approved from the test-area run). Restore: `python ~/.claude_tmp_verify/round4/promote_a3a4.py --root /pscratch/sd/w/wcmca1/hackathon --archive _prev_production_A3A4_20260921 --rollback` (CFS: `--root /global/cfs/cdirs/wcm_shr/hk25 --archive-only --rollback`, possible while the destination names are free).
The repository notebooks now hold the outputs of that execution (the figures embedded as before; the figure folder of the execution was a scratch folder, so the CFS figure folders of the notebooks are unchanged) and the missing line break of the 1-source notebook (`/ 100track_stats.head()`) is fixed.

## GSMaP added to the extreme-precipitation-threshold analysis; per-year percentiles and interannual IQR (2026-09-22, `8f94ff6`)

**GSMaP (GSMaPv8) added as a second observational product**, for cross-checking IMERG in `plot_extreme_rain_threshold_map.ipynb`'s 7-source comparison: a new `GSMAP` entry in `config_sources.yaml` and a dedicated loader branch in `calc_extreme_precip_thresholds.py` that reads a local 6-hourly HEALPix store directly (`/pscratch/sd/w/wcmca1/GSMaP/healpix/GSMaPv8_6H_zoom8_20100101_20241231.zarr`). Not part of the COF pipeline registry (`config_pipeline.yaml`) — no MCS/AR/ETC tracking for GSMaP, so it is run standalone through `calc_extreme_precip_thresholds.py` only, over the same 2019-2021 window IMERG's 6-hourly store is limited to. Checked for the Arctic-ocean retrieval artifact the 1-hourly version of this analysis found and masked (GSMaP over the Arctic ocean, P95 > 6 mm/h, lat > 75, landfrac < 0.5): not present in the 6-hourly (temporally averaged) data — those cells are lower, not higher, than the rest of the domain (P95 mean 1.05 mm/h there vs. 2.74 mm/h domain-wide) — so no masking was added.

**Per-calendar-year percentiles and their interannual IQR**, ported from `calc_extreme_precip_thresholds_1h.py`: for any source with at least 2 qualifying calendar years (`--min_year_coverage_days`, default 300 distinct days; `--no_annual` turns it off), the threshold file gains a `year` coordinate plus `pr_annual_p90`/`pr_annual_p95` and their interannual `pr_q25_p*`/`pr_q75_p*`/`pr_iqr_p*`. Purely additive — IMERG's existing `pr_p90`/`pr_p95` verified bit-identical before and after the port, and a `scream` spot check (1 qualifying year) confirmed it gains no new variables, as expected. New test file `tests/test_extreme_precip_annual_iqr.py` (7 tests, synthetic gamma-distributed data): coverage gate, threshold behavior, IQR-vs-`numpy.quantile` match, `write_netcdf` with/without annual data.

**11-year (2014-2024) vs. 3-year (2019-2021) threshold comparison, at the user's request** — built an 11-year IMERG 6-hourly store with `coarsen_healpix.py` (1-hourly → 6-hourly, same-zoom temporal-only resample of `IMERG_V7_1H_zoom8_20010101_20241231.zarr`, ~69 s) and computed both IMERG and GSMaP thresholds over 2014-2024. The interannual IQR is meaningfully larger at 11 years than 3 (domain-mean `pr_iqr_p95`: IMERG 0.485 → 0.668 mm/h, GSMaP 0.642 → 0.870 mm/h, both ~1.36-1.38x) and the shaded uncertainty band goes from barely visible to clearly visible in the zonal-mean panel — but the P90/P95 maps themselves barely move (~1-2% domain-mean difference). **Decision (user, 2026-09-22): keep the 3-year threshold for production; do not re-run the storm-type attribution (`calc_stormtype_extreme_precip_spatial.py`) against the 11-year data.** The 11-year files are kept as archival references only (`{source}_precip_percentiles_6h_hp8_v1_2014_2024.nc`, both `IMERGv7` and `GSMaPv8`, on pscratch and CFS), outside the pipeline registry and not read by `calc_stormtype_extreme_precip_spatial.py` or `run_cof_pipeline.py`.

**Notebook restructuring**, following from that decision: production's IMERG threshold file predates this port and carries no annual/IQR variables, so `plot_extreme_rain_threshold_map.ipynb`'s zonal-mean panels gained a `shade_interannual_iqr` toggle (default `False`) and a guard (`years_dict[key] is not None and annual_var in ds.data_vars`) instead of assuming every source has annual data. The interannual mean+IQR figure itself moved out of that notebook into a new, dedicated one, `plot_extreme_rain_threshold_map_obs_interannual.ipynb` (IMERG + GSMaP only, no GSRM/models, reading the 11-year archival files directly) — keeping the "does production have annual data" question fully decoupled from "is the interannual comparison available at all". `plot_maps_grid` (duplicated in both notebooks, per this repo's existing per-notebook-copy convention) also lost a hidden `wspace`/`hspace` floor of 0.12 that silently overrode anything passed below it whenever a zonal panel was present, and gained a `zonal_height_frac` parameter (mirrors the existing `zonal_width_frac`) as an additional, independent clearance lever.

**Production**: `GSMaPv8_precip_percentiles_6h_hp8_v1.nc` (3-year) added to `extreme_precip/` on pscratch and CFS; the two 11-year archival files added likewise; `IMERGv7_precip_percentiles_6h_hp8_v1.nc` (the plain 3-year file) left untouched throughout. All copies verified (size, variable set, `xr.open_dataset` succeeds); both notebooks re-run end-to-end against their own committed (CFS) `rootdir` with no errors, producing figures byte-identical to the test-area runs already reviewed and approved.

## Follow-ups after this phase

Held on purpose while the pipeline runner and the second full round are done; none of them changes the results of Analyses 1 and 2 as they stand.

1. **`extract_etc_2d_vars.py` and the SCREAM 6-hourly file.** *Done 2026-09-20 (Analysis 3 section above).* The extraction no longer reads `scream_pr6h_z8.zarr`, the aligned file or the catalog's 6-hour mean: `pr` of the models is Step 1's `tot_pr`
   (`--pr_source cof_tot_pr`, the default); `--pr_source legacy` restores the old files. `scream_pr6h_z8_aligned.zarr` is not read by anything.
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
7. **IFS** (`ifs_tco3999_rcbmf`) is in `config_sources.yaml` and has an MCS tracking config (`config_mcs_tbpf_ifs_tco3999_rcbmf.yml`, confirmed 2026-09-22 to be the only IFS tracking config that exists) — no AR, TC or ETC tracking is configured for it, so it cannot go through the co-occurrence pipeline and is not enabled in `config/config_pipeline.yaml`.
8. **Step 3 exit status.** *Done 2026-09-20.* `make_cooccurrence_masks.py` exited with status 0 when a frame had no result (the frame stayed NaN and only a log line said so), when its input could not be read and when its store could not
   be written. It now lists the frames without a result and exits 1 in all three cases (`stream_process_to_zarr(return_missing=True)`, `all_time_steps_written`; `tests/test_step3_fail_loud.py`); UM on 48 real frames is bit-identical
   to production and exits 0. All six round-2 Step 3 runs processed every frame, so nothing that ran before would fail. Found on the way and not changed: the sequential fallback of `stream_process_to_zarr` sits inside the
   `if parallel and client` branch, so `--no-parallel` never processes a frame (the store stays all-NaN); it now fails with exit status 1 instead of exiting 0.
9. **Analysis 3: radii of the spatial statistics.** *Done 2026-09-21 (2d2b4e5).* The radii are in degrees now (x and y converted with `lon_res` and `lat_res`); statistics made before that (the March production files, now in `_prev_production_A3A4_20260921`) use 2.5° and 3.75°.
10. **Analysis 3: the `all` composites include points without a frame.** Not changed, needs a decision: their mask frequencies and exclusive precipitation are low by the share of those points (ICON 16.7% NH / 12.6% SH; UM 3.7% / 2.8%; SCREAM 0.6% / 0.5%);
    a `has_frame` selection of about three lines in `create_composites` fixes it (Analysis 3 section, "Not changed on 2026-09-20 and what became of it").
11. **Analysis 3 outputs in test areas.** *Done 2026-09-21 (third round above).* The outputs of the core-tier run (`tmp/etc_round1`), the statistics of the radius and latitude re-run (`tmp/etc_radius_fix/new/stats`) and the composites at coverage 0.7 (`tmp/allfixes_20260921/composites`) were promoted by rename with the regenerated Analysis 4 files (145 items, 263.0 GB);
    the previous production (the March data) is in `_prev_production_A3A4_20260921` on pscratch and CFS. The test areas keep the run markers, logs and figures, the March-store re-run (`tmp/etc_radius_fix/march*`) and the symlinks to the environment stores. The runner protects the production outputs (they have no markers, so a run with `--data-root` set to the production tree stops unless `--force`).
12. **Analysis 3: polar ETCs.** *Done 2026-09-21 (234e04b; default of the composites 0.7 in the third round).* Composites and statistics use only the ETCs whose box lies within |latitude| <= 60° for at least `--min-lat-coverage` of its rows (0.7 for the composites, 0.5 for the statistics; adjustable). The small samples that remain
    (northern ETC-AR: ICON 5, SCREAM 6, NICAM 13 tracks; northern `ar_only` composites of 58-197 points for ICON, SCREAM, NICAM and ERA5 at 0.7) are a property of the data (third round above). Composites and statistics made before it contain the polar ETCs.
13. **Analysis 3: the window of `pr` relative to the track time.** Not changed, needs a decision: `pr` of every source is the mean of [T, T+6 h), on average about 1° downstream of the instantaneous fields at T; near the centre the composites depend on this at the level of tens of percent (second round above).
    A window centred on T would need hourly precipitation (the inputs of Step 1), not the 6-hourly `tot_pr`.
14. **Globus copy of the round-2 items.** *Done (user, 2026-09-21), verified.* All 54 items of `NEW_ITEMS_for_globus.txt` exist on CFS with the same number of files and the same byte size of every file as on pscratch (192.6 GB each).
15. **Analysis 4 regenerated and promoted.** *Done 2026-09-21 (third round above).* The per-source files and the combined parquet come from the round-2 Step 3 stores (the scripts are unchanged); the type shares differ from the March files by up to 2.5 points (Isolated).
16. **Notebook labels and sample sizes per year.** *Done 2026-09-21.* The COF labels of six notebooks are joined by "-"; the box plots and lifecycle plots of the three notebooks with sample sizes show the average count per year. "MCS+AR" remains in the description attribute of `overlap_flag` (`combine_etc_2d_vars.py`) and in scripts and READMEs; the stale markdown of the MCS-COF notebook (1.48 M rows, 74 columns) is unchanged.
17. **Globus copy of the Analysis 3 and 4 items.** *Done (user, 2026-09-22), verified.* The 145 new items (263.0 GB) listed in `NEW_ITEMS_for_globus_A3A4.txt` in the pscratch archive folder were copied to CFS. Verified 2026-09-22 by comparing every one of the 145 items between pscratch and CFS (file count and total byte size per item): 8,490 files, 263.01 GB, identical on both sides, 0 mismatches or missing items.
18. **Analysis 2, IMERG 11-year vs. 3-year threshold.** *Decided 2026-09-22 (GSMaP section above).* Kept the 3-year threshold for production; the interannual IQR is measurably larger at 11 years (~1.36-1.38x domain-mean) but the P90/P95 maps themselves change only ~1-2%, and the storm-type attribution was not re-run. The 11-year IMERG and GSMaP threshold files are kept as archival references (`_2014_2024` suffix), outside the pipeline registry, for a future revisit of this decision.
19. **`plot_extreme_rain_threshold_map.ipynb`'s `shade_interannual_iqr` toggle.** *Off by default (GSMaP section above), needs a decision only if revisited.* Production's IMERG threshold file has no annual/IQR variables, so the toggle is off and the zonal panel shows plain lines for every source. Flipping it on needs a source in `rootdir` whose file has `pr_annual_p*`/`year` (GSMaP already qualifies; IMERG would need to be swapped for the 3- or 11-year version this round regenerated, which was deliberately not done to production).
