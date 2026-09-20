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
2. **Extreme thresholds.** `calc_extreme_precip_thresholds.py` still uses each source's 6-hourly product while the attribution uses `tot_pr`.
   One-month comparison (thresholds from `tot_pr` / from the product): UM, IMERG and CASESM2 identical (median ratio 1.000); NICAM same mean, per-cell
   correlation 0.94; ICON 0.56-0.88 at high and northern mid-latitudes (snow counted twice in the product); SCREAM 1.21-1.25 at the median. A
   full-record comparison decides; `compare_thresholds.py` (audit tooling, `thresholds/`) writes to a new folder and never overwrites.
3. **SCREAM 6-hourly product** (`scripts/avg_pr_healpix.py` -> `scream_pr6h_z8.zarr`, still read by the threshold script and by
   `extract_environments/extract_etc_2d_vars.py`):
   - It is built from `scream_ne120`, 3-hourly averages that are labelled at the end of each interval (the first label is 2019-08-01 03:00).
     `resample(time='6h', label='left', closed='left')` therefore puts the means labelled T and T+3 h in the window labelled T, which covers
     [T-3 h, T+3 h) instead of [T, T+6 h). With `closed='right'` the correlation with Step 1's hourly window mean rises from 0.74 to 0.87
     (16 windows from 2020-07-11, |lat| < 20). The window centred on T can suit pairing with an instantaneous track time, but not Step 1's windows.
   - It holds liquid only: the `prs` line in `vars_to_include` is commented out (the 3-hourly `pr` is `precip_liq_surf_mass_flux`).
   - The noleap calendar mapped to standard dates leaves four all-NaN frames on 2020-02-29; the first (2019-08-01T00) and last (2020-09-01T00)
     windows hold a single 3-hourly mean; `resample().mean()` does not check that a window is complete.
   - It is a different data stream from Step 1's: 3-hourly averages saved at 25 km and interpolated to HEALPix level 8, against hourly
     instantaneous fields remapped to level 10 and coarse-grained to level 8. The mean rate is the same, but in the tropics the 3-hourly field has
     20% wet cells (at least 0.1 mm/h) against 16%, and P99 of 6.8 against 9.8 mm/h. This, not the time alignment (which leaves the distribution
     unchanged) and not sampling (two instantaneous 3-hourly samples per window add 13-19%), accounts for most of the +23-37% tropical threshold difference.
4. **IMERG.** Its AR, TC and ETC masks come from ERA5, so its Step 2 is `combine_era5_imerg_tracking_masks.py` (no command-line options; it carries
   `tot_pr`); `slurm/tasks_all_IR_IMERG.txt` lists only Step 3 and the monthly script.
5. **Partial windows.** Windows with fewer than six hourly steps are written from the steps present, including single-step windows at gaps and record
   ends (SCREAM 2020-04-20T00 and 2020-09-01T00, CASESM2 2021-03-01T00). No minimum number of steps is enforced.
