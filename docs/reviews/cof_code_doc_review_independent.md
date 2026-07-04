# Independent Review of the Codex COF Code-Documentation Review

**Author:** Claude, reviewed and tested by Zhe Feng

This report is an independent, code-driven verification of
[cof_code_doc_review.md](cof_code_doc_review.md) (the "Codex report"). I re-read the four scripts,
four procedure docs, and both plotting notebooks
named in that report directly, traced every cited line, and re-derived each conclusion rather than
trusting the citations. Where I agree, I say so briefly and add line references of my own. Where I
found nuance, additional context, or something the Codex report missed, I call it out explicitly.

The first version of this report (below, unchanged) was a review-only artifact that did not modify
any code or documentation. After user review (see **Revision note** below), every finding and item
in this report was fixed in the code/docs/notebooks; each section now also records what was
actually changed, under a **"Fix applied"** heading.

**Scope reviewed:** `scripts/make_mcs_swath_masks.py`, `scripts/make_cooccurrence_masks.py`,
`scripts/calc_monthly_rainmap_by_cof.py`, `scripts/calc_stormtype_extreme_precip_spatial.py`, their
four corresponding procedure docs, and `notebooks/plot_cof_extreme_raintype_rank_map.ipynb` /
`notebooks/plot_cof_raintype_rank_map.ipynb`.

**Revision note:** this report was revised after user feedback (Zhe Feng, 2026-07-02). The user
corrected one claim in Finding 4, elevated Finding 5's priority based on an observed real-world
anomaly, resolved Items B and C, and requested an explicit/defensive AR/ETC TC-overlap filter be
added. All findings and items in this report (2, 4, 5, 8, A, B, C, and the new Item D) have now
been **fixed in the code/docs/notebooks**, as noted inline and in the priority table below. Findings
3, 6, 7 required no further action (already correct from the previous pass).

**Second revision note (Zhe Feng, 2026-07-04):** after the fixes above were applied, the user reran
the full pipeline and the four affected plotting notebooks and confirmed the fixes' actual impact on
the scientific results — see the new **"Post-fix validation"** section below. The paper-methodology
paragraphs and the "Notes on how this review was conducted" section were also updated to reflect the
current, fixed state of the code and the completed rerun.

---

## Summary verdict

I agree with all 8 Codex findings — every one reproduces from the current code. On three of them
(Findings 1, 4, 5) I have a materially different read of severity or completeness, detailed below.
I also identify three additional low-severity items the Codex report did not mention. **All items
below are now fixed** (see the per-finding sections and the consolidated priority table for what
changed and where); severity below reflects the review's original assessment, revised per user
feedback.

| # | Codex finding | My verdict | Severity (mine, revised) | Status |
|---|---|---|---|---|
| 1 | No-MCS 6-h window crashes swath processing | Agree — real bug, but practically rare for global data | Low–Medium (Codex: High) | Fixed |
| 2 | Extreme script swaps ST and ND cloud types | Agree — fully confirmed, highest-impact item | **P0** | Fixed |
| 3 | Monthly `*_count` are time-step counts, not hours | Agree — already fixed in docs/metadata | Low (resolved) | Done (prior pass) |
| 4 | Monthly notebook has a residual formula hazard | Agree — blast radius larger than stated, and a second, worse bug was found during implementation | Medium–High | Fixed |
| 5 | Cloud-type TC filtering is class-level, not object-level | Agree; user connects this to an observed CASESM2 anomaly and elevated it regardless of magnitude | **P0** (user-elevated) | Fixed |
| 6 | COF threshold metadata was misleading | Agree — already fixed | Low (resolved) | Done (prior pass) |
| 7 | Extreme doc claimed unsaved `{type}_precip` variables | Agree — already fixed | Low (resolved) | Done (prior pass) |
| 8 | Extreme notebook residual creates artificial 100% at no-event cells | Agree; also plausibly connected to the CASESM2 anomaly | **P0** (elevated) | Fixed |

---

## Finding-by-finding verification

### Finding 1 — No-MCS 6-hour window can crash swath processing

**Codex claim:** an all-zero `mcs_mask` in a 6-hour window makes
`create_track_swaths_and_coverage()` return `{}`, and `combine_swaths_with_priority()` then does
`list(track_swaths_dict.keys())[0]`, raising `IndexError`.

**My verification — confirmed.** In [make_mcs_swath_masks.py](../../scripts/make_mcs_swath_masks.py):
- `unique_tracks = unique_tracks[unique_tracks > 0]` (line 73) leaves `track_swaths_dict = {}` when
  no pixel is positive.
- `combine_swaths_with_priority()` immediately does `first_key = list(track_swaths_dict.keys())[0]`
  (line 114) — this is an unguarded index into an empty list and will raise `IndexError`.
- The call site is `process_timechunk_swath()`, lines 487-488.

**Where I differ from Codex: severity.** The Codex report rates this "High" with the framing "any
model/source/time period with an MCS-free aggregation window can abort the pipeline." I'd temper
this: the pipeline runs on **global** HEALPix grids (`main()` processes the whole domain each
chunk, `make_mcs_swath_masks.py:1049-1073`). A 6-hour window with literally zero MCS pixels
anywhere on Earth is very unlikely in the production GSRM/observational datasets this pipeline
targets. The realistic triggers are: `--test-steps` debug runs, regional subsets someone creates
later, or a source with unusually sparse/short MCS tracking output. I'd call this Low–Medium in
practice, not High — but it is a real crash risk and the guard is essentially free, so I agree it
should be fixed regardless of how rarely it fires.

**Fix assessment: correct and sufficient.** The proposed guard —
```python
if not mcs_swaths_dict:
    combined_mcs_swath = np.zeros(mcs_mask.shape[1:], dtype=int)
else:
    combined_mcs_swath = combine_swaths_with_priority(mcs_swaths_dict, mcs_coverage_dict)
```
correctly reproduces what the rest of the function already expects for an MCS-free cell (a swath
value of 0), and every downstream consumer of `combined_mcs_swath` in
`process_timechunk_swath()` (lines 532, 536-540, 554) already treats 0 as "no MCS," so no other
code path needs to change.

**Fix applied.** `scripts/make_mcs_swath_masks.py` now guards the empty-dict case exactly as
proposed, with a brief comment noting the MCS-free-window scenario.

---

### Finding 2 — Extreme attribution swaps stratiform and non-deep-convective cloud types

**Codex claim:** the upstream swath script defines type 2 = stratiform, type 3 = non-deep
convective, but the extreme-attribution script reverses these labels, and the swap propagates
uncorrected into the extreme plotting notebook's maps and regional summaries.

**My verification — confirmed, and this is the most important finding in either report.**

Upstream convention, [make_mcs_swath_masks.py](../../scripts/make_mcs_swath_masks.py):
```
199	    1 = Deep convective (tb < tb_thresh & pr >= pr_threshold)
200	    2 = Stratiform (tb < tb_thresh & pr < pr_threshold)
201	    3 = Non-deep convective (tb >= tb_thresh & pr >= pr_threshold)
202	    4 = Drizzle (tb >= tb_thresh & pr < pr_threshold)
```
and reaffirmed in the output metadata (`flag_meanings`, line 930) and the procedure doc
([mcs_swath_cloud_type.md:131-137](../procedures/mcs_swath_cloud_type.md)). The monthly script
respects this convention throughout (`calc_monthly_rainmap_by_cof.py` uses the precomputed `st_pr`
for type 2 and `nd_pr` for type 3 without re-deriving the classes).

The extreme script inverts it,
[calc_stormtype_extreme_precip_spatial.py:531-543](../../scripts/calc_stormtype_extreme_precip_spatial.py):
```python
# ND = 2
mask_nd = (cloud_types_cleaned == 2)   # this is stratiform upstream
...
# ST = 3
mask_st = (cloud_types_cleaned == 3)   # this is non-deep convective upstream
```
So the `nd_frac`/`nd_count` output variables actually hold **stratiform** extreme precipitation,
and `st_frac`/`st_count` actually hold **non-deep convective** extreme precipitation.

**Downstream propagation — I traced the full chain and confirm no compensating swap exists
anywhere:**
- `extreme_precip_by_stormtype.md` lines 161-162 document the *code's* (swapped) convention —
  `ND | cloud_types == 2` and `ST | cloud_types == 3` — rather than the upstream convention, so
  the doc is internally consistent with the bug, not a source of truth for what's physically
  correct.
- `plot_cof_extreme_raintype_rank_map.ipynb` reads the swapped fields verbatim:
  `nd_frac = 100 * ds.nd_frac` (line 401), `st_frac = 100 * ds.st_frac` (line 402). The notebook's
  own markdown header table (lines 40-42) and its legend `cof_labels` list (lines 1406, 1408,
  and repeated at lines 1217/2310/2709) both label these as "ND" and "ST" with no relabeling —
  the notebook author clearly believed the variables meant what their names say.
- The ranking dictionary assigns `nd_frac → code 10`, `st_frac → code 11` (lines 970-971,
  1137-1138), which feed directly into the rank/dominant-type maps and the regional/tropics/
  land-ocean averaging cells (lines 547-916) and pie-chart arrays (~2585-2608).

**Net effect:** any map, table, or summary statistic in the released notebook labeled "ND" is
actually showing stratiform precipitation, and everything labeled "ST" is actually non-deep
convective. The combined ST+ND total is correct (unaffected by the swap); only the per-category
attribution and its physical interpretation are reversed. Given this feeds a JGR paper draft, I
consider this the single highest-priority fix in either report.

**Fix assessment: correct, and appropriately conservative.** Codex's suggested fix swaps only the
comparison values (`==2`↔`==3`) inside the existing `nd`/`st` blocks, leaving the dict keys and
therefore the output variable names (`nd_frac`, `st_frac`, `nd_count`, `st_count`) unchanged. That
is the right minimal-diff choice: it fixes the *content* of the variables without touching their
*names*, so the already-written notebook code (`ds.nd_frac`, `ds.st_frac`, the `cof_labels` legend
order) keeps working correctly once the script is rerun — no notebook edits are required, only a
rerun against freshly regenerated NetCDF files. I'd add one thing the Codex report doesn't
mention: `extreme_precip_by_stormtype.md` lines 161-162 need the same correction (swap which
`cloud_types` value maps to ND vs. ST) — Codex does flag updating "lines 160-161" for this, so we
agree on the need, just off by one line in the citation (the actual table rows are 161-162 in the
current file).

**Fix applied.** `scripts/calc_stormtype_extreme_precip_spatial.py` now assigns `st` from
`cloud_types_cleaned == 2` and `nd` from `cloud_types_cleaned == 3`, matching the upstream
convention; dict keys and block order unchanged, exactly as recommended above.
`docs/procedures/extreme_precip_by_stormtype.md`'s priority table (rows 10-11) was updated to
match. **This changes the extreme-attribution NetCDF output and must be rerun before the extreme
notebook is rerun.**

---

### Finding 3 — Monthly count variables are time-step counts, not hours

**Codex claim:** `*_count` and `*_precipitation_count` in the monthly output are counts of
sub-daily time steps, not elapsed hours; this has already been corrected in the script attributes
and procedure doc.

**My verification — confirmed, and independently confirmed as already resolved.** In
[calc_monthly_rainmap_by_cof.py](../../scripts/calc_monthly_rainmap_by_cof.py), precipitation sums
multiply by `time_interval` (e.g. line 304 `chunk_totprecip = (precipitation * time_interval).sum(dim='time')`),
while counts do not (e.g. line 308 `chunk_mcs_count_sum = (mcs_mask > 0).sum(dim='time')` and line
309 `chunk_mcs_pcp_count_sum = (precipitation.where(mcs_mask > 0) > pcp_thresh).sum(dim='time')`).
The variable attributes I read (e.g. lines 821-824, `'Number of time steps MCS is present'` /
`'Number of time steps MCS precipitation exceeds threshold'`) already describe these correctly as
step counts, and `monthly_precip_by_cof.md` lines 190-201/240-241 match. Agree with Codex; no
further action needed beyond what's already been done.

---

### Finding 4 — Monthly notebook has a residual formula hazard

**Codex claim:** the notebook's first residual computation double-applies a `24/nhours` factor to
already-normalized rates, producing a near-zero, dimensionally wrong residual; a later cell
recomputes the residual correctly, and the ranking dictionary uses the later value, so a full
top-to-bottom run is not affected for the ranking map — but running cells out of order, or using
the first residual plot, shows the wrong field.

**My verification — confirmed, and I found the blast radius is larger than the Codex report
states.** In `notebooks/plot_cof_raintype_rank_map.ipynb`:
- Lines 267-291 compute each category as an already-normalized mm/day rate, e.g. line 270:
  `tot_pcp_all = (24. * ds.precipitation.sum(dim='time') / nhours).compute()`.
- Lines 293-295 then reapply the same factor:
  `res_pcp_all = (24. * (tot_pcp_all - (... + dc_pcp_all + nd_pcp_all + st_pcp_all + dz_pcp_all)) / nhours).compute()`,
  and line 318 sets `res_pcpfrac_all = 100. * res_pcp_all / tot_pcp_all` from that value — this
  drives a standalone map at line 328.
- A corrected version appears later: lines 999-1011 build
  `total_iso_pcp_all`, `total_cof_pcp_all`, `total_cloudtypes_pcp_all`, sum them into
  `total_feature_pcp_all`, and redefine `res_pcp_all = tot_pcp_all - total_iso_pcp_all -
  total_cof_pcp_all - total_cloudtypes_pcp_all` (line 1007) and
  `res_pcpfrac_all = 100. * res_pcp_all / tot_pcp_all` (line 1010). The ranking dictionary at line
  1295 (`pcpfrac_dict = {..., 'res_pcpfrac_all': res_pcpfrac_all, ...}`) executes after this
  redefinition, so in a full top-to-bottom run, the ranking map does use the corrected value, as
  Codex states.

*(Line numbers above reflect the notebook's state at the time of the original review pass; the
notebook was edited locally during this review — see the "second bug found during implementation"
note further below, where fixes are described against the current, verified content by cell ID.)*

**What I found that Codex's report doesn't mention:** the notebook's regional/latitude-band
aggregation cells — `res_frac_trop` (line 496), `res_frac_trop_l`/`res_frac_trop_o` (lines
515, 534), `res_frac_nh_l`/`res_frac_sh_l`/`res_frac_nh_o`/`res_frac_sh_o` (lines 569, 588, 607,
626), and their `.mean(dim='cell')` summaries later in the notebook — all read `res_pcpfrac_all`
by variable name, and **all of these cells sit between the first (buggy) definition at line 318
and the corrected redefinition at line 1007**. In a standard top-to-bottom run, every one of these
regional residual numbers is computed from the near-zero, dimensionally-wrong first version, not
the corrected one. So this isn't only a "you have to run cells out of order" risk as Codex frames
it — in the notebook's own natural cell order, an entire block of regional residual statistics
(tropics/NH/SH, land/ocean splits) uses the broken value even in a clean top-to-bottom execution.
Only the final ranked-map cell happens to come after the fix.

**Correction (per user feedback):** I originally wrote that "TC is not excluded from AR/ETC
footprints in this mask product," implying the residual could never be an exact closure term even
after the arithmetic bug was fixed. The user (Zhe Feng) corrected this: TC-overlap removal from AR
and ETC masks is already expected to happen **upstream of this repository**, in the AR/ETC/TC mask
generation step maintained by colleague Bryce Harrop (Bryce.Harrop@pnnl.gov) — so this was not a
live residual-math concern in the actual data, and my original caveat overstated the risk. As a
defensive measure regardless, an explicit pixel-level TC-overlap filter has now been added for AR
and ETC inside `make_cooccurrence_masks.py` (see **New Item D** below), so the mutual-exclusivity
guarantee is enforced within this repo too, not only assumed from upstream.

**A second, more severe bug found during implementation (not in the original Codex report or my
first pass):** when I re-verified this cell's exact current content immediately before editing it,
its residual formula had changed from what I originally read — the subtracted terms were no longer
the mutually-exclusive isolated/co-occurrence decomposition, but `mcs_pcp_all + ar_pcp_all +
etc_pcp_all + tc_pcp_all` (the "all instances of this feature" totals, i.e. before splitting into
isolated vs. co-occurring). Since `mcs_pcp_all` and `ar_pcp_all` both include precipitation from
MCS-AR co-occurring cells, this double-counts any co-occurring cell in the subtraction — compounding
the pre-existing double-scaling bug with a double-counting bug. This is likely the result of the
notebook having been actively edited locally during this review (a live/uncommitted file, not
something either report analyzed). I fixed the version actually on disk, not my stale notes — see
below.

**Fix applied.** In `notebooks/plot_cof_raintype_rank_map.ipynb` (cell id `5d9d75b1`), the residual
computation was replaced in place with the mutually-exclusive decomposition (isolated + co-occurrence
+ cloud types — the same ingredients the later, already-correct block at cell `119d5001` uses) and
the redundant `24./nhours` rescaling was removed, since all inputs are already in mm/day:
```python
res_pcp_all = tot_pcp_all - (
    mcs_iso_pcp_all + ar_iso_pcp_all + etc_iso_pcp_all + tc_pcp_all +
    mcs_ar_pcp_all + mcs_etc_pcp_all + ar_etc_pcp_all + mcs_ar_etc_pcp_all +
    dc_pcp_all + nd_pcp_all + st_pcp_all + dz_pcp_all)
```
No cell reordering was needed — every ingredient was already computed earlier in the same cell, so
this single in-place replacement fixes the standalone plot, all regional/tropics/land-ocean
aggregation cells, and keeps the later duplicate block (left untouched) consistent with it. Verified
structurally: only this one cell's source changed relative to the file's prior state.

---

### Finding 5 — Cloud-type TC filtering operates on whole classes, not spatial objects

**Codex claim:** `filter_mcs_tc_overlaps()` is reused for both MCS tracks and cloud types; for
cloud types it receives `_ds.cloud_types` (class labels 1-4, not object IDs), so it computes
overlap fractions per class and can remove an entire cloud class domain-wide rather than only the
TC-adjacent pixels; documentation doesn't explain this; this can inflate `unassigned` extreme
precipitation.

**My verification — the mechanism is confirmed exactly as described.** In
[make_cooccurrence_masks.py](../../scripts/make_cooccurrence_masks.py):
```
1034	    cloudtypes_filtering_results = filter_mcs_tc_overlaps(
1035	        mcs_mask=_ds.cloud_types,
1036	        tc_mask=_ds.tc_mask,
1037	        overlap_threshold=CLOUD_TYPE_TC_FILTER_THRESHOLD,
1038	        verbose=verbose
1039	    )
```
Inside `filter_mcs_tc_overlaps()`, the "track ID" extraction is generic:
`tracks_in_overlap, overlap_counts = np.unique(mcs_mask.where(overlap_condition), ...)` (lines
418-420) — when `mcs_mask` is actually `cloud_types`, this returns the small set of class values
{1,2,3,4} present near TC pixels, not per-object IDs. The removal step,
`tracks_to_remove_mask = np.isin(mcs_mask, tracks_exceeding)` (line 518), then zeroes **every**
pixel of that class in the entire time-step domain if the class's global TC-overlap fraction meets
the 1% threshold — not just the pixels physically adjacent to a TC. `cof_identification.md` lines
86-89 describes "the same procedure is applied to the cloud-type mask," implying object-level
behavior parallel to the MCS case, without flagging this class-vs-object distinction.

**Where I differ from Codex: likely impact.** I'd rate the practical impact lower than the Codex
report implies, for two reasons I verified directly:
1. The overlap fraction denominator is the **entire domain's pixel count for that class at that
   time step**. Common classes (stratiform, drizzle) span huge areas globally, so it would take an
   enormous, implausible amount of TC-adjacent area in that single class to push the domain-wide
   fraction past 1%. In the overwhelmingly common case, this filter likely removes **nothing at
   all** for cloud types (the function's own early-exit at lines 500-512, "No tracks exceed
   overlap threshold, returning original mask," is probably hit most of the time for the
   cloud-type call). When it does trigger, the effect is the coarse all-or-nothing class removal
   Codex describes — so the bug is real, but it is likely rare-but-severe rather than
   continuously present.
2. Both downstream consumers already independently exclude `tc_mask` from cloud-type accounting:
   the monthly script's `all_feature_mask` includes `(tc_mask > 0)` (line 360), and the extreme
   script assigns TC pixels to the `tc` category (priority 4, lines 500-505) **before** it reaches
   the cloud-type block (priority 5, lines 507-550) and additionally removes `tc_mask > 0` pixels
   from `cloud_types_cleaned` (line 516). So even with the class-level bug doing nothing useful (or
   occasionally doing too much), TC-cloud-type double counting is already prevented elsewhere in
   the pipeline. The main risk from this bug is specifically the rare full-class wipeout, not
   routine TC-adjacent leakage into `unassigned`.

I'd characterize this as a genuine correctness defect worth fixing (the function does not do what
its name/threshold implies), but I would not, without the targeted before/after class-area
diagnostic Codex itself proposes, assume it's a primary driver of any specific large-residual
anomaly (e.g., CASESM2 P95) — that still needs the empirical check both reports call for.

**User elevation (per feedback):** the user connects this finding directly to an observed anomaly —
CASESM2 extreme-precipitation maps show large regions where "residual"/unassigned ranks as the top
extreme-precipitation category, which the user considers physically implausible given the categories
are designed to be mutually exclusive, and given "unassigned" is not a physically coherent storm
type in its own right. The user asked for this fixed **regardless of quantified magnitude**, which I
agree with independent of my own impact estimate above: a full-class wipeout is a genuine
correctness bug even if it fires rarely, and "rarely" is not the same as "never" for a 6-hourly,
multi-year, multi-source dataset. I'll also flag a hypothesis worth checking together with this fix:
Finding 8 (below) independently produces a *spurious* 100% residual at every extreme-precipitation-free
cell in the extreme notebook, unrelated to any real feature attribution. The CASESM2 "large
residual" symptom could plausibly be a mix of both effects, and disentangling them requires both
fixes to be applied before rerunning the diagnostic Codex proposed (plot `unassigned_frac` /
`total_extreme_count` / cloud-type area removed by TC filtering, per Codex's Finding 5 recommendation).
This is my own hypothesis, not a confirmed diagnosis — I have no data access in this session to
verify it.

**Fix applied.** In `scripts/make_cooccurrence_masks.py`, Step 1 of
`process_single_timestep_overlaps()` (previously calling `filter_mcs_tc_overlaps(mcs_mask=_ds.cloud_types,
...)`) now uses a genuine pixel-level spatial exclusion:
```python
cloud_types = xr.where(tc_present, 0, _ds.cloud_types)
```
where `tc_present = _ds.tc_mask > 0`. This is applied together with the new AR/ETC filtering (see
**New Item D**) in the same code block, with an optional verbose diagnostic print of how many
pixels were removed from each field. The now-unused `CLOUD_TYPE_TC_FILTER_THRESHOLD` constant was
removed, and the `tc_filtering_threshold` output metadata attribute was updated to describe
cloud-type filtering as pixel-level rather than a percentage threshold.
**This changes the `cofmasks` output and must be rerun before downstream monthly/extreme scripts
and notebooks are rerun.**

---

### Finding 6 — COF threshold metadata was misleading

**Codex claim:** already resolved — thresholds are now defined as named constants and used
consistently in both the overlap calls and the output metadata.

**My verification — confirmed, already resolved.** Constants at
[make_cooccurrence_masks.py:37-44](../../scripts/make_cooccurrence_masks.py) (`MCS_OVERLAP_THRESHOLD
= 0.20`, `AR_OVERLAP_THRESHOLD = 0.10`, `ETC_THREEWAY_THRESHOLD = 0.00`, `MCS_AR_AR_THRESHOLD =
0.00`, `AR_ETC_ETC_THRESHOLD = 0.01`, `MCS_ETC_ETC_THRESHOLD = 0.00`, `MCS_TC_FILTER_THRESHOLD =
0.10`, `CLOUD_TYPE_TC_FILTER_THRESHOLD = 0.01`) match `cof_identification.md`'s threshold summary
table (lines 209-221) exactly, and are used both in the processing calls (lines 1015-1020,
1027-1039, 1081-1091) and rebuilt into the output `thresholds`/`tc_filtering_threshold` attributes
from the same constants (lines 1404-1415), so metadata cannot drift from behavior again. Agree,
nothing further to do.

---

### Finding 7 — Extreme documentation claimed unsaved per-type precipitation variables

**Codex claim:** already resolved — the script only ever wrote `total_extreme_precip`,
`{type}_count`, and `{type}_frac`; the doc previously claimed `{type}_precip`/`{type}_fraction`
and has been corrected.

**My verification — confirmed, already resolved.** In
[calc_stormtype_extreme_precip_spatial.py](../../scripts/calc_stormtype_extreme_precip_spatial.py),
`ds_out['total_extreme_precip']` (line 683) is the only precipitation-sum field written; the
per-type loop (lines 692-708) writes only `{key}_count` and `{key}_frac`. `{type}_precip` never
appears as an output variable name in the current script — the per-type precipitation sums live
only in the internal `accumulated_precip` dict (lines 611, 632-637). `extreme_precip_by_stormtype.md`
lines 190-199 now list exactly `total_extreme_count`, `total_extreme_precip`, `{type}_count`,
`{type}_frac`, matching. Agree, nothing further to do.

---

### Finding 8 — Extreme notebook residual should use `unassigned_frac` and mask no-event cells

**Codex claim:** the script already computes and saves `unassigned_frac` correctly (0.0 at
zero-total cells like every other category); the notebook instead computes
`res_frac = 100 - sum(12 named fractions)` without masking zero-total cells, so every cell with no
extreme precipitation shows as 100% residual, which can dominate the rank map.

**My verification — confirmed, mechanism and math both check out.**
- `storm_type_keys` includes `'unassigned'` (lines 606-608), and the fraction loop
  (lines 661-668) applies the same rule to every key, `unassigned` included:
  `frac = xr.where(total_extreme_precip > 0, accumulated_precip[key] / total_extreme_precip, 0.0)`.
  So at a cell with `total_extreme_precip == 0`, **all 13** `{type}_frac` fields — including
  `unassigned_frac` — are set to exactly 0.0, not NaN.
- The notebook never reads `ds.unassigned_frac`; instead it computes
  `res_frac = 100 - (mcs_iso_frac + ar_iso_frac + etc_iso_frac + tc_frac + mcs_ar_frac +
  mcs_etc_frac + ar_etc_frac + mcs_ar_etc_frac + dc_frac + nd_frac + st_frac + dz_frac)`
  (lines 405-409). At a zero-total cell all twelve named terms are 0, so this formula evaluates to
  exactly 100, and that 100 is not distinguishable from a genuine cell that is "100% unassigned
  extreme precipitation." This residual then enters the ranking dictionary as code 13
  (lines 973, 1128, 1138) and can dominate rank-1 in regions with few or no extreme events.
- I verified the underlying identity holds: since the 13 categories are mutually exclusive and
  exhaustive by construction (`assigned_mask` accumulates through the whole priority chain,
  `process_single_timestep()` lines 451-564), at any cell with `total_extreme_precip > 0` the 13
  fractions sum to exactly 1, so `100 - sum(12 named) == 100 * ds.unassigned_frac` at valid cells.
  The *only* defect is the missing `total_extreme_precip > 0` mask for zero-total cells — the
  arithmetic Codex's fix relies on is otherwise sound.

**Fix assessment: correct.** Masking every fraction (including a residual term) with
`ds.total_extreme_precip > 0` before ranking is the right fix, and reading `unassigned_frac`
directly is a clean way to get the residual without re-deriving it from a sum (it also sidesteps
needing every named term to be present/correctly summed, which is one less place for a future
edit to silently break the identity). Comparing `100*ds.unassigned_frac` against the current
`res_frac` before/after masking, as Codex suggests, is a good sanity check to run once fixed.

**Elevated to P0 (per user feedback and my own cross-reference):** the user connected the CASESM2
"large residual" symptom to Finding 5. Independently, I'd flag that *this* bug is a second,
unrelated candidate explanation for the same symptom — a spurious 100% residual at every
zero-extreme-precipitation cell has nothing to do with real storm-type attribution or the Finding 5
cloud-type bug, and would show up identically regardless of whether Finding 5 fires. Both should be
fixed together before the user tries to interpret the CASESM2 pattern, so the two effects aren't
conflated. (My own hypothesis; not verified against real data in this session.)

**Fix applied.** In `notebooks/plot_cof_extreme_raintype_rank_map.ipynb` (cell id `6fb70234`), added
`valid_extreme = ds.total_extreme_precip > 0`, applied `.where(valid_extreme)` to all 12 named
fraction variables, and replaced the `res_frac = 100 - (...)` formula with
`res_frac = (100 * ds.unassigned_frac).where(valid_extreme)`. Single self-contained cell edit; no
other cell needed to change since everything downstream references these Python variables by name.
Verified structurally: only this one cell's source changed relative to the file's prior state.

---

## Additional items identified in this review (not raised by Codex)

**A — Low / documentation. Extreme-attribution output filename doesn't match the doc. Fixed.**
`extreme_precip_by_stormtype.md` documented the Stage 2 output as
`{source}_stormtype_extreme_precip_{Pxx}_hp8_v1.nc` (at line 74 in the file at review time — I
originally mis-cited this as line 76; corrected during implementation by re-verifying against the
live file), but
[calc_stormtype_extreme_precip_spatial.py:989](../../scripts/calc_stormtype_extreme_precip_spatial.py)
actually writes `{source_name}_stormtype_spatial_{pname.lower()}{date_suffix}.nc` — different
stem (`stormtype_spatial` vs. `stormtype_extreme_precip`), lowercase percentile, and an optional
date-range suffix the doc didn't mention. (The Stage 1 threshold-file name matched already.)
**Fix applied:** the doc's Output line now matches the actual filename pattern, including the
conditional `date_suffix`.

**B — Low / cosmetic. `total_extreme_precip` units label is imprecise. Fixed (resolved per user
feedback: keep units, sharpen description).** `ds_out['total_extreme_precip'].attrs['units']` is
`'mm/h'`, but the variable is a **sum of instantaneous mm/h samples across all extreme time
steps** — not a rate, and not multiplied by a time step to become an accumulated depth in mm. This
doesn't affect any `{type}_frac` ratio (the units cancel), but could mislead a reader taking the
raw value at face value. The user's decision: keep `units: 'mm/h'` unchanged (don't invent a new
unit convention), and instead sharpen the `description` attribute to state this explicitly.
**Fix applied:** the `description` attribute in `calc_stormtype_extreme_precip_spatial.py` now
reads "Sum of instantaneous precipitation rate (mm/h) over extreme time steps (pr > threshold).
This is a sum over samples, not time-integrated by the sampling interval, so it is not a physical
accumulated depth." Mirrored as a parenthetical in `extreme_precip_by_stormtype.md`'s Output
Variables table. I also added a downstream-usage note in that doc's Step 6 pointing consumers at
`unassigned_frac` directly and at masking `total_extreme_count == 0` — directly relevant to
Finding 8's bug, to help prevent the same mistake recurring in a future consumer.

**C — Low / dead code. `dominant_type` is documented and encoded but never created. Fixed (resolved
per user feedback: remove, don't implement).** `save_spatial_results()` claimed in its own
attribute text "Use dominant_type to identify dominant contributor at each location" and set up
NetCDF encoding for a `dominant_type` variable, but `process_timeseries_dask()` never added a
`dominant_type` variable to `ds_out` anywhere in its body — harmless (the encoding branch was
simply never exercised) but confusing to a future reader. The user chose removal over
implementation, to avoid confusion. **Fix applied:** both the attribute-text mention and the dead
`elif var == 'dominant_type':` encoding branch were removed from
`calc_stormtype_extreme_precip_spatial.py`.

**D — New, at user's request. Explicit/defensive AR and ETC TC-overlap filtering. Fixed.** The
user corrected my original Finding 4 claim that TC overlap was unhandled for AR/ETC: it is already
expected to be removed **upstream of this repository**, in the AR/ETC/TC mask generation step
maintained by colleague Bryce Harrop. To make this guarantee explicit and defensive within this
repo too (rather than only assumed from an external process), the user asked for an AR/ETC
TC-overlap filter to be added to `make_cooccurrence_masks.py`, choosing pixel-level exclusion (no
percentage threshold) — the same mechanism as the Finding 5 cloud-type fix — over a track-level,
threshold-based filter like MCS's. **Fix applied:** immediately after the existing MCS/cloud-type
TC filtering in Step 1 of `process_single_timestep_overlaps()`:
```python
_ds = _ds.assign(
    ar_mask=xr.where(tc_present, 0, _ds.ar_mask),
    etc_mask=xr.where(tc_present, 0, _ds.etc_mask),
)
```
`_ds.assign()` returns a new local dataset (no mutation of the caller's object), and reassigning the
local `_ds` name means every later reference to `_ds.ar_mask`/`_ds.etc_mask` in the rest of the
function — Steps 2-7, `extract_etc_overlap_info()`, and the returned masks — automatically picks up
the filtered version with no other call sites needing to change. `_ds.ar_mask`/`_ds.etc_mask` are
genuine track-ID masks (confirmed: the pairing/promotion logic elsewhere in the file extracts
unique AR/ETC track IDs from them, e.g. `_ds.ar_mask.isin(ar_tracks_paired_with_mcs)`), so — unlike
the cloud-type case in Finding 5 — this is filtering real per-object masks, just without a
percentage threshold. `cof_identification.md`'s Step 1 section and threshold table were updated to
document both the MCS (track-level, 10%) and cloud-type/AR/ETC (pixel-level) mechanisms, and to
note the AR/ETC filter is a defensive safety net given upstream filtering. **This changes the
`cofmasks` output and must be rerun before downstream monthly/extreme scripts and notebooks are
rerun** (bundled with the Finding 5 rerun, since both are in the same function).

---

## Consolidated fix priority

Revised per user feedback: Findings 5 and 8 elevated to P0 (tied to the CASESM2 anomaly the user is
actively investigating — see the Synthesis note under Finding 8). **All items are now fixed.**

| Priority | Item | Why | Status |
|---|---|---|---|
| **P0** | Finding 2 (ST/ND swap) + matching doc fix | Reverses scientific interpretation of two categories in a paper-bound figure | Fixed |
| **P0** | Finding 5 (class-level TC filter → pixel-level) | User-elevated; possible contributor to the CASESM2 large-residual anomaly regardless of quantified magnitude | Fixed |
| **P0** | Finding 8 (extreme notebook residual no-event mask) | Second, independent candidate explanation for the same CASESM2 symptom (spurious 100% residual at no-event cells) | Fixed |
| **P1** | Finding 4 (monthly notebook residual arithmetic) | Wrong regional residual stats even in a clean top-to-bottom run; a second, worse double-counting bug was found during implementation | Fixed |
| **P1** | Finding 1 (empty-swath guard) | Cheap, prevents a hard pipeline crash on edge-case inputs | Fixed |
| **P1** | New Item D (explicit AR/ETC TC filtering) | Defensive safety net requested by user; upstream (outside this repo) already expected to handle this | Fixed |
| **P2** | Item A (doc filename mismatch) | Cosmetic/documentation only | Fixed |
| **P2** | Item B (`total_extreme_precip` description) | Cosmetic; units kept as `mm/h` per user decision, description sharpened | Fixed |
| **P2** | Item C (remove `dominant_type` dead code) | Cosmetic; removed per user decision | Fixed |
| Done | Findings 3, 6, 7 | Already corrected in the prior pass per both reports | No action needed |

**Re-run note (applies to Findings 2, 5, and New Item D):** these change the values written to the
`cofmasks` and extreme-attribution NetCDF/zarr outputs. The affected scripts must be rerun
end-to-end, in pipeline dependency order (`make_mcs_swath_masks.py` → `make_cooccurrence_masks.py`
→ `calc_stormtype_extreme_precip_spatial.py` / `calc_monthly_rainmap_by_cof.py`), before the
plotting notebooks are rerun, or the notebooks will keep visualizing stale data. Findings 1, 4, and
8 are notebook/upstream-script fixes that take effect on the next rerun of their respective
script/notebook without further code changes needed.

**Verification performed in this session:** every script edit was checked with `python -m
py_compile` (no syntax errors) and grepped for dangling references to anything removed; every
notebook edit was verified via a structural, cell-by-cell diff confirming only the intended cell's
source changed. **Not performed in this session (no data/cluster access available):** actually
executing any script against real data, or running either notebook end-to-end. That verification —
confirming the CASESM2 anomaly actually shrinks or resolves, and that no output shape/dtype
regressions occur — was the user's planned next step, and has since been completed; see
**"Post-fix validation"** below for the confirmed results.

---

## Post-fix validation (pipeline + notebooks rerun by the user)

After the fixes above were applied, the user reran the full pipeline end-to-end (COF mask
creation → monthly total precipitation statistics → extreme precipitation attribution) and reran
the four affected plotting notebooks (`plot_cof_annual_map.ipynb`,
`plot_cof_total_raintype_rank_map.ipynb`, `plot_cof_extreme_raintype_rank_map.ipynb`,
`plot_cof_total_rain_seasonal_zonalmean_maps.ipynb`) against the newly regenerated output. This is
the empirical check that neither report could perform without data/cluster access. Comparing the
new figures against the pre-fix versions, the user confirms:

1. **Total precipitation ranking by feature type:** the change is minimal or negligible. This is
   expected — none of the fixes altered the monthly total-precipitation script or its ST/AR/ETC
   accounting; the total-precipitation pipeline was not a target of any finding.

2. **Extreme precipitation attribution:** the change is more substantial, and directly confirms two
   of this report's findings:
   - **Confirms Findings 5 and 8 (P0):** the CASESM2 anomaly — large regions where "residual"/
     unassigned was previously the top extreme-precipitation category — is resolved; those regions
     now correctly attribute to defined feature types. Both hypothesized contributors (Finding 5's
     class-level cloud-type TC filter, capable of a full-class wipeout, and Finding 8's spurious
     100% residual at every no-extreme-event cell) were fixed together, and the anomaly's
     disappearance is consistent with one or both actually having driven the CASESM2 symptom, as
     hypothesized in this report's original Finding 8 synthesis note.
   - **Confirms Finding 2 (P0):** the previously mislabeled **ST (stratiform)** category is now
     correctly **ND (non-deep convective)**, and vice versa, as intended by the ST/ND swap fix.
   - Pie charts for most models and regions remain qualitatively similar; the two changes above
     (CASESM2's anomaly resolving, and the systematic ST↔ND relabeling) are the dominant visible
     differences.

3. **Feature frequency maps:** magnitudes increased after the fix, but the spatial pattern is
   unchanged. This is consistent with a separate, earlier fix in this session (not one of the 8
   Codex findings): `plot_cof_annual_map.ipynb`'s frequency calculation was corrected to divide
   mask time-step counts by the number of time steps (`ntimes`) rather than by elapsed hours
   (`nhours`) — see the notebook-accuracy review conducted earlier in this session. Since each time
   step spans a 6-hour window (`nhours = 6 × ntimes`), dividing by `ntimes` instead of `nhours`
   scales every frequency value up by a factor of about 6, while the relative spatial pattern
   (which cells have relatively more or less feature presence) is unaffected — exactly matching the
   user's observation.

4. **Zonal-mean plots** (precipitation and precipitation fraction by feature type): unaffected. This
   is expected, since none of the fixes altered the underlying precipitation accumulation by
   feature type used in `plot_cof_total_rain_seasonal_zonalmean_maps.ipynb` — only the ST/ND
   labeling, TC filtering, and extreme-attribution residual logic changed, none of which this
   notebook's zonal-mean quantities depend on.

**Net assessment:** the fixes behaved as intended — they materially corrected the extreme
precipitation attribution (resolving the CASESM2 anomaly and the ST/ND mislabeling) while leaving
total precipitation and zonal-mean statistics essentially untouched, consistent with the fixes'
scope. No unexpected regressions were reported.

---

## Methodology paragraphs for the paper draft

The paragraphs below are my own rewrite of the "Methodology Paragraphs For Paper Draft" section of
the Codex report, written for a general atmospheric-science audience (e.g., graduate students) and
verified against the current code. They assume the Finding 2 (ST/ND) fix has been applied, since
these paragraphs describe the *intended* science, and are written using the **upstream cloud-type
convention**: type 2 = stratiform, type 3 = non-deep convective.

**MCS swaths and non-tracked cloud types.** Hourly, object-tracked MCS masks, cloud-top brightness
temperature, and precipitation are grouped into sub-daily windows (6 hours in this study, aligned
to 00/06/12/18 UTC). Within each window, every tracked MCS is collapsed into a *swath*: the union
of all grid cells it touched at any hour during the window, so that a moving system is represented
by its full footprint rather than a single instantaneous snapshot. Where two systems' swaths
overlap the same cell, the cell is assigned to whichever system occupied it for more hours during
the window. Outside the MCS swaths, every hourly grid cell is sorted into one of four cloud types
using its brightness temperature relative to a latitude-dependent threshold (250 K in the tropics,
decreasing linearly to 230 K by 60° latitude and staying flat poleward of that, to account for the
naturally colder tropopause at high latitudes) together with whether precipitation exceeds a
0.5 mm/h cutoff: deep convective (cold and raining), stratiform (cold with little or no rain),
non-deep convective (warm and raining), and drizzle (warm with light rain). Each non-MCS cell then
receives one label for the whole window by priority — deep convective outranks stratiform, which
outranks non-deep convective, which outranks drizzle — so a cell that is deep convective even
briefly during the window is labeled deep convective for the whole window. Each cloud type's
window-mean precipitation rate is computed as the average rate during the hours the cell held that
label, scaled by how often that label occurred; the four cloud types' contributions add up exactly
to the simple time-mean precipitation for that window. Both the cloud-type label and its
precipitation are forced to zero wherever the MCS swath is present, so MCS and non-MCS cloud-type
precipitation are always counted separately with no double-counting.

**Co-occurrence feature (COF) identification.** At every time step, the gridded masks of
mesoscale convective systems (MCS), atmospheric rivers (AR), and extratropical cyclones (ETC) are
compared to identify where two or more feature types spatially coincide. Before any comparison, MCS
objects that overlap a tropical cyclone by 10% or more of their area are set aside, so that
TC-embedded convection is not later miscategorized as an independent MCS-AR or MCS-ETC system.
Any AR, ETC, or non-tracked cloud-type grid cell that directly coincides with a TC grid cell is
also removed, regardless of how small that overlap is — unlike the MCS case, this is a per-pixel
exclusion rather than a percentage-of-area threshold, since AR/ETC track objects and cloud-type
classes are large enough that a fractional threshold could let TC-embedded pixels slip through.
Co-occurrence is then accepted using area-overlap thresholds that reflect each feature type's
typical physical size: because an MCS is much smaller than an AR, which is in turn much smaller
than an ETC, any co-occurrence involving an MCS requires that at least 20% of the MCS's own area be
involved, an AR must contribute at least 10% of its own area whenever an ETC is part of the
co-occurrence, and an ETC's threshold is set near zero in every case, since even a small fraction of
a synoptic-scale cyclone's area is still large in absolute terms. (An MCS paired only with an AR, with
no ETC involved, does not additionally require a minimum AR fraction — the MCS-side threshold alone
is enough to accept that pairing.) Three-way
co-occurrences (MCS, AR, and ETC all overlapping) are identified first and require the three
objects to genuinely share common grid cells, not merely each appear somewhere in the same broad
region. Remaining objects are then tested pairwise for MCS-AR, MCS-ETC, and AR-ETC overlap using
the same size-aware thresholds. Finally, if a single ETC is found to independently pair with both an
AR and an MCS — even if all three never share a grid cell simultaneously — that group is promoted
to a three-way co-occurrence, since the ETC effectively bridges the AR and MCS into one physically
connected system; this promotion is applied transitively so that all connected objects end up in the
same category. Every feature ends up assigned to exactly one mutually exclusive category — isolated,
one of the three two-way pairings, or the three-way group — and these masks, together with the TC
mask and the non-tracked cloud-type fields, form the basis for all downstream precipitation
attribution.

**Monthly precipitation by feature type.** The co-occurrence masks are matched in time with a
model or observational precipitation dataset, and for each calendar month the precipitation
falling under each feature category is summed at every grid cell, with each time step's
contribution scaled by the length of that time step so the result is a physical accumulation in
millimeters. Because the co-occurrence step records each feature's own share of a pairing
separately (for example, the MCS's footprint within an MCS-AR pair and the AR's footprint within
that same pair are stored as two different fields), these perspectives are first merged into one
combined footprint per co-occurrence type before summing, so a co-occurring grid cell is counted
once rather than twice. Precipitation from the four non-tracked cloud types is taken from the
pre-computed, frequency-weighted cloud-type fields described above and is only counted at grid
cells that fall outside every tracked feature (MCS, AR, ETC, TC, and all of their co-occurrence
combinations), preserving mutual exclusivity across the whole classification. Alongside the
precipitation totals, two occurrence counts are recorded for each category: how many sub-daily
time steps that category was present, and how many of those time steps had precipitation above a
0.1 mm/h threshold; both are stored as step counts and can be converted to elapsed hours by
multiplying by the time-step length. Percentage contributions to total precipitation are not
computed by this script — they are derived afterward, as a category's accumulated precipitation
divided by the total accumulated precipitation over the same set of grid cells.

**Extreme precipitation by feature type.** A location-specific extreme-precipitation threshold is
first computed at every grid cell as a high percentile (for example, the 90th or 95th) of that
cell's own precipitation time series, after first discarding near-zero precipitation values so that
frequent light rain does not artificially lower the threshold. At every time step, any grid cell
whose precipitation exceeds its local threshold is flagged as an extreme event and assigned to the
storm-type category it belongs to. Because the co-occurrence categories are constructed to be
mutually exclusive at each grid cell in the large majority of cases (isolated vs. two-way vs.
three-way partitioning at the track level, with non-tracked cloud types explicitly excluded from
every tracked-feature category), this assignment is usually unambiguous. A fixed priority order is
applied only as a deterministic tie-break for the small residual of cases where a cell could
otherwise match more than one category, so that no event is ever double-counted: three-way
co-occurrences are considered first, then each of the two-way co-occurrences, then isolated MCS,
AR, and ETC, then tropical cyclones, then the four non-tracked cloud types (deep convective,
non-deep convective, stratiform, and drizzle, in that order), and finally an "unassigned" category
for any extreme event that does not fall into any of the preceding groups. This tie-break order is
a fixed convention, not a claim that one storm type physically takes precedence over another where
their footprints coincide. Counts and precipitation amounts for each category are accumulated
across the full analysis period,
and each category's contribution at a grid cell is expressed as that category's share of the cell's
total extreme precipitation. Because every extreme event is assigned to exactly one category, these
shares sum to one (100%) at any cell that experiences at least one extreme event. When building
summary maps, the stored "unassigned" share should be used directly for the residual category, and
grid cells that never experience an extreme event over the analysis period should be masked out
rather than treated as fully unassigned, since a cell with no extreme precipitation has an
undefined — not 100% residual — attribution.

---

## Notes on how this review was conducted

- Every code citation in this report was re-derived directly from the current files in the
  workspace (not copied from the Codex report), and cross-checked line-by-line during writing.
  Two small citation corrections relative to the Codex report are noted inline (Finding 2's doc
  lines are 161-162, not 160-161; Finding A's doc line is 74, not 76 — I originally mis-cited this
  as 76 in the first pass of this report and corrected it by re-verifying against the live file
  immediately before editing it) — these are one-line discrepancies from careful re-verification,
  not disagreements about substance.
- I did not execute any script or notebook as part of the original review pass (no cluster/data
  access was used); all conclusions there are from static code reading, including the
  synthetic-input reasoning (e.g., all-zero `mcs_mask`, `cloud_types_cleaned` values of 2 vs. 3)
  that both this report and the Codex report rely on. The recommended synthetic unit tests and the
  CASESM2 P95 diagnostic listed in the Codex report's "Recommended Independent Review Checks"
  section remain good follow-up actions and are not superseded by this review.
- **During implementation of the fixes** (this revision), I re-verified every target file's exact
  current content immediately before editing it, rather than trusting the line numbers/content
  recorded in the original review pass. This caught two things worth flagging: (1) the doc-line
  citation correction noted above, and (2) `notebooks/plot_cof_raintype_rank_map.ipynb` had been
  edited locally since the original review (most likely by the user, working in the same session) —
  its residual formula had drifted to a *different and more broken* state than what Finding 4
  originally described (a double-counting bug on top of the double-scaling bug; see Finding 4's
  "second bug found during implementation" note). The fix applied targets this current, verified
  content, not the original (now-stale) analysis.
- Notebook edits were made by direct, minimal JSON manipulation (loading the `.ipynb` as JSON,
  replacing only the target cell's `source` field, writing back with matching formatting) rather
  than through a full notebook read/write cycle, because both notebooks exceed this session's file
  read size limit. Each edit was verified with a structural, cell-by-cell diff against the
  git-committed version confirming only the intended cell's source changed. Note for the user: both
  notebooks' on-disk execution outputs were already substantially out of sync with git history
  *before* this session's edits began (many cells' `execution_count`/`outputs` differed from the
  last commit) — this predates and is unrelated to my edits, and will be superseded regardless once
  the notebooks are rerun.
- No script or notebook was executed against real data in this session (no cluster/data access
  available). Verification performed instead: `python -m py_compile` on every edited script,
  `json.load` validation on every edited notebook, and `grep` checks for dangling references to
  anything removed (e.g., `CLOUD_TYPE_TC_FILTER_THRESHOLD`, `dominant_type`). Confirming the actual
  numerical/scientific impact of these fixes — including whether the CASESM2 anomaly shrinks or
  resolves — was the user's planned next step (pipeline rerun, then notebook rerun); the user has
  since completed this rerun and reported confirmed results, recorded in **"Post-fix validation"**
  above. That validation was performed by the user against real data outside this session; it is
  not something I executed or independently re-verified numerically.
