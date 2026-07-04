# COF Code-Documentation Review Findings

This report records a code-driven review of the COF workflow against the procedure documentation and downstream plotting notebooks. Line references are from the current workspace after the low-risk documentation and metadata corrections applied in this pass.

Behavior-changing fixes were not applied here. In particular, the ST/ND extreme-attribution swap, no-MCS swath edge case, notebook residual logic, and cloud-type TC filtering behavior are documented below for independent review and later reruns.

## Low-Risk Corrections Applied In This Pass

- `scripts/make_cooccurrence_masks.py`: threshold values are now defined as constants at lines 37-44, reused in the processing calls at lines 1014-1038 and 1081-1090, and written into output metadata at lines 1404-1415. This replaces the previous hard-coded and misleading `ETC: 5% (2-way)` metadata.
- `scripts/calc_monthly_rainmap_by_cof.py`: monthly count variable attributes now describe count variables as time-step counts, not hours, at lines 811-934.
- `docs/procedures/monthly_precip_by_cof.md`: Step 6 and the output table now describe `*_count` and `*_precipitation_count` as sub-daily time-step counts, with elapsed hours obtained by multiplying by `time_interval`; see lines 190-201 and 240-241.
- `docs/procedures/extreme_precip_by_stormtype.md`: the output table now matches the script output: `total_extreme_count`, `total_extreme_precip`, `{type}_count`, and `{type}_frac`; see lines 190-199. The zero-total fraction behavior is also corrected to 0.0 at line 182.

## Finding 1 - High - No-MCS 6-Hour Windows Can Crash MCS Swath Processing

**Files reviewed:** `scripts/make_mcs_swath_masks.py` vs. `docs/procedures/mcs_swath_cloud_type.md`.

**Confirmed code behavior:** `process_timechunk_swath()` creates per-track swath dictionaries from the 6-hour `mcs_mask` at lines 486-488. If a 6-hour window contains no MCS pixels, `create_track_swaths_and_coverage()` returns empty dictionaries. `combine_swaths_with_priority()` then tries to get the first dictionary key at lines 113-115, which raises `IndexError`.

**Documentation mismatch:** The swath documentation describes each 6-hour output as a valid aggregated swath and cloud-type product even when no MCS occurs. It does not mention a no-MCS failure mode.

**Scientific impact:** Any model/source/time period with an MCS-free aggregation window can abort the pipeline, preventing complete swath products. This is especially plausible for sparse MCS regions or short subsets.

**Synthetic check:** For `mcs_mask = zeros((6, ncell))`, there are no positive track IDs. The swath dictionary is `{}`, so `list(track_swaths_dict.keys())[0]` fails.

**Suggested fix:** In `process_timechunk_swath()`, handle empty dictionaries before calling `combine_swaths_with_priority()`:

```python
mcs_swaths_dict, mcs_coverage_dict = create_track_swaths_and_coverage(mcs_mask)
if not mcs_swaths_dict:
    combined_mcs_swath = np.zeros(mcs_mask.shape[1:], dtype=int)
else:
    combined_mcs_swath = combine_swaths_with_priority(mcs_swaths_dict, mcs_coverage_dict)
```

## Finding 2 - High - Extreme Attribution Swaps Stratiform And Non-Deep Convective Cloud Types

**Files reviewed:** `scripts/calc_stormtype_extreme_precip_spatial.py`, `scripts/make_mcs_swath_masks.py`, `docs/procedures/mcs_swath_cloud_type.md`, `docs/procedures/extreme_precip_by_stormtype.md`, and `notebooks/plot_cof_extreme_raintype_rank_map.ipynb`.

**Confirmed code behavior:** The upstream swath script defines cloud type 2 as stratiform and cloud type 3 as non-deep convective at `make_mcs_swath_masks.py` lines 199-201 and assigns precipitation components consistently at lines 305-308. The extreme-attribution script reverses those labels: it assigns `nd` from `cloud_types_cleaned == 2` at lines 531-535 and `st` from `cloud_types_cleaned == 3` at lines 538-542.

**Documentation mismatch:** `docs/procedures/mcs_swath_cloud_type.md` is consistent with the upstream code, but `docs/procedures/extreme_precip_by_stormtype.md` repeats the swapped mapping at lines 160-161. This means the extreme procedure documentation currently documents the bug rather than the upstream cloud-type convention.

**Downstream impact:** `plot_cof_extreme_raintype_rank_map.ipynb` reads the already-swapped output variables directly: `nd_frac = 100 * ds.nd_frac` and `st_frac = 100 * ds.st_frac` at lines 400-402. The ranking dictionary uses `nd_frac` and `st_frac` at lines 1136-1139, the map legend labels code 10 as ND and code 11 as ST at lines 1400-1408, and regional summaries average these same variables at lines 775-777 and include them in pie-chart arrays at lines 2585-2608. Therefore the spatial maps and regional mean fractions label stratiform extreme precipitation as ND and non-deep convective extreme precipitation as ST.

**Scientific impact:** The combined contribution of ST+ND is unchanged, but the interpretation of each category is reversed in extreme precipitation maps and regional summaries. This can alter conclusions about the relative roles of stratiform versus non-deep convective precipitation in extremes.

**Synthetic check:** If `cloud_types_cleaned = [2, 3]` and both cells are extreme, the current extreme script counts the type-2 cell under `nd` and the type-3 cell under `st`. Upstream definitions require the opposite.

**Suggested exact code fix:** In `calc_stormtype_extreme_precip_spatial.py`, change only the comparisons/comments in the existing blocks so notebook variable names remain valid after rerunning the analysis:

```python
# ND = 3
mask_nd = (cloud_types_cleaned == 3)
extreme_nd_mask = (extreme_mask & mask_nd & ~assigned_mask)
results['nd'] = extreme_nd_mask.astype(float)
results['nd_pr'] = xr.where(extreme_nd_mask, pr_t, 0.0)
assigned_mask = assigned_mask | extreme_nd_mask

# ST = 2
mask_st = (cloud_types_cleaned == 2)
extreme_st_mask = (extreme_mask & mask_st & ~assigned_mask)
results['st'] = extreme_st_mask.astype(float)
results['st_pr'] = xr.where(extreme_st_mask, pr_t, 0.0)
assigned_mask = assigned_mask | extreme_st_mask
```

After rerunning the extreme attribution files, the existing notebook references to `ds.nd_frac` and `ds.st_frac` should produce corrected maps and regional summaries. Also update `docs/procedures/extreme_precip_by_stormtype.md` lines 160-161 to state `ST = cloud_types == 2` and `ND = cloud_types == 3`.

## Finding 3 - Medium - Monthly Count Variables Are Time-Step Counts, Not Hours

**Files reviewed:** `scripts/calc_monthly_rainmap_by_cof.py`, `docs/procedures/monthly_precip_by_cof.md`, and `notebooks/plot_cof_raintype_rank_map.ipynb`.

**Confirmed code behavior:** The monthly script multiplies precipitation rates by `time_interval` before summing precipitation amounts, for example at lines 303-307, 311-315, 320-332, 338-353, and 372-387. In contrast, occurrence counts and precipitating counts are plain boolean sums over `time`, for example lines 308-309, 312-317, 321-334, 339-355, and 373-388. Thus `*_count` and `*_precipitation_count` are counts of sub-daily samples.

**Documentation/metadata mismatch:** Before this pass, the script attributes and documentation described these fields as hours. The script attributes have now been corrected at lines 811-934, and the procedure documentation has been corrected at `monthly_precip_by_cof.md` lines 190-201 and 240-241.

**Downstream impact:** The current precipitation contribution workflow in `plot_cof_raintype_rank_map.ipynb` is not affected by this count-unit issue. The notebook computes `nhours = ntimes * time_interval` at lines 261-267, converts accumulated precipitation to mm/day using `24 * sum(precipitation) / nhours` at lines 269-291, and computes contribution fractions from precipitation amounts at lines 297-318. It does not use `*_count` variables for the current ranked precipitation contribution plots.

**Synthetic check:** For four 6-hour samples with a feature present at every sample, the code writes `*_count = 4`, not 24. To obtain hours, downstream code must compute `4 * time_interval = 24`.

**Suggested fix status:** Metadata and docs have been corrected. Do not change the numeric count fields unless a future data-interface decision is made to store elapsed hours instead of sample counts.

## Finding 4 - Medium - Monthly Notebook Has A Residual Formula Hazard

**Files reviewed:** `notebooks/plot_cof_raintype_rank_map.ipynb`.

**Confirmed notebook behavior:** The first residual calculation uses already normalized precipitation rates (`tot_pcp_all`, `mcs_pcp_all`, etc.) and then applies another `24 / nhours` factor at lines 292-295. That is dimensionally incorrect. Later, the notebook recomputes `res_pcp_all` more sensibly as `tot_pcp_all - total_iso_pcp_all - total_cof_pcp_all - total_cloudtypes_pcp_all` at lines 997-1011, and the ranking dictionary uses the later `res_pcpfrac_all` at lines 1281-1299 if the notebook is executed top-to-bottom.

**Documentation mismatch:** The procedure documentation describes percentage contribution as a category precipitation total divided by total precipitation. It does not document the notebook's first, dimensionally inconsistent residual calculation.

**Scientific impact:** In a full top-to-bottom run, the later cell should overwrite the earlier residual and reduce the risk for ranked maps. However, executing selected cells out of order, or using the first residual plot at lines 328-329, can show a residual field that is too small by an extra factor of `24 / nhours`. The later residual still assumes that all categories in the subtraction are mutually exclusive; TC overlap with other non-MCS features should be checked if residual closure is used scientifically.

**Suggested fix:** Remove or correct the first residual calculation. Prefer one residual definition after all mutually exclusive category totals are computed:

```python
total_feature_pcp_all = (
    mcs_iso_pcp_all + ar_iso_pcp_all + etc_iso_pcp_all + tc_pcp_all
    + mcs_ar_pcp_all + mcs_etc_pcp_all + ar_etc_pcp_all + mcs_ar_etc_pcp_all
    + dc_pcp_all + nd_pcp_all + st_pcp_all + dz_pcp_all
)
res_pcp_all = tot_pcp_all - total_feature_pcp_all
res_pcpfrac_all = 100.0 * res_pcp_all / tot_pcp_all
```

If TC is not mutually exclusive with AR/ETC in the mask product, either remove overlapping TC pixels from the other categories before closure or describe residual as an approximate diagnostic rather than a strict complement.

## Finding 5 - Medium - Cloud-Type TC Filtering Is Class-Level, Not Object-Level

**Files reviewed:** `scripts/make_cooccurrence_masks.py` and `docs/procedures/cof_identification.md`.

**Confirmed code behavior:** The same helper `filter_mcs_tc_overlaps()` is used for both MCS tracks and cloud types. For MCS, the input values are track IDs, which matches the helper's assumptions at lines 354-381 and 417-487. For cloud types, the call passes `_ds.cloud_types` as `mcs_mask` at lines 1033-1038. Since `cloud_types` contains class labels 1-4 rather than object IDs, the helper computes overlap fractions by cloud class and removes every pixel whose class value exceeds the TC overlap threshold. The returned `cloud_types` field is then used downstream at lines 1041-1042.

**Documentation mismatch:** `docs/procedures/cof_identification.md` says the cloud-type mask is filtered using a 1% TC threshold at lines 86-89 and lists the threshold at lines 211-212, but it does not explain that the current implementation operates on whole cloud classes within each time step rather than on spatially connected cloud objects or only TC-overlapping pixels.

**Scientific impact:** This can inflate unassigned or residual extreme precipitation. Cells inside `tc_mask` are assigned to TC before cloud types in the extreme script, so they should not become residual. But cells near TC, outside the TC mask, or in cloud classes removed globally for that time step can lose their cloud-type label. If those cells are not covered by MCS/AR/ETC/TC masks and have intense precipitation, the extreme script assigns them to `unassigned` at lines 561-564. This is a plausible contributor to large residual areas such as the CASESM2 P95 issue, although it must be separated from the notebook residual issue in Finding 8.

**Suggested targeted review:** For CASESM2 and at least one control source, compare cloud-type area by class before and after TC filtering for each time step. Then plot `100 * ds.unassigned_frac`, `100 * ds.tc_frac`, `ds.total_extreme_count`, and the notebook-computed `res_frac`. If large residual regions appear where `total_extreme_count == 0`, the notebook is the likely cause. If residual is large where `unassigned_frac` is large and cloud-type area was removed by TC filtering, the class-level TC filter is a likely source-side contributor.

**Suggested code direction for later:** Replace the cloud-type call to the track-removal helper with a spatial TC exclusion, such as setting `cloud_types = xr.where(_ds.tc_mask > 0, 0, _ds.cloud_types)`, or implement connected-object labeling for cloud-type regions before applying overlap fractions. This should not be changed without testing because it changes the science product.

## Finding 6 - Low - COF Threshold Metadata Was Misleading

**Files reviewed:** `scripts/make_cooccurrence_masks.py` and `docs/procedures/cof_identification.md`.

**Confirmed code behavior:** Actual thresholds are: MCS 20%, AR 10%, ETC 0% for three-way, AR threshold 0% in MCS-AR, ETC threshold 1% in AR-ETC, and ETC threshold 0% in MCS-ETC. These are now defined at lines 37-44 and used in the overlap calls at lines 1014-1020 and 1081-1090.

**Documentation mismatch:** The procedure documentation already lists the nuanced pair thresholds correctly at lines 213-221. The script metadata previously advertised `ETC: 5% (2-way)`, which did not match the actual pair calls.

**Fix applied:** The output metadata now builds the `thresholds` and `tc_filtering_threshold` attributes from the constants at lines 1404-1415. No behavioral change was made.

## Finding 7 - Medium - Extreme Documentation Claimed Per-Type Precipitation Variables That Are Not Saved

**Files reviewed:** `scripts/calc_stormtype_extreme_precip_spatial.py`, `docs/procedures/extreme_precip_by_stormtype.md`, and `notebooks/plot_cof_extreme_raintype_rank_map.ipynb`.

**Confirmed code behavior:** The extreme script accumulates per-type precipitation internally at lines 631-637, computes per-type fractions at lines 661-668, and writes `{type}_count` plus `{type}_frac` at lines 691-708. It writes only one precipitation-sum field, `total_extreme_precip`, at lines 683-689. It does not save `{type}_precip` variables.

**Documentation mismatch:** Before this pass, `extreme_precip_by_stormtype.md` claimed `{type}_precip` and `{type}_fraction` outputs. The table is now corrected at lines 190-199 to show `{type}_count` and `{type}_frac`, and to note that per-type precipitation sums are internal only.

**Downstream impact:** The extreme plotting notebook uses the saved fraction fields directly at lines 389-403 and does not require `{type}_precip` variables. Therefore this documentation mismatch does not break the current visualization notebook.

**Suggested follow-up:** If future analysis needs regional precipitation-weighted sums rather than means of fractions, consider saving `{type}_precip` variables explicitly. Otherwise keep the documentation aligned with the current compact output.

## Finding 8 - Medium - Extreme Notebook Residual Should Use `unassigned_frac` And Mask No-Event Cells

**Files reviewed:** `scripts/calc_stormtype_extreme_precip_spatial.py` and `notebooks/plot_cof_extreme_raintype_rank_map.ipynb`.

**Confirmed code behavior:** The extreme script includes `unassigned` in `storm_type_keys` at lines 606-608 and writes `unassigned_frac` through the same `{type}_frac` loop at lines 691-708. At cells with no total extreme precipitation, all fractions are set to 0.0 by lines 661-667.

**Notebook behavior:** The notebook does not read `ds.unassigned_frac`. Instead, it computes `res_frac = 100 - sum(named fractions)` at lines 405-409. Since zero-total cells have all saved category fractions equal to 0.0, this formula creates an artificial 100% residual wherever `total_extreme_precip == 0` unless those cells are masked. The residual is then included in the ranking dictionary at lines 1126-1144 and can become the dominant map category.

**Scientific impact:** A large residual region in an extreme rank map can represent real unassigned extreme precipitation, no-event cells being displayed as 100% residual, or both. This is especially relevant to the CASESM2 P95 concern.

**Suggested notebook fix:** Read residual from the script output and mask no-event cells before ranking:

```python
valid_extreme = ds.total_extreme_precip > 0

mcs_iso_frac = (100 * ds.mcs_isolated_frac).where(valid_extreme)
ar_iso_frac = (100 * ds.ar_isolated_frac).where(valid_extreme)
etc_iso_frac = (100 * ds.etc_isolated_frac).where(valid_extreme)
tc_frac = (100 * ds.tc_frac).where(valid_extreme)
mcs_ar_frac = (100 * ds.mcs_ar_2way_frac).where(valid_extreme)
mcs_etc_frac = (100 * ds.mcs_etc_2way_frac).where(valid_extreme)
ar_etc_frac = (100 * ds.ar_etc_2way_frac).where(valid_extreme)
mcs_ar_etc_frac = (100 * ds.mcs_ar_etc_3way_frac).where(valid_extreme)
dc_frac = (100 * ds.dc_frac).where(valid_extreme)
nd_frac = (100 * ds.nd_frac).where(valid_extreme)
st_frac = (100 * ds.st_frac).where(valid_extreme)
dz_frac = (100 * ds.dz_frac).where(valid_extreme)
res_frac = (100 * ds.unassigned_frac).where(valid_extreme)
```

As a diagnostic, compare the current `res_frac` against `100 * ds.unassigned_frac`. Large differences indicate arithmetic/no-event artifacts rather than source-side unassigned precipitation.

## Methodology Paragraphs For Paper Draft

**MCS swath and non-tracked cloud-type aggregation.** Hourly MCS masks are grouped into sub-daily windows, typically 6 hours, and each MCS track is converted to a swath representing all grid cells occupied during the window. Where multiple MCS tracks overlap, the track with the largest number of hourly occurrences at that grid cell is retained. For non-MCS regions, each hourly grid cell is classified into one of four cloud types from brightness temperature and precipitation rate. The 6-hour cloud-type label is then assigned by priority, while cloud-type precipitation is computed separately from all hourly occurrences of each type as a frequency-weighted mean precipitation rate. Both the final cloud-type labels and the cloud-type precipitation components are set to zero wherever the 6-hour MCS swath is present.

**COF mask identification.** At each time step, MCS, AR, and ETC masks are converted to binary feature-presence masks after removing MCS regions that overlap tropical cyclones above the configured threshold. Candidate MCS-AR-ETC, MCS-AR, MCS-ETC, and AR-ETC overlaps are identified from spatial intersections and accepted when feature-specific overlap fractions meet the scale-dependent thresholds. Three-way overlaps are identified first, and ETC systems that separately link an AR and an MCS are promoted to the three-way category. The output masks assign MCS, AR, and ETC features to isolated, two-way, or three-way categories, with additional TC and cloud-type fields carried forward for precipitation attribution.

**Monthly precipitation attribution.** Monthly precipitation totals are computed by aligning the COF mask data with the precipitation data, grouping by month, and summing precipitation rate multiplied by the detected time interval under each category mask. For two-way and three-way COFs, the script first unions the perspective-specific masks so each co-occurrence category is represented by a single spatial footprint. Cloud-type precipitation uses the precomputed frequency-weighted cloud-type precipitation fields and excludes all tracked feature footprints. Percentage contributions are not written by the monthly script; downstream analysis computes them as each category's accumulated precipitation divided by the total accumulated precipitation over the same grid cell or region.

**Extreme precipitation attribution.** Extreme precipitation thresholds are computed independently at each grid cell from the selected percentile of the precipitation time series after excluding trace precipitation. For each time step, cells exceeding the local threshold are assigned to exactly one category using a priority order that places higher-order COFs first, followed by isolated MCS, AR, ETC, TC, non-tracked cloud types, and finally unassigned cells. The script accumulates event counts and precipitation values for each category, then writes per-cell category fractions as category extreme precipitation divided by total extreme precipitation at that cell. Downstream maps and regional summaries should use the saved `unassigned_frac` field for residual contributions and mask cells with no extreme precipitation.

## Recommended Independent Review Checks

- Add a small unit test or standalone synthetic check for a 6-hour MCS-free swath window to confirm the proposed zero-swath fix.
- Add a synthetic extreme-attribution test with `cloud_types_cleaned == 2` and `== 3` to verify that ST and ND are assigned according to the upstream cloud-type convention.
- For monthly outputs, confirm that `*_count * time_interval` gives elapsed hours and that precipitation fractions do not use count variables.
- For CASESM2 P95, plot `total_extreme_count`, `unassigned_frac`, notebook-computed residual, and cloud-type area removed by TC filtering. This separates no-event residual artifacts from source-side unassigned precipitation.
- After any behavior-changing fix, rerun the relevant script before rerunning the notebooks; otherwise the notebooks will continue to visualize stale NetCDF variables.
