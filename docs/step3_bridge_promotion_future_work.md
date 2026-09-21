# Step 3 Bridge Promotion: Not Adopted (Plan and Rationale Kept for Future Work)

**Reference script:** `scripts/make_cooccurrence_masks.py` (`process_single_timestep_overlaps`, `promote_dual_etc_overlaps_to_3way`)

**Author:** Zhe Feng | zhe.feng@pnnl.gov

---

## Status

**Decision of 2026-09-20: the bridge promotion described below is not adopted. Step 3 is unchanged.** The existing rule (an ETC that touches both an MCS and an AR promotes the group to 3-way) stays as it is, and the priority order stays in the consumers of the COF masks (monthly maps, extreme-precipitation attribution). Nothing in the code was changed for this proposal.

This page records what was planned, how large the effect is, why it was not adopted, and what a later revisit would have to do. Of the three options in `AUDIT_FINDINGS.md` (Open follow-ups, item 1): option A (bridge promotion) is the subject of this page and is not adopted; option C (leave Step 3, keep the priority order in the consumers) stays; option B (make "isolated" exclusive of 3-way footprints in Step 3) was not part of this decision and remains open.

## How Step 3 assigns categories today

Step 3 works on one 6-hourly frame at a time:

1. `find_true_3way_overlaps` finds tracks that share a pixel where MCS, AR and ETC all overlap (the true 3-way tracks).
2. `find_overlapping_tracks_and_pairs` builds three pair lists (MCS-AR, AR-ETC, MCS-ETC). A pair is two tracks whose overlap exceeds the thresholds (share of the track's pixels inside the overlap: MCS 0.20 in MCS-AR and MCS-ETC, AR 0.10 in AR-ETC, ETC 0.01 in AR-ETC; the other partners 0.0).
3. Pairs that contain a true 3-way track are dropped from the pair lists.
4. `promote_dual_etc_overlaps_to_3way` promotes an ETC that sits in both the AR-ETC and the MCS-ETC list to a 3-way group (the "core" MCS-AR-ETC triplets), adds every track that has a pair with a core track (one hop), and removes every pair that contains a promoted track.
5. The masks are built per track: a track that is in a list gets its **whole footprint** into that category (`isin` on track IDs); an isolated object is a track that is in no list and in no 3-way set.

Nothing stops one track from being in two lists, and an isolated object can lie inside the footprint of a group. The monthly and extreme scripts resolve such overlaps per pixel with a priority order (MCS-AR-ETC > MCS-AR > MCS-ETC > AR-ETC > isolated MCS > isolated AR > isolated ETC > TC, then the cloud types), so the published totals count every pixel once. Only per-category quantities that are not resolved that way, such as the track flags of Analyses 3 and 4, see the categories as Step 3 defines them.

### What the existing promotion does on toy input

Run with the existing function on hand-made pair lists (IDs: MCS 1, 2; AR 10, 11; ETC 100, 101), 2026-09-20:

| Case | Pairs | Result today |
|------|-------|--------------|
| 1. AR bridge | MCS 1-AR 10, AR 10-ETC 100 | Nothing changes. AR 10 stays in both the MCS-AR and the AR-ETC list. |
| 2. ETC bridge | MCS 1-ETC 100, AR 10-ETC 100 | All three are promoted to 3-way and all pairs are removed. The only case handled. |
| 3. MCS bridge | MCS 1-AR 10, MCS 1-ETC 100 | Nothing changes. MCS 1 stays in both lists. |
| 4. Case 2 plus a second AR | as case 2, plus MCS 1-AR 11 and AR 11-ETC 101 | MCS 1, AR 10, AR 11 and ETC 100 are promoted. ETC 101 is not, and its only pair (AR 11-ETC 101) is removed, so ETC 101 becomes an isolated ETC although it overlaps AR 11. |

The ETC-bridge rule (case 2) never fires on the six datasets (0 frames in the full records), so it, and the case-4 behaviour, have never run on real data.

## How large the effect is

E is the precipitation that would be counted twice if the eight feature categories were simply added (0% in the published totals, because of the priority order). Measured on the round-2 COF stores over the full records with the audit tooling (`overlap_full.py`, not in the repository):

| Source | E at 60S-60N | E at 30-60N / 30-60S | AR in two pair lists (% of AR track-frames) | MCS in two pair lists (% of MCS track-frames) | Share of E that a bridge promotion removes (lower bound) |
|--------|--------------|----------------------|------------------------------------------|-----------------------------------------------|-------------------------------------------------------|
| SCREAM | 2.18% | 5.8% / 5.0% | 4.4% | 0.13% | 23% |
| ICON | 2.71% | 6.2% / 6.3% | 4.5% | 0.25% | 31% |
| NICAM | 2.37% | 5.3% / 5.3% | 6.5% | 0.15% | 50% |
| UM | 1.63% | 3.9% / 3.4% | 7.0% | 0.07% | 64% |
| CASESM2 | 0.65% | 1.6% / 1.7% | 3.9% | 0.14% | 54% |
| IMERG | 2.56% | 5.6% / 5.7% | 5.9% | 0.15% | 42% |

Two different things sit inside E. One track in two pair lists (an AR that overlaps an MCS in one place and an ETC in another; an MCS that overlaps an AR and an ETC) is 7-23% of E and is what a promotion addresses. Different objects overlapping, mostly an isolated MCS inside the footprint of a 3-way or 2-way group (no valid pair between them), is 77-93% of E and cannot be touched by a promotion.

## The change that was planned (not implemented)

**Rule (option A).** Generalise the existing ETC rule to a track of any type that sits in two pair lists, a "bridge": an ETC in AR-ETC and MCS-ETC (existing), an AR in MCS-AR and AR-ETC (new), an MCS in MCS-AR and MCS-ETC (new). The bridge and its partners in the two lists form the core triplets; the existing one-hop expansion and the removal of every pair with a promoted member follow.

```python
ar_bridges  = {a for _, a in mcs_ar} & {a for a, _ in ar_etc}    # AR in MCS-AR and AR-ETC
mcs_bridges = {m for m, _ in mcs_ar} & {m for m, _ in mcs_etc}   # MCS in MCS-AR and MCS-ETC
etc_bridges = {e for _, e in ar_etc} & {e for _, e in mcs_etc}   # ETC in AR-ETC and MCS-ETC (the only one handled today)
```

**Code.** One pure function on integer track sets (no xarray; IDs normalised to `int`, sorted output), the old function kept as `legacy`; an option `--promotion legacy|bridge` (default `legacy` until a difference report is approved), written into the store attributes. No threshold, pixel operation or mask construction changes. The mode has to reach the Dask workers explicitly: they call `process_single_timestep_overlaps(_ds, verbose)` through `process_timestep_wrapper_zarr` in `src/zarr_tools.py` (and the sequential fallback in `stream_process_to_zarr`), and no option reaches them today. Each frame would return the mode it used, and the main process would check that all frames agree.

**Why it can be checked.** After the promotion no track can remain in two lists: a track in two lists is a bridge, so it is in the core set, is promoted, and every pair containing it is removed (no iteration needed). If a frame has no AR or MCS bridge the output is identical to today's, and the union of all tracked pixels never changes; only the category of a bridged track and its partners does.

**Open design choice: what happens to a track that touches a promoted group (case 4).** (a) One hop, like the existing rule: the touching track is cut loose and becomes isolated (an orphan). (b) Closure: connected components of the pair graph, every track in exactly one system and no orphans, but it also moves tracks that today are dropped from 2-way lists because their partner is a true 3-way track, a larger change. The choice was to be made with data (stage 0 below), not assumed.

**Stages.**
0. Read-only measurement, no repository change: dump per-frame tracks, pairs and per-track precipitation for the six sources once, then compare `legacy`, one-hop and closure offline (bridges, orphans, tracks whose valid pair was dropped, predicted movement between categories, remaining E).
1. The pure function with tests (`tests/test_promotion.py`: the toy cases above, a randomised cross-check against an independent union-find implementation, order independence, integer/float IDs), then the option.
2. Docs.
3. A full-record run of the six sources in a test root with the runner (`--step-args s3 "--promotion bridge"`, about 1 h), with a difference report: masks and invariants (union of tracked pixels, `tot_pr`, cloud types, TC identical; E not increased; no track in two lists), monthly categories (closure stays 0.0000%), extreme shares (unassigned 0.000%), the ETC overlap parquet and the MCS track flags, determinism across worker counts.
4. Decision: adopt (default to `bridge`, re-run) or keep `legacy`.

**Downstream effects to expect.** Precipitation would move between the mcs_ar, ar_etc and mcs_ar_etc categories (total precipitation and the residual unchanged). The ETC `overlap_flag` of an ETC whose AR is promoted would change from "AR only" (2) to "MCS and AR" (3) (Analysis 3), and the 3-way lifetime fractions of MCS tracks would rise slightly (Analysis 4).

**What it would not fix.** The isolated-inside-a-footprint part of E (roughly a fifth to a third of E, depending on the source). That part is already handled by the priority order in the consumers.

## Rationale for not adopting it

- **An ETC is the main dynamical driver.** An ETC is likely to influence both the AR and the MCS it touches. Treating an MCS and an AR that separately touch the same ETC as one 3-way system is therefore physically justified. This is the existing rule and it stays.
- **An AR (or MCS) bridge is less physical.** Using an AR as the bridge to group an MCS and an ETC that do not touch each other, or an MCS as the bridge between an AR and an ETC, does not have that dynamical justification.
- **Case 4 stays an edge case.** Adding the extra ETC (ETC 101) to the group sounds logical geometrically, but the resulting group is harder to explain dynamically. The current algorithm deliberately does not include it, and that stays as a known edge case rather than being closed by a closure rule.

## What stays true today

- Track-level double membership exists: an AR is in both the MCS-AR and the AR-ETC list in 3.9-7.0% of AR track-frames (0.18-0.46 per 6-hourly frame), an MCS in both MCS-AR and MCS-ETC in 0.07-0.25% of MCS track-frames. Published totals are not affected (priority order); per-category track counts and the Analysis 3 and 4 flags see the categories as defined in the section above.
- The COF categories are exclusive after the priority order of the consumers, not in the Step 3 masks themselves. `cof_identification.md` and `mcs_cof_trackstats.md` say so.
- Step 3 exits with status 0 even when a frame fails (the frame stays NaN and only a log line says so). This is independent of the decision above and is handled separately.

## If this is revisited

1. Decide with a physical argument first which groupings are meaningful (ETC-centred only, or also an AR or MCS bridge; one hop or closure), before writing code. The decision above says the ETC-centred rule is the one that is justified.
2. Run stage 0 (measure only) and look at real example frames with a bridge before choosing a rule.
3. Implement behind the `--promotion` option with `legacy` as the default, test as described above, produce the difference report on the full records in a test root, and only then decide on the default.
4. Re-run the monthly and extreme steps and compare; check the Analysis 3 overlap flags and the Analysis 4 flags against the previous ones.
