#!/usr/bin/env python
"""
Check the tracking inputs of a source before a (re-)processing round, and say how its ETC environment stores can be made.

Input: the ETC stitched-node track file of the new tracking (--new-track-file), the file of the old tracking that the existing environment stores were
extracted for (--old-track-file, default: the registry's production file), optionally the file of the previous round (--prev-track-file), and the
TempestExtremes mask netCDF of the source (ETC_test_tracks_<source_te>_<source_res>.<YYYYMM>.nc, from config_sources.yaml, or --mask-dir/--mask-basename).

Prints and checks:
  1. the new track file against the old one: points and storms, points of the new file that are not in the old one ((time, lon, lat) as the extraction
     sees them), storm order, how the storms changed (unchanged, trimmed at an end, split, removed entirely). Storm IDs are running counters, so a removed
     or split storm renumbers all later ones: only (time, lon, lat) identifies a point.
  2. the ROUTE for the environment stores (they depend only on the ETC points, not on TC or AR, nor on precipitation):
        link     the new file is identical to the previous round's (--prev-track-file, same md5): link that round's stores (run_etc_pipeline.py --env-from)
        reuse    every new point is in the old file: select the rows from the old stores (run_etc_pipeline.py --tracks-dir ... --reuse-env)
        extract  new points exist: extract from the catalogs (--extract-env), slow for online catalogs
  3. the ETC masks against the new track file, for a few months: every (ETC ID, time) of the mask must be a point of the track file (the mask was built
     from that file), and the ETC ID sets must be equal.
  4. the TC side, with --tc-track-file (the TC stitched-node file of the new tracking, e.g. tc_stitched_nodes.qs_filter_r15_d96.txt) and the raw TC file
     (--tc-raw-file, default: tc_stitched_nodes.txt next to it): how the TC storms changed (input | unchanged | trimmed | split | removed | output, the terms of
     Bryce's quasi-stationary logs), and the TC masks (TC_test_tracks_<source_te>_<source_res>.<YYYYMM>.nc, variable TC_int_tag) against the file: every (TC ID,
     time) of the mask must be a point of the file. The months checked include the first month with a removed point (storm IDs are running counters: after the
     first removal all later IDs shift, so a mask built from another version of the file fails there). The file NAME says little (r30_d48 files were once made
     with the wrong TC criteria): trust this check, not the name.

Exit status: 0 consistent; 1 a mask does not match its track file; 2 usage error.

Examples:
  python scripts/check_tracking_inputs.py --source icon --new-track-file /path/icon_d3hp003_hp8.etc_stitched_nodes.filtered_out_tcs.qs_filter_r30_d48.txt
  python scripts/check_tracking_inputs.py --source casesm2 --new-track-file NEW.txt --prev-track-file ROUND1/etc_tracks/casesm2_..._tcs.txt --months 202007 202101
  python scripts/check_tracking_inputs.py --source casesm2 --tc-track-file /path/casesm2_10km_nocumulus_hp8.tc_stitched_nodes.qs_filter_r15_d96.txt --time-stride 6

Author: Zhe Feng | zhe.feng@pnnl.gov
"""
import argparse
import glob
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
for p in (REPO / "scripts", REPO / "extract_environments", REPO):
    sys.path.insert(0, str(p))
from subset_etc_env_store import point_list, point_keys        # noqa: E402  (the extraction's point list and keys)

SOURCES_YAML = {"scream": "scream_ne120"}                      # ETC registry name -> key of config_sources.yaml (others are the same)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def load_all_points(path, structured=False):
    """All parsed points of a track file (with the keys), without the extraction's latitude limit."""
    import contextlib
    import io
    from src.env_extract_utilities import parse_etc_track_file
    with contextlib.redirect_stdout(io.StringIO()):
        df = parse_etc_track_file(path, unstructured_mesh=not structured)
    df["key"] = point_keys(df)
    return df


def relate(old, new):
    """How the new track file (DataFrame with 'key' and 'storm_id') relates to the old one."""
    old_keys = {}
    for i, (k, sid) in enumerate(zip(old.key, old.storm_id)):
        old_keys.setdefault(k, []).append((i, sid))
    not_in_old = sum(1 for k in new.key if k not in old_keys)
    shared = sum(1 for v in old_keys.values() if len({s for _, s in v}) > 1)
    old_pts = {sid: list(g.key) for sid, g in old.groupby("storm_id", sort=False)}
    old_sets = {sid: set(v) for sid, v in old_pts.items()}
    kinds = {"identical": 0, "trimmed": 0, "split-part": 0, "other": 0}
    children, order_kept = {}, True
    for nsid, g in new.groupby("storm_id", sort=False):
        ks = list(g.key)
        cand = None
        for c in {s for k in ks for _, s in old_keys.get(k, [])}:
            if set(ks) <= old_sets[c]:
                cand = c if cand is None else min(cand, c)
        if cand is None:
            kinds["other"] += 1
            continue
        children.setdefault(cand, []).append(nsid)
        pos = {k: i for i, k in enumerate(old_pts[cand])}
        idx = [pos[k] for k in ks]
        order_kept &= all(b > a for a, b in zip(idx, idx[1:]))
        if len(ks) == len(old_pts[cand]):
            kinds["identical"] += 1
        elif idx == list(range(idx[0], idx[0] + len(idx))) and (idx[0] == 0 or idx[-1] == len(old_pts[cand]) - 1):
            kinds["trimmed"] += 1
        else:
            kinds["split-part"] += 1
    return {"old_points": len(old), "new_points": len(new), "new_points_not_in_old": not_in_old, "old_storms": int(old.storm_id.nunique()),
            "new_storms": int(new.storm_id.nunique()), "storms": kinds, "old_storms_removed_entirely": int(old.storm_id.nunique() - len(children)),
            "old_storms_split": sum(1 for v in children.values() if len(v) > 1), "points_shared_by_2plus_storms_in_old": shared, "order_kept": bool(order_kept)}


def route(new_file, prev_file, old_df_filtered, new_df_filtered):
    """'link' / 'reuse' / 'extract' and the reason, from the files as the extraction sees them."""
    if prev_file and os.path.exists(prev_file) and md5(new_file) == md5(prev_file):
        return "link", "identical (same md5) to the previous round's track file: link that round's environment stores"
    old_keys = set(point_keys(old_df_filtered))
    missing = sum(1 for k in point_keys(new_df_filtered) if k not in old_keys)
    if missing == 0:
        return "reuse", "every point (as the extraction sees them) is in the old track file: select the rows of the old stores"
    return "extract", f"{missing} points of the new file are not in the old one: extract the environment variables from the catalogs"


def check_masks(mask_dir, basename, new_df, months, stride=1, var="ETC_int_tag"):
    """Compare the masks (variable var: ETC_int_tag or TC_int_tag) of some months with the points of the track file. Returns (ok, list of result dicts)."""
    import xarray as xr
    results, ok = [], True
    t_all = pd.to_datetime(new_df["base_time"]).values.astype("datetime64[ns]")
    for m in months:
        path = os.path.join(mask_dir, f"{basename}.{m}.nc")
        if not os.path.exists(path):
            results.append({"month": m, "file": path, "error": "file not found"})
            ok = False
            continue
        ds = xr.open_dataset(path)
        times = ds.time.values.astype("datetime64[ns]")[::stride]
        sel = (t_all >= times.min()) & (t_all <= times.max()) & np.isin(t_all, times)
        track_pairs = set(zip(new_df.storm_id.values[sel].astype(int).tolist(), t_all[sel].astype("int64").tolist()))
        mask_pairs = set()
        for i0 in range(0, len(times), 24):                               # 24 steps at a time (memory)
            idx = np.arange(i0, min(i0 + 24, len(times))) * stride
            block = ds[var].isel(time=idx).values
            for j in range(block.shape[0]):
                ids = np.unique(block[j][np.isfinite(block[j]) & (block[j] > 0)]).astype(int)
                tt = int(times[i0 + j].astype("int64"))
                mask_pairs.update((int(k), tt) for k in ids)
        ids_mask, ids_track = {p[0] for p in mask_pairs}, {p[0] for p in track_pairs}
        r = {"month": m, "mask_pairs": len(mask_pairs), "track_pairs": len(track_pairs), "in_mask_not_track": len(mask_pairs - track_pairs),
             "in_track_not_mask": len(track_pairs - mask_pairs), "ids_mask": len(ids_mask), "ids_track": len(ids_track), "ids_mask_not_track": len(ids_mask - ids_track)}
        if "TC_int_tag" in ds:
            r["tc_cells_sampled"] = int((ds["TC_int_tag"].isel(time=slice(0, None, max(1, len(times) // 8))).values > 0).sum())
        r["ok"] = r["in_mask_not_track"] == 0 and r["ids_mask_not_track"] == 0 and (r["mask_pairs"] > 0 or r["track_pairs"] == 0)   # a month without a storm is fine
        ok &= r["ok"]
        results.append(r)
    return ok, results


def storm_actions(old, new):
    """What happened to each old storm, in the terms of Bryce's logs: unchanged (all points kept, one storm), trimmed (one piece, points lost),
    split (two or more pieces), removed (no point kept). old and new: DataFrames with 'key' and 'storm_id'."""
    new_keys = {sid: set(g.key) for sid, g in new.groupby("storm_id", sort=False)}
    key_to_new = {}
    for sid, ks in new_keys.items():
        for k in ks:
            key_to_new.setdefault(k, set()).add(sid)
    counts = {"unchanged": 0, "trimmed": 0, "split": 0, "removed": 0}
    for sid, g in old.groupby("storm_id", sort=False):
        ks = list(g.key)
        kset = set(ks)
        kept = [k for k in ks if k in key_to_new]
        if not kept:
            counts["removed"] += 1
            continue
        pieces = {n for k in kept for n in key_to_new[k] if new_keys[n] <= kset}          # new storms made of points of this old storm only
        if len(pieces) >= 2:
            counts["split"] += 1
        elif len(kept) == len(ks):
            counts["unchanged"] += 1
        else:
            counts["trimmed"] += 1
    return counts


def derive_raw_tc_file(tc_file):
    """The raw TC file next to a tc_stitched_nodes.qs_filter_<criteria>.txt file (None when the name does not follow that pattern)."""
    import re
    raw = re.sub(r"\.tc_stitched_nodes\.qs_filter_[^/]*\.txt$", ".tc_stitched_nodes.txt", tc_file)
    return raw if raw != tc_file else None


def pick_tc_months(mask_dir, basename, old, new):
    """The first and the last month of the masks, and the first month with a removed point (all later IDs shift from there); at most three."""
    files = sorted(glob.glob(os.path.join(mask_dir, f"{basename}.??????.nc")))
    have = [Path(f).name.split(".")[-2] for f in files]
    if not have:
        return []
    picks = list(dict.fromkeys([have[0], have[-1]]))
    new_keys = set(new.key)
    removed = [k[0] for k in old.key if k not in new_keys]
    if removed:
        first = str(np.datetime64(int(min(removed)), "ns"))[:7].replace("-", "")
        if first in have and first not in picks:
            picks.insert(1, first)
    return picks


def resolve(source):
    """(registry name, old track file, mask dir, ETC mask basename or None, structured mesh, TC mask basename or None) of a source of the ETC registry."""
    import run_etc_pipeline as rp
    defaults, sources = rp.load_registry(str(REPO / "config" / "config_etc_pipeline.yaml"))
    name = rp.resolve_sources([source], sources)[0]
    mc, a3, _ = rp.source_info(sources[name])
    key = SOURCES_YAML.get(name, name)
    cfg = yaml.safe_load(open(REPO / "config" / "config_sources.yaml")).get(key, {})
    basename = f"ETC_test_tracks_{cfg['source_te']}_{cfg['source_res']}" if cfg.get("dir_te") else None
    tc_basename = f"TC_test_tracks_{cfg['source_te']}_{cfg['source_res']}" if cfg.get("dir_te") else None
    return name, mc["track_file"], cfg.get("dir_te"), basename, bool(mc.get("structured_mesh")), tc_basename


def pick_months(mask_dir, basename):
    files = sorted(glob.glob(os.path.join(mask_dir, f"{basename}.??????.nc")))
    months = [Path(f).name.split(".")[-2] for f in files]
    return [months[0], months[-1]] if len(months) > 1 else months


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Examples")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True, help="source of the ETC registry (name or alias: scream, icon, nicam, um, casesm2, era5)")
    ap.add_argument("--new-track-file", default=None, help="ETC stitched-node track file of the new tracking (this and/or --tc-track-file)")
    ap.add_argument("--old-track-file", default=None, help="the file the existing environment stores belong to (default: the registry's production file)")
    ap.add_argument("--prev-track-file", default=None, help="the file of the previous round (identical means the route 'link')")
    ap.add_argument("--mask-dir", default=None, help="folder of the mask netCDF (default: dir_te of config_sources.yaml)")
    ap.add_argument("--mask-basename", default=None, help="ETC_test_tracks_<source_te>_<source_res> (default: from config_sources.yaml)")
    ap.add_argument("--tc-track-file", default=None, help="TC stitched-node track file of the new tracking: checks the TC storms against the raw file and the TC masks against it")
    ap.add_argument("--tc-raw-file", default=None, help="the raw TC stitched-node file (default: tc_stitched_nodes.txt next to --tc-track-file)")
    ap.add_argument("--tc-mask-basename", default=None, help="TC_test_tracks_<source_te>_<source_res> (default: from config_sources.yaml)")
    ap.add_argument("--months", nargs="*", default=None, help="YYYYMM of the ETC masks to check (default: the first and the last month found)")
    ap.add_argument("--time-stride", type=int, default=1, help="use every n-th time step of the masks (default 1: all)")
    ap.add_argument("--radius", type=float, default=20.0, help="extraction radius (the extraction keeps |lat| <= 90 - radius), default 20")
    ap.add_argument("--json", default=None, help="also write the result to this JSON file")
    args = ap.parse_args()

    if not args.new_track_file and not args.tc_track_file:
        print("ERROR: give --new-track-file (ETC) and/or --tc-track-file (TC)")
        return 2
    try:
        name, reg_old, dir_te, basename, structured, tc_basename = resolve(args.source)
    except SystemExit as exc:
        print(exc)
        return 2
    old_file = args.old_track_file or reg_old
    files = ([("new", args.new_track_file), ("old", old_file)] if args.new_track_file else []) + ([("TC", args.tc_track_file)] if args.tc_track_file else [])
    for label, f in files:
        if not os.path.exists(f):
            print(f"ERROR: {label} track file not found: {f}")
            return 2
    mask_dir, basename = args.mask_dir or dir_te, args.mask_basename or basename
    tc_basename = args.tc_mask_basename or tc_basename
    tc_raw = None
    if args.tc_track_file:
        tc_raw = args.tc_raw_file or derive_raw_tc_file(args.tc_track_file)
        if not tc_raw or not os.path.exists(tc_raw):
            print(f"ERROR: raw TC file not found ({tc_raw}); give --tc-raw-file")
            return 2

    print(f"== {name} ==")
    ok, masks, rel, kind, why = True, [], None, None, None
    if args.new_track_file:
        print(f"new track file: {args.new_track_file}\nold track file: {old_file}" + (f"\nprevious round: {args.prev_track_file}" if args.prev_track_file else ""))
        old_all, new_all = load_all_points(old_file, structured), load_all_points(args.new_track_file, structured)
        rel = relate(old_all, new_all)
        print(f"points: old {rel['old_points']}, new {rel['new_points']} (removed {rel['old_points'] - rel['new_points']}), new points not in old: {rel['new_points_not_in_old']}, "
              f"storm order kept: {rel['order_kept']}")
        print(f"storms: old {rel['old_storms']}, new {rel['new_storms']}; new storms unchanged {rel['storms']['identical']}, trimmed at an end {rel['storms']['trimmed']}, "
              f"split-part {rel['storms']['split-part']}, other {rel['storms']['other']}; old storms removed entirely {rel['old_storms_removed_entirely']}, "
              f"split into 2+ {rel['old_storms_split']}; points shared by 2+ old storms {rel['points_shared_by_2plus_storms_in_old']}")
        old_f, new_f = point_list(old_file, args.radius, not structured), point_list(args.new_track_file, args.radius, not structured)
        kind, why = route(args.new_track_file, args.prev_track_file, old_f, new_f)
        print(f"environment route: {kind.upper()} - {why}")

        if mask_dir and basename:
            months = args.months or pick_months(mask_dir, basename)
            if months:
                ok, masks = check_masks(mask_dir, basename, new_all, months, args.time_stride)      # the masks come from the whole file
                for r in masks:
                    print(f"mask {r['month']}: " + (r["error"] if "error" in r else
                          f"(ETC id, time) pairs: mask {r['mask_pairs']}, track {r['track_pairs']}, in mask but not track {r['in_mask_not_track']}, in track but not mask "
                          f"{r['in_track_not_mask']}; ids: mask {r['ids_mask']}, track {r['ids_track']}, mask-not-track {r['ids_mask_not_track']}"
                          + (f"; TC cells (sampled steps) {r['tc_cells_sampled']}" if "tc_cells_sampled" in r else "") + f"  -> {'OK' if r['ok'] else 'MISMATCH'}"))
            else:
                print(f"masks: no files {basename}.YYYYMM.nc in {mask_dir}: not checked")
        else:
            print("masks: no mask directory for this source (ERA5: remapped masks): not checked")
    tc_ok, tc_masks, tc_rel = True, [], None
    if args.tc_track_file:
        print(f"TC track file: {args.tc_track_file}\nraw TC file:   {tc_raw}")
        tc_old, tc_new = load_all_points(tc_raw, structured), load_all_points(args.tc_track_file, structured)
        tc_rel = dict(storm_actions(tc_old, tc_new), input_storms=int(tc_old.storm_id.nunique()), output_storms=int(tc_new.storm_id.nunique()),
                      old_points=len(tc_old), new_points=len(tc_new), new_points_not_in_raw=sum(1 for k in tc_new.key if k not in set(tc_old.key)))
        print(f"TC storms vs the raw file: input {tc_rel['input_storms']} | unchanged {tc_rel['unchanged']} | trimmed {tc_rel['trimmed']} | split {tc_rel['split']} | "
              f"removed {tc_rel['removed']} | output {tc_rel['output_storms']}; points {tc_rel['old_points']} -> {tc_rel['new_points']}, new points not in the raw file: {tc_rel['new_points_not_in_raw']}")
        if mask_dir and tc_basename:
            months = pick_tc_months(mask_dir, tc_basename, tc_old, tc_new)
            if months:
                tc_ok, tc_masks = check_masks(mask_dir, tc_basename, tc_new, months, args.time_stride, var="TC_int_tag")
                for r in tc_masks:
                    print(f"TC mask {r['month']}: " + (r["error"] if "error" in r else
                          f"(TC id, time) pairs: mask {r['mask_pairs']}, file {r['track_pairs']}, in mask but not file {r['in_mask_not_track']}, in file but not mask "
                          f"{r['in_track_not_mask']}; ids: mask {r['ids_mask']}, file {r['ids_track']}, mask-not-file {r['ids_mask_not_track']}  -> {'OK' if r['ok'] else 'MISMATCH'}"))
            else:
                print(f"TC masks: no files {tc_basename}.YYYYMM.nc in {mask_dir}: not checked")
        else:
            print("TC masks: no mask directory for this source: not checked")
    ok = ok and tc_ok
    print("RESULT:", "CONSISTENT" if ok else "MISMATCH between the masks and the track file(s)")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"source": name, "relation": rel, "route": kind, "reason": why, "masks": masks, "tc": tc_rel, "tc_masks": tc_masks, "ok": ok}, f, indent=1)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
