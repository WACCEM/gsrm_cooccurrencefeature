#!/usr/bin/env python
"""
Compare the products of two processing rounds (two data roots), exactly, source by source.

Per source, when the files exist in both roots:
  stores    mcs_masks/{src}_mcs_masks_hp8.zarr, all_masks/{src}_allmasks_hp8_v1.zarr, cof_masks/{src}_cofmasks_hp8_v1.zarr: per variable the number of
            cell-times that differ (NaN-aware) and the footprint (> 0) of the old and the new store
  netcdf    cof_masks/stats/monthly/{src}_monthly_rainmap_cof_hp8_v1.nc, extreme_precip/{src}_precip_percentiles_6h_hp8_v1.nc,
            extreme_precip/{src}_stormtype_spatial_p{90,95}.nc, cof_masks/stats/{src}_mcs_cof_tracks_2d.nc, the Analysis 3 composites and statistics
            (etc_data/stats/{a3}/etc_2d_composite_*.nc, etc_data/stats/etc_spatial_stats_{a3}.nc): per file whether every variable is identical
  parquet   cof_masks/stats/{src}_etc_overlap_tracking.parquet, {src}_mcs_cof_flags.parquet, {src}_mcs_trackstats_cof.parquet, etc_tracks/{a3}_etc_cof_data.parquet

The smoke test of a procedure (same inputs, so the products must be identical) uses --expect identical: the exit status is 1 when anything differs.
Without it the differences are only reported (a round with new inputs).

Usage:
  python scripts/compare_processing_rounds.py --old-root ROUND1 --new-root ROUND2 --sources icon_d3hp003 casesm2_10km_nocumulus [--expect identical] [--only stores netcdf parquet]

Author: Zhe Feng | zhe.feng@pnnl.gov
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

A3_NAME = {"IMERGv7": "era5"}                                    # the Analysis 3 name of a source (the others are the same)
NAN_SAFE = lambda a, b: (a == b) | (a.isnull() & b.isnull())      # noqa: E731


def compare_store(old_path, new_path, workers=24):
    """[(variable, n_differing_cell_times, footprint_old, footprint_new)] of two zarr stores with a time dimension; also (note) lines."""
    import dask
    import xarray as xr
    dask.config.set(num_workers=workers)
    o, n = xr.open_zarr(old_path), xr.open_zarr(new_path)
    rows, notes = [], []
    if o.sizes.get("time") != n.sizes.get("time") or not (o.time.values == n.time.values).all():
        notes.append(f"time axes differ ({o.sizes.get('time')} vs {n.sizes.get('time')} steps)")
    for v in sorted(set(o.data_vars) | set(n.data_vars)):
        if v not in o.data_vars or v not in n.data_vars:
            notes.append(f"variable {v} only in the {'old' if v in o.data_vars else 'new'} store")
            continue
        a, b = o[v], n[v]
        if "time" not in a.dims or a.shape != b.shape:
            continue
        diff = ((a != b) & ~(a.isnull() & b.isnull())).sum()
        fa, fb = (a > 0).sum(), (b > 0).sum()
        d, ca, cb = dask.compute(diff, fa, fb)
        rows.append((v, int(d), int(ca), int(cb)))
    return rows, notes


def compare_netcdf(old_path, new_path):
    """(identical, [variables that differ or exist in one file only])."""
    import xarray as xr
    o, n = xr.open_dataset(old_path), xr.open_dataset(new_path)
    bad = sorted(set(o.data_vars) ^ set(n.data_vars))
    for v in sorted(set(o.data_vars) & set(n.data_vars)):
        a, b = o[v], n[v]
        if a.shape != b.shape or not bool(NAN_SAFE(a, b).all()):
            bad.append(v)
    o.close(), n.close()
    return not bad, bad


def compare_parquet(old_path, new_path):
    """(identical, description)."""
    a, b = pd.read_parquet(old_path), pd.read_parquet(new_path)
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        return False, f"shape {a.shape} vs {b.shape}"
    if a.equals(b):
        return True, ""
    cols = [c for c in a.columns if not a[c].equals(b[c])]
    return False, f"columns that differ: {cols[:8]}"


def products(old_root, new_root, src):
    """(kind, relative path) of the products of a source that exist in both roots."""
    a3 = A3_NAME.get(src, src)
    rels = [("stores", f"mcs_masks/{src}_mcs_masks_hp8.zarr"), ("stores", f"all_masks/{src}_allmasks_hp8_v1.zarr"),
            ("stores", f"cof_masks/{src}_cofmasks_hp8_v1.zarr"),
            ("netcdf", f"cof_masks/stats/monthly/{src}_monthly_rainmap_cof_hp8_v1.nc"),
            ("netcdf", f"extreme_precip/{src}_precip_percentiles_6h_hp8_v1.nc"),
            ("netcdf", f"extreme_precip/{src}_stormtype_spatial_p90.nc"), ("netcdf", f"extreme_precip/{src}_stormtype_spatial_p95.nc"),
            ("netcdf", f"cof_masks/stats/{src}_mcs_cof_tracks_2d.nc"),
            ("netcdf", f"etc_data/stats/etc_spatial_stats_{a3}.nc"),
            ("parquet", f"cof_masks/stats/{src}_etc_overlap_tracking.parquet"), ("parquet", f"cof_masks/stats/{src}_mcs_cof_flags.parquet"),
            ("parquet", f"cof_masks/stats/{src}_mcs_trackstats_cof.parquet"), ("parquet", f"etc_tracks/{a3}_etc_cof_data.parquet")]
    for pat in sorted(glob.glob(os.path.join(new_root, f"etc_data/stats/{a3}/etc_2d_composite_*.nc"))):
        rels.append(("netcdf", os.path.relpath(pat, new_root)))
    out = []
    for kind, rel in rels:
        o, n = os.path.join(old_root, rel), os.path.join(new_root, rel)
        if os.path.exists(o) and os.path.exists(n):
            out.append((kind, rel))
        elif os.path.exists(o) != os.path.exists(n):
            out.append(("missing", rel + (" (only in the old root)" if os.path.exists(o) else " (only in the new root)")))
    return out


def run(old_root, new_root, sources, only=("stores", "netcdf", "parquet"), workers=24, out=print):
    """Compare and print; returns the number of differing items (missing products count as differing)."""
    n_diff = 0
    for src in sources:
        out(f"\n===== {src}")
        for kind, rel in products(old_root, new_root, src):
            o, n = os.path.join(old_root, rel), os.path.join(new_root, rel)
            if kind == "missing":
                out(f"  MISSING   {rel}")
                n_diff += 1
            elif kind == "stores" and "stores" in only:
                rows, notes = compare_store(o, n, workers)
                bad = [r for r in rows if r[1] > 0]
                out(f"  {'IDENTICAL' if not bad and not notes else 'DIFFERS  '} {rel}: {len(rows)} variables, {len(bad)} differ")
                for v, d, ca, cb in bad:
                    out(f"      {v:30s} cell-times differing {d:>12,d}  footprint old {ca:>12,d} new {cb:>12,d} ({100 * (cb / max(ca, 1) - 1):+.2f}%)")
                for note in notes:
                    out(f"      {note}")
                n_diff += 1 if (bad or notes) else 0
            elif kind == "netcdf" and "netcdf" in only:
                same, bad = compare_netcdf(o, n)
                out(f"  {'IDENTICAL' if same else 'DIFFERS  '} {rel}" + ("" if same else f": {len(bad)} variables differ: {bad[:6]}"))
                n_diff += 0 if same else 1
            elif kind == "parquet" and "parquet" in only:
                same, why = compare_parquet(o, n)
                out(f"  {'IDENTICAL' if same else 'DIFFERS  '} {rel}" + ("" if same else f": {why}"))
                n_diff += 0 if same else 1
    return n_diff


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-root", required=True)
    ap.add_argument("--new-root", required=True)
    ap.add_argument("--sources", nargs="+", required=True, help="source names of the products (scream, icon_d3hp003, nicam_gl11, um_glm_n2560_RAL3p3, casesm2_10km_nocumulus, IMERGv7)")
    ap.add_argument("--only", nargs="*", default=["stores", "netcdf", "parquet"], choices=["stores", "netcdf", "parquet"])
    ap.add_argument("--expect", choices=["identical"], default=None, help="exit status 1 when anything differs")
    ap.add_argument("--workers", type=int, default=24, help="dask workers for the stores")
    args = ap.parse_args()
    n = run(args.old_root.rstrip("/") + "/", args.new_root.rstrip("/") + "/", args.sources, tuple(args.only), args.workers)
    print(f"\nRESULT: {n} product(s) differ or are missing")
    return 1 if (args.expect == "identical" and n) else 0


if __name__ == "__main__":
    sys.exit(main())
