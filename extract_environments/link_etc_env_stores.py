#!/usr/bin/env python
"""
Link the unchanged environment-variable stores of Analysis 3 into a run's single_vars folder.

The environment variables (huss, tas, psl, winds, humidity, geopotential, ...) do not depend on the COF products or on the
precipitation, so a re-run of Analysis 3 does not extract them again. The stores of an earlier extraction are linked (symbolic
links, no copy) into <run>/etc_data/<source>/single_vars/, next to the newly extracted pr and COF-mask stores, where the combine step
(combine_etc_2d_vars.py) finds all of them. Every linked store is checked first: complete chunk by chunk (check_zarr_store.py, which
also fails on a store whose metadata was purged), and, when --expected-points is given, with exactly that number of storm points
(the point list of the current track file; the combine step then verifies that the storm times and IDs agree too).

Usage:
  python link_etc_env_stores.py --src-dir /pscratch/.../etc_data/scream/single_vars --dst-dir <run>/etc_data/scream/single_vars \\
      --exclude pr mcs_ar_etc_overlap_mask ... --expected-points 22301 --min-stores 15

Exit status 1 when a store is incomplete, has another number of points, or fewer than --min-stores stores are linked.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""
import argparse
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from check_zarr_store import scan_store, summarize  # noqa: E402


def variable_of(store_name, suffix):
    """'etc_2d_<variable>_<suffix>.zarr' -> '<variable>' (None when the name does not fit)."""
    m = re.fullmatch(rf"etc_2d_(.+)_{re.escape(suffix)}\.zarr", store_name)
    return m.group(1) if m else None


def n_points(store):
    """Number of storm points of a single-variable store (the length of its time coordinate)."""
    import zarr
    return int(zarr.open(store, mode="r")["time"].shape[0])


def main():
    p = argparse.ArgumentParser(description="Link the unchanged environment-variable stores into a run's single_vars folder")
    p.add_argument("--src-dir", required=True, help="folder with the environment stores of the earlier extraction (etc_2d_<variable>_<suffix>.zarr)")
    p.add_argument("--dst-dir", required=True, help="the run's single_vars folder (created if needed)")
    p.add_argument("--exclude", nargs="*", default=[], help="variables that are NOT linked (pr and the COF masks: they are extracted again)")
    p.add_argument("--suffix", default="all_all", help="store name suffix (default all_all)")
    p.add_argument("--expected-points", type=int, default=None, help="every linked store must have this many storm points")
    p.add_argument("--min-stores", type=int, default=1, help="fail when fewer stores than this are linked")
    a = p.parse_args()

    src, dst = Path(a.src_dir), Path(a.dst_dir)
    if not src.is_dir():
        print(f"ERROR: source folder does not exist: {src}")
        return 1
    names = sorted(n for n in os.listdir(src) if variable_of(n, a.suffix) and (src / n).is_dir())
    to_link = [n for n in names if variable_of(n, a.suffix) not in set(a.exclude)]
    print(f"{len(names)} stores in {src}; {len(names) - len(to_link)} excluded (extracted again); {len(to_link)} to link")
    problems = []
    for n in to_link:
        store = str(src / n)
        try:
            complete, lines = summarize(store, scan_store(store))
            if not complete:
                problems.append(f"{n}: INCOMPLETE: " + "; ".join(t.strip() for pr_, t in lines if pr_))
                continue
            if a.expected_points is not None and n_points(store) != a.expected_points:
                problems.append(f"{n}: {n_points(store)} storm points, expected {a.expected_points}")
        except Exception as exc:  # noqa: BLE001
            problems.append(f"{n}: cannot be checked ({type(exc).__name__}: {exc})")
    if problems:
        print("ERROR: the following stores cannot be used:")
        for pr_ in problems:
            print("  -", pr_)
        return 1
    if len(to_link) < a.min_stores:
        print(f"ERROR: only {len(to_link)} stores to link, at least {a.min_stores} expected")
        return 1
    dst.mkdir(parents=True, exist_ok=True)
    made = kept = 0
    for n in to_link:
        target, link = (src / n).resolve(), dst / n
        if os.path.lexists(link):
            if link.is_symlink() and link.resolve() == target:
                kept += 1
                continue
            print(f"ERROR: {link} exists and is not a link to {target}")
            return 1
        os.symlink(target, link)
        made += 1
    print(f"linked {made} stores ({kept} were already linked) into {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
