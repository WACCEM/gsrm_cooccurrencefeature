#!/usr/bin/env python
"""
Check Zarr (v2) stores for chunk files that are missing or empty.

A missing chunk is not an error when the store is read: zarr returns the fill value (NaN) for it, so a store that lost chunks
(for example to the scratch purge, or to an interrupted copy) still opens and silently contains NaN blocks. This script compares
the chunk keys that the array metadata says should exist with the files that do, key by key (not by count: a stray extra file
must not hide a missing one), and reports the missing ones and the zero-byte ones. It reads directory listings and file sizes
only, so it is quick even for large stores. Both chunk layouts are handled ('.' separated keys such as 12.3 in one directory, and
'/' separated keys as nested directories).

Usage:
  python check_zarr_store.py STORE [STORE ...] [--arrays mcs_mask tot_pr] [--quiet]

Exit status 1 when any chunk is missing or has zero bytes; unexpected extra chunk files are only reported.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import argparse
import itertools
import json
import os
import sys


def _read_json(path):
    with open(path) as f:
        return json.load(f)


def array_specs(store):
    """
    Shape, chunk shape and key separator of every array of a store.

    Reads the consolidated metadata when the store has it and the .zarray file of each top-level array otherwise.

    Returns:
    --------
    dict : {array name: {'shape': [...], 'chunks': [...], 'sep': '.' or '/'}}
    """
    specs = {}
    zmeta = os.path.join(store, ".zmetadata")
    if os.path.exists(zmeta):
        for key, meta in _read_json(zmeta)["metadata"].items():
            if key.endswith("/.zarray"):
                specs[key[:-len("/.zarray")]] = meta
    else:
        for name in sorted(os.listdir(store)):
            zarray = os.path.join(store, name, ".zarray")
            if os.path.exists(zarray):
                specs[name] = _read_json(zarray)
    return {name: {"shape": list(m["shape"]), "chunks": list(m["chunks"]), "sep": m.get("dimension_separator") or "."}
            for name, m in specs.items()}


def expected_keys(shape, chunks, sep):
    """The chunk keys an array of this shape and chunking must have ('0' for a scalar, none for an array with a zero-length dimension)."""
    if not shape:
        return {"0"}
    if 0 in shape:
        return set()
    grid = [-(-s // c) for s, c in zip(shape, chunks)]
    return {sep.join(str(i) for i in idx) for idx in itertools.product(*[range(g) for g in grid])}


def _present(path, prefix, sep, levels_left, found):
    """Collect {chunk key: size in bytes} below an array directory, following nested directories for the '/' layout."""
    for entry in os.scandir(path):
        if entry.name.startswith("."):
            continue
        key = f"{prefix}{sep}{entry.name}" if prefix else entry.name
        if sep == "/" and levels_left > 0 and entry.is_dir():
            _present(entry.path, key, sep, levels_left - 1, found)
        else:
            found[key] = entry.stat().st_size


def scan_store(store, arrays=None):
    """
    Compare the chunk files of a store with what its metadata says should exist.

    Parameters:
    -----------
    store : str
        Path of the Zarr store
    arrays : list of str, optional
        Only these arrays (default: all)

    Scalar (0-d) arrays such as 'crs' are skipped: zarr does not write a chunk that only holds the fill value, so a scalar
    that equals its fill value has no chunk file although nothing is lost.

    Returns:
    --------
    dict : {array name: {'expected': int, 'present': int, 'missing': [keys], 'zero': [keys], 'unexpected': [keys]}}
        'present' counts the expected chunks that exist; the lists are sorted.
    """
    if not os.path.isdir(store):
        raise FileNotFoundError(f"Not a directory: {store}")
    result = {}
    for name, spec in sorted(array_specs(store).items()):
        if (arrays and name not in arrays) or not spec["shape"]:
            continue
        adir = os.path.join(store, name)
        exp = expected_keys(spec["shape"], spec["chunks"], spec["sep"])
        found = {}
        if os.path.isdir(adir):
            _present(adir, "", spec["sep"], max(len(spec["shape"]) - 1, 0), found)
        result[name] = {
            "expected": len(exp),
            "present": len(exp & set(found)),
            "missing": sorted(exp - set(found)),
            "zero": sorted(k for k in exp & set(found) if found[k] == 0),
            "unexpected": sorted(set(found) - exp),
        }
    return result


def summarize(store, scan, max_list=5):
    """Lines describing a scan result; the bool says whether the store is complete (nothing missing, nothing empty)."""
    lines, complete = [], True
    for name, r in scan.items():
        problem = bool(r["missing"] or r["zero"])
        complete &= not problem
        text = f"  {name}: {r['present']}/{r['expected']} chunks"
        if r["missing"]:
            text += f", MISSING {len(r['missing'])} (e.g. {', '.join(r['missing'][:max_list])})"
        if r["zero"]:
            text += f", ZERO-BYTE {len(r['zero'])} (e.g. {', '.join(r['zero'][:max_list])})"
        if r["unexpected"]:
            text += f", {len(r['unexpected'])} unexpected extra file(s) (e.g. {', '.join(r['unexpected'][:max_list])})"
        lines.append((problem, text))
    return complete, lines


def main():
    parser = argparse.ArgumentParser(description="Check Zarr stores for missing or zero-byte chunk files")
    parser.add_argument("stores", nargs="+", help="Zarr store(s) to check")
    parser.add_argument("--arrays", nargs="*", default=None, help="only these arrays (default: all)")
    parser.add_argument("--quiet", action="store_true", help="print only the arrays with a problem")
    args = parser.parse_args()

    all_ok = True
    for store in args.stores:
        try:
            scan = scan_store(store, args.arrays)
        except Exception as exc:  # noqa: BLE001
            print(f"{store}: cannot be checked ({type(exc).__name__}: {exc})")
            all_ok = False
            continue
        complete, lines = summarize(store, scan)
        all_ok &= complete
        print(f"{store}: {'complete' if complete else 'INCOMPLETE'}")
        for problem, text in lines:
            if problem or not args.quiet:
                print(text)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
