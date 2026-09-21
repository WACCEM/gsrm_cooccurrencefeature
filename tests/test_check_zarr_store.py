"""Tests of scripts/check_zarr_store.py: complete, missing chunk, hollow (purged) store, metadata lost, requested array absent.

Run:  python tests/test_check_zarr_store.py      (or pytest tests/)      Needs zarr 2 and numpy.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import zarr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from check_zarr_store import scan_store, summarize  # noqa: E402


def make_store(path, consolidate=True):
    g = zarr.open_group(str(path), mode="w")
    g.create_dataset("time", data=np.arange(10), chunks=(5,))
    g.create_dataset("pr", data=np.ones((10, 8), dtype="float32"), chunks=(5, 8))
    g.create_dataset("crs", data=np.int32(0))                      # scalar: no chunk file needed
    if consolidate:
        zarr.consolidate_metadata(str(path))


def run_cli(*stores):
    return subprocess.run([sys.executable, str(REPO / "scripts" / "check_zarr_store.py"), *map(str, stores)],
                          capture_output=True, text=True).returncode


def test_all():
    tmp = Path(tempfile.mkdtemp(prefix="czs_"))
    try:
        # 1) complete store, with and without consolidated metadata
        for cons in (True, False):
            s = tmp / f"ok_{cons}.zarr"; make_store(s, cons)
            complete, _ = summarize(str(s), scan_store(str(s)))
            assert complete, "a complete store must pass"
            assert run_cli(s) == 0
        # 2) one chunk missing -> incomplete, key reported
        s = tmp / "missing.zarr"; make_store(s)
        os.remove(s / "pr" / "1.0")
        scan = scan_store(str(s)); complete, lines = summarize(str(s), scan)
        assert not complete and scan["pr"]["missing"] == ["1.0"], scan["pr"]
        assert run_cli(s) == 1
        # 3) hollow store: the purge left only empty directories (no .zgroup, .zmetadata, .zarray, chunks)
        s = tmp / "hollow.zarr"; make_store(s)
        for dp, dn, fn in os.walk(s):
            for f in fn:
                os.remove(os.path.join(dp, f))
        assert (s / "pr").is_dir() and not list((s / "pr").iterdir())
        try:
            scan_store(str(s)); raise AssertionError("a hollow store must not pass")
        except ValueError as e:
            assert "no array metadata" in str(e)
        assert run_cli(s) == 1
        # 4) metadata lost but chunk files still there (the IR_IMERG zoom-9 case)
        s = tmp / "nometa.zarr"; make_store(s, consolidate=False)
        for name in (".zgroup", ".zattrs", ".zmetadata"):
            if (s / name).exists():
                os.remove(s / name)
        for a in ("time", "pr", "crs"):
            if (s / a / ".zarray").exists():
                os.remove(s / a / ".zarray")
        assert (s / "pr" / "0.0").exists()
        try:
            scan_store(str(s)); raise AssertionError("a store without metadata must not pass")
        except ValueError:
            pass
        assert run_cli(s) == 1
        # 5) requested array that is not in the store
        s = tmp / "ok_True.zarr"
        try:
            scan_store(str(s), ["tot_pr"]); raise AssertionError("a requested array that is absent must not pass")
        except ValueError as e:
            assert "tot_pr" in str(e)
        assert set(scan_store(str(s), ["pr"])) == {"pr"}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_all()
    print("test_check_zarr_store: all checks passed")
