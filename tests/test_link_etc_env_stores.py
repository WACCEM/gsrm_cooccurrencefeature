"""Tests of extract_environments/link_etc_env_stores.py: what is linked, what is refused.

  - the environment stores are linked (symbolic links), the excluded variables (pr, COF masks) are not
  - a second run changes nothing
  - a hollow store (metadata purged), a store with another number of storm points, and too few stores make it fail (exit status 1)

Run:  python tests/test_link_etc_env_stores.py      (or pytest tests/)
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
SCRIPT = str(REPO / "extract_environments" / "link_etc_env_stores.py")


def make_store(path, n=6):
    g = zarr.open_group(str(path), mode="w")
    g.create_dataset("time", data=np.arange(n), chunks=(4,))
    g.create_dataset("v", data=np.ones((n, 3, 3), dtype="float32"), chunks=(4, 3, 3))
    zarr.consolidate_metadata(str(path))


def run(src, dst, *extra):
    return subprocess.run([sys.executable, SCRIPT, "--src-dir", str(src), "--dst-dir", str(dst), *extra], capture_output=True, text=True)


def test_all():
    tmp = Path(tempfile.mkdtemp(prefix="lnk_"))
    try:
        src, dst = tmp / "src", tmp / "dst"
        src.mkdir()
        for v in ("psl", "tas", "hus_850hPa", "pr", "mcs_etc_overlap_mask"):
            make_store(src / f"etc_2d_{v}_all_all.zarr")
        make_store(src / "notes.zarr")                                  # not an etc_2d store: ignored
        r = run(src, dst, "--exclude", "pr", "mcs_etc_overlap_mask", "--expected-points", "6", "--min-stores", "3")
        assert r.returncode == 0, r.stdout + r.stderr
        linked = sorted(p.name for p in dst.iterdir())
        assert linked == ["etc_2d_hus_850hPa_all_all.zarr", "etc_2d_psl_all_all.zarr", "etc_2d_tas_all_all.zarr"], linked
        assert all((dst / n).is_symlink() and (dst / n).resolve() == (src / n).resolve() for n in linked)
        r = run(src, dst, "--exclude", "pr", "mcs_etc_overlap_mask", "--expected-points", "6")          # second run: nothing new
        assert r.returncode == 0 and "linked 0 stores (3 were already linked)" in r.stdout, r.stdout
        # wrong number of storm points
        r = run(src, tmp / "dst2", "--exclude", "pr", "mcs_etc_overlap_mask", "--expected-points", "7")
        assert r.returncode == 1 and "expected 7" in r.stdout
        # too few stores
        r = run(src, tmp / "dst3", "--exclude", "pr", "mcs_etc_overlap_mask", "--min-stores", "4")
        assert r.returncode == 1 and "at least 4" in r.stdout
        # a hollow store (the purge left empty directories)
        hollow = src / "etc_2d_psl_all_all.zarr"
        for dp, dn, fn in os.walk(hollow):
            for f in fn:
                os.remove(os.path.join(dp, f))
        r = run(src, tmp / "dst4", "--exclude", "pr", "mcs_etc_overlap_mask")
        assert r.returncode == 1 and "psl" in r.stdout and not (tmp / "dst4").exists(), r.stdout
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_all()
    print("test_link_etc_env_stores: all checks passed")
