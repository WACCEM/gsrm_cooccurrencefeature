"""Tests of extract_environments/subset_etc_env_store.py (the environment store of a new ETC track file, selected from the old store).

  - the new store has the rows of the points that are still in the new track file, in the new file's order, with the NEW storm IDs, and the layout that
    save_to_zarr() gives (time, storm/grid metadata, x/y coordinates, attributes); a point shared by two old storms uses one row of the old store
  - points whose track time has no frame in the given time axis: NaN slabs, counted, listed in the attributes, warned about; above --max-missing-fraction
    the run is refused; --no-frame-check keeps the old values
  - structured meshes (lon_id/lat_id, ERA5) work the same way
  - refused (exit status 1, nothing left behind): a new point that is not in the old file, an old store that does not belong to the old track file,
    an existing destination

Run:  python tests/test_subset_etc_env_store.py      (or pytest tests/)
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "extract_environments"))
sys.path.insert(0, str(REPO))
import subset_etc_env_store as sub                     # noqa: E402
from extract_etc_2d_vars import save_to_zarr           # noqa: E402

SCRIPT = str(REPO / "extract_environments" / "subset_etc_env_store.py")
RADIUS = 0.25                                          # 3 x 3 points around a storm


def write_track_file(path, storms):
    """storms: list of lists of (grid_id, lon, lat, (Y, M, D, H))."""
    with open(path, "w") as f:
        for pts in storms:
            y, m, d, h = pts[0][3]
            f.write(f"start\t{len(pts)}\t{y}\t{m}\t{d}\t{h}\n")
            for gid, lon, lat, (y, m, d, h) in pts:
                f.write(f"\t{gid}\t{lon}\t{lat}\t1.0\t2.0\t3.0\t4.0\t{y}\t{m}\t{d}\t{h}\n")


def make_old_store(path, track_file):
    """A store as the extraction writes it, for the points of track_file; the values identify the row (row number + 0.5 * y + 0.01 * x)."""
    df = sub.point_list(track_file, RADIUS)
    n = len(df)
    yy, xx = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
    data = (np.arange(n)[:, None, None] + 0.5 * yy + 0.01 * xx).astype("float32")
    save_to_zarr(data, np.array(pd.to_datetime(df["base_time"])), df["storm_id"].values, df["grid_id"].values, df["lat"].values, df["lon"].values,
                 np.array([-1, 0, 1]), np.array([-1, 0, 1]), "tas", str(path)[:-5], RADIUS, 0.25, 0.25, chunk_size=4, unstructured_mesh=True,
                 var_attrs={"units": "K", "n_points": n})
    return df


def setup(tmp):
    A, B, C, D, E = ((100 + i, 10.0 + i, 20.0 + i) for i in range(5))
    t = lambda h: (2020, 1, 1 + h // 24, h % 24)                          # noqa: E731
    # old file: storm 1 = 5 points, storm 2 = 3 points, the last point of storm 2 is the same point as the first of storm 1 (shared, one time)
    s1 = [(g, lo, la, t(6 * i)) for i, (g, lo, la) in enumerate((A, B, C, D, E))]
    s2 = [(200, 50.0, -30.0, t(0)), (201, 51.0, -31.0, t(6)), (A[0], A[1], A[2], t(0))]
    write_track_file(tmp / "old.txt", [s1, s2])
    # new file: the quasi-stationary middle of storm 1 (C) is removed, so storm 1 is split in two; storm 2 lost its first point; IDs are renumbered
    n1, n2, n3 = [s1[0], s1[1]], [s1[3], s1[4]], [s2[1], s2[2]]
    write_track_file(tmp / "new.txt", [n1, n2, n3])
    return s1, s2


def run(*args):
    return subprocess.run([sys.executable, SCRIPT, *map(str, args)], capture_output=True, text=True)


def test_selects_rows_of_the_new_points_with_new_ids():
    tmp = Path(tempfile.mkdtemp(prefix="subset_"))
    try:
        setup(tmp)
        old_df = make_old_store(tmp / "old_tas.zarr", tmp / "old.txt")
        assert len(old_df) == 8
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "new_tas.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new.txt")
        assert r.returncode == 0, r.stdout + r.stderr
        new, old = xr.open_zarr(tmp / "new_tas.zarr"), xr.open_zarr(tmp / "old_tas.zarr")
        # new points: A, B (storm 1), D, E (storm 2), s2[1] (storm 3), s2[2] = A again (storm 3); old rows: A 0, B 1, D 3, E 4, s2[1] 6, A (first occurrence) 0
        assert dict(new.sizes) == {"time": 6, "y": 3, "x": 3}
        assert list(new.storm_id.values) == [1, 1, 2, 2, 3, 3], "the storm IDs are those of the new file"
        rows = [0, 1, 3, 4, 6, 0]
        assert (new["tas"].values == old["tas"].values[rows]).all()
        assert (new.grid_id.values == [100, 101, 103, 104, 201, 100]).all()
        assert np.allclose(new.storm_lon.values, [10, 11, 13, 14, 51, 10]) and np.allclose(new.storm_lat.values, [20, 21, 23, 24, -31, 20])
        assert (new.time.values == pd.to_datetime(sub.point_list(tmp / "new.txt", RADIUS)["base_time"]).values).all()
        assert list(new.x.values) == [-1, 0, 1] and new.attrs["radius"] == RADIUS and new["tas"].attrs["units"] == "K"
        assert new["tas"].attrs["n_points"] == 6 and "not checked" in new["tas"].attrs["time_match"], "no catalog given: the old values are kept and it says so"
        assert new.attrs["reused_n_points"] == 6 and new.attrs["reused_from"].endswith("old_tas.zarr")
        assert new["tas"].encoding["chunks"][1:] == (3, 3)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_points_without_a_frame_become_nan_and_are_recorded():
    tmp = Path(tempfile.mkdtemp(prefix="subset_"))
    try:
        setup(tmp)
        make_old_store(tmp / "old_tas.zarr", tmp / "old.txt")
        new_df = sub.point_list(tmp / "new.txt", RADIUS)                       # times: A 00, B 06 (storm 1); D 18, E 24 (storm 2); s2[1] 06, A 00 (storm 3)
        times = pd.to_datetime(new_df["base_time"])
        # a frame axis without 2020-01-01 06:00 (rows 1 and 4 of the new file) and without 2020-01-02 00:00 (row 3)
        axis = [t for t in pd.date_range("2020-01-01T00:00", "2020-01-02T18:00", freq="6h") if t not in (pd.Timestamp("2020-01-01T06:00"), pd.Timestamp("2020-01-02T00:00"))]
        (tmp / "axis.txt").write_text("\n".join(str(t) for t in axis))
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "n.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new.txt",
                "--time-axis-file", tmp / "axis.txt", "--max-missing-fraction", 0.6)
        assert r.returncode == 0, r.stdout + r.stderr
        assert "WARNING" in r.stdout and "3 of 6" in r.stdout and "2020-01-01 06" in r.stdout
        new, old = xr.open_zarr(tmp / "n.zarr"), xr.open_zarr(tmp / "old_tas.zarr")
        gap = np.array([t in (pd.Timestamp("2020-01-01T06:00"), pd.Timestamp("2020-01-02T00:00")) for t in times])
        assert gap.sum() == 3 and np.isnan(new["tas"].values[gap]).all() and not np.isnan(new["tas"].values[~gap]).any(), "NaN slabs exactly where no frame exists"
        rows = np.array([0, 1, 3, 4, 6, 0])
        assert (new["tas"].values[~gap] == old["tas"].values[rows][~gap]).all(), "the other rows keep the old values"
        a = new["tas"].attrs
        assert a["n_points_time_missing"] == 3 and a["n_missing_track_times"] == 2 and a["n_points"] == 6 and "exact" in a["time_match"]
        assert set(map(pd.Timestamp, a["missing_track_times"])) == {pd.Timestamp("2020-01-01T06:00"), pd.Timestamp("2020-01-02T00:00")}
        # more gaps than the limit: refused, nothing left
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "m.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new.txt",
                "--time-axis-file", tmp / "axis.txt", "--max-missing-fraction", 0.4)
        assert r.returncode == 1 and "ERROR" in r.stdout and not (tmp / "m.zarr").exists()
        # --no-frame-check keeps the old values
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "k.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new.txt",
                "--time-axis-file", tmp / "axis.txt", "--no-frame-check")
        assert r.returncode == 0 and not np.isnan(xr.open_zarr(tmp / "k.zarr")["tas"].values).any()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_structured_mesh_store():
    """ERA5-like: lon_id/lat_id instead of grid_id (track file and store)."""
    tmp = Path(tempfile.mkdtemp(prefix="subset_"))
    try:
        t = lambda h: (2020, 1, 1 + h // 24, h % 24)                          # noqa: E731
        pts = [(50 + i, 60 + i, 10.0 + i, 20.0 + i, t(6 * i)) for i in range(5)]
        def write(path, storms):
            with open(path, "w") as f:
                for st in storms:
                    y, m, d, h = st[0][4]
                    f.write(f"start\t{len(st)}\t{y}\t{m}\t{d}\t{h}\n")
                    for lonid, latid, lon, lat, (y, m, d, h) in st:
                        f.write(f"\t{lonid}\t{latid}\t{lon}\t{lat}\t1.0\t2.0\t3.0\t4.0\t{y}\t{m}\t{d}\t{h}\n")
        write(tmp / "old.txt", [pts[:3], pts[3:]])
        write(tmp / "new.txt", [pts[:2], pts[4:]])                            # the third point and the fourth were removed
        old_df = sub.point_list(tmp / "old.txt", RADIUS, unstructured=False)
        assert list(old_df.columns[:4]) == ["storm_id", "num_timesteps", "lon_id", "lat_id"]
        yy, xx = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
        data = (np.arange(len(old_df))[:, None, None] + 0.5 * yy + 0.01 * xx).astype("float32")
        save_to_zarr(data, np.array(pd.to_datetime(old_df["base_time"])), old_df["storm_id"].values, list(zip(old_df["lon_id"], old_df["lat_id"])),
                     old_df["lat"].values, old_df["lon"].values, np.array([-1, 0, 1]), np.array([-1, 0, 1]), "tas", str(tmp / "old_tas"), RADIUS, 0.25, 0.25,
                     chunk_size=4, unstructured_mesh=False)
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "new_tas.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new.txt")
        assert r.returncode == 0, r.stdout + r.stderr
        new = xr.open_zarr(tmp / "new_tas.zarr")
        assert list(new.lon_id.values) == [50, 51, 54] and list(new.lat_id.values) == [60, 61, 64] and list(new.storm_id.values) == [1, 1, 2]
        assert (new["tas"].values == xr.open_zarr(tmp / "old_tas.zarr")["tas"].values[[0, 1, 4]]).all()
        assert "grid_id" not in new.data_vars
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_refuses_and_leaves_nothing():
    tmp = Path(tempfile.mkdtemp(prefix="subset_"))
    try:
        s1, s2 = setup(tmp)
        make_old_store(tmp / "old_tas.zarr", tmp / "old.txt")
        # a new point that is not in the old file
        write_track_file(tmp / "new_bad.txt", [[s1[0], (300, 70.0, 10.0, (2020, 1, 9, 0))]])
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "x.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new_bad.txt")
        assert r.returncode == 1 and "not in the old track file" in r.stdout and not (tmp / "x.zarr").exists()
        # an old store that belongs to another track file (other number of points)
        write_track_file(tmp / "other.txt", [s1])
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "y.zarr", "--old-track-file", tmp / "other.txt", "--new-track-file", tmp / "new.txt")
        assert r.returncode == 1 and "old track file gives" in r.stdout and not (tmp / "y.zarr").exists()
        # an existing destination is never overwritten
        (tmp / "z.zarr").mkdir()
        r = run("--src-store", tmp / "old_tas.zarr", "--dst-store", tmp / "z.zarr", "--old-track-file", tmp / "old.txt", "--new-track-file", tmp / "new.txt")
        assert r.returncode == 1 and "exists" in r.stdout
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_match_points_uses_time_and_position():
    tmp = Path(tempfile.mkdtemp(prefix="subset_"))
    try:
        s1, s2 = setup(tmp)
        old, new = sub.point_list(tmp / "old.txt", RADIUS), sub.point_list(tmp / "new.txt", RADIUS)
        assert list(sub.match_points(old, new)) == [0, 1, 3, 4, 6, 0]
        # the same position at another time is another point
        shifted = [[(g, lo, la, (2020, 1, 5, h)) for (g, lo, la, (y, m, d, h)) in s1]]
        write_track_file(tmp / "shifted.txt", shifted)
        try:
            sub.match_points(old, sub.point_list(tmp / "shifted.txt", RADIUS))
            raise AssertionError("shifted times must not match")
        except ValueError as exc:
            assert "not in the old track file" in str(exc)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_selects_rows_of_the_new_points_with_new_ids(); test_points_without_a_frame_become_nan_and_are_recorded(); test_structured_mesh_store()
    test_refuses_and_leaves_nothing(); test_match_points_uses_time_and_position()
    print("test_subset_etc_env_store: all checks passed")
