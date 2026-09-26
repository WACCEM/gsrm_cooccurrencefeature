"""Tests of scripts/check_tracking_inputs.py (the checks made before a re-processing round).

  - how a new ETC track file relates to the old one: unchanged / trimmed at an end / split storms, storms removed entirely, IDs renumbered, points that are
    not in the old file, order kept
  - the environment route: link (identical to the previous round's file), reuse (every point is in the old file), extract (new points)
  - the ETC masks against the track file: consistent, or a mask ID / time that the track file does not have (exit status 1)

Run:  python tests/test_check_tracking_inputs.py      (or pytest tests/)
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import check_tracking_inputs as chk        # noqa: E402

SCRIPT = str(REPO / "scripts" / "check_tracking_inputs.py")
t = lambda h: (2020, 1, 1 + h // 24, h % 24)                                       # noqa: E731


def P(i, hour):
    return (100 + i, 10.0 + i, 20.0 + i, t(hour))                                    # grid_id, lon, lat, time


def write(path, storms):
    with open(path, "w") as f:
        for st in storms:
            y, m, d, h = st[0][3]
            f.write(f"start\t{len(st)}\t{y}\t{m}\t{d}\t{h}\n")
            for gid, lon, lat, (y, m, d, h) in st:
                f.write(f"\t{gid}\t{lon}\t{lat}\t1.0\t2.0\t3.0\t4.0\t{y}\t{m}\t{d}\t{h}\n")


def files(tmp):
    s1 = [P(i, 6 * i) for i in range(5)]                    # A..E, split in the new file (C removed)
    s2 = [P(10 + i, 6 * i) for i in range(3)]               # F, G, H: F removed (trimmed at the start)
    s3 = [P(20 + i, 6 * i) for i in range(2)]               # removed entirely
    s4 = [P(30 + i, 6 * i) for i in range(2)]               # unchanged
    write(tmp / "old.txt", [s1, s2, s3, s4])
    write(tmp / "new.txt", [s1[:2], s1[3:], s2[1:], s4])    # IDs 1..4 in the new file (the old storm 4 was number 4 too, but 3 disappeared)
    return s1, s2, s3, s4


def test_relation_of_the_new_file_to_the_old_one():
    tmp = Path(tempfile.mkdtemp(prefix="chk_"))
    try:
        files(tmp)
        old, new = chk.load_all_points(str(tmp / "old.txt")), chk.load_all_points(str(tmp / "new.txt"))
        r = chk.relate(old, new)
        assert (r["old_points"], r["new_points"], r["new_points_not_in_old"]) == (12, 8, 0)
        assert (r["old_storms"], r["new_storms"]) == (4, 4) and r["order_kept"]
        assert r["storms"] == {"identical": 1, "trimmed": 3, "split-part": 0, "other": 0}, r["storms"]
        assert r["old_storms_removed_entirely"] == 1 and r["old_storms_split"] == 1
        # the storm IDs are running counters: the last storm is number 4 in the old file and number 4 in the new one, but the storms in between changed
        assert list(old.groupby("storm_id").size()) == [5, 3, 2, 2] and list(new.groupby("storm_id").size()) == [2, 2, 2, 2]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_route_link_reuse_or_extract():
    tmp = Path(tempfile.mkdtemp(prefix="chk_"))
    try:
        s1, s2, s3, s4 = files(tmp)
        shutil.copy(tmp / "new.txt", tmp / "prev.txt")
        old_f, new_f = chk.point_list(str(tmp / "old.txt"), 20.0), chk.point_list(str(tmp / "new.txt"), 20.0)
        assert chk.route(str(tmp / "new.txt"), str(tmp / "prev.txt"), old_f, new_f)[0] == "link"
        assert chk.route(str(tmp / "new.txt"), None, old_f, new_f)[0] == "reuse"
        assert chk.route(str(tmp / "new.txt"), str(tmp / "old.txt"), old_f, new_f)[0] == "reuse"          # a previous file that is not identical
        write(tmp / "newer.txt", [s1[:2], [(999, 55.5, 44.0, t(30))]])                                       # a point that the old file does not have
        kind, why = chk.route(str(tmp / "newer.txt"), None, old_f, chk.point_list(str(tmp / "newer.txt"), 20.0))
        assert kind == "extract" and "1 points" in why
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def write_masks(tmp, ids_by_time, name="ETC_test_tracks_x_y.202001.nc"):
    times = np.array([np.datetime64("2020-01-01T00:00") + np.timedelta64(6 * h, "h") for h in range(len(ids_by_time))], dtype="datetime64[ns]")
    etc = np.zeros((len(times), 6), dtype="float32")
    for k, ids in enumerate(ids_by_time):
        for c, i in enumerate(ids):
            etc[k, c] = i
    xr.Dataset({"ETC_int_tag": (("time", "cell"), etc), "TC_int_tag": (("time", "cell"), np.ones_like(etc))}, coords={"time": times}).to_netcdf(tmp / name)


def test_masks_against_the_track_file():
    tmp = Path(tempfile.mkdtemp(prefix="chk_"))
    try:
        files(tmp)
        new = chk.load_all_points(str(tmp / "new.txt"))
        # new file: storm 1 = A(t0) B(t6); storm 2 = D(t18) E(t24); storm 3 = G(t6) H(t12); storm 4 = K(t0) L(t6)
        good = [[1, 4], [1, 3, 4], [3], [2], [2], []]
        write_masks(tmp, good)
        ok, res = chk.check_masks(str(tmp), "ETC_test_tracks_x_y", new, ["202001"])
        assert ok and res[0]["in_mask_not_track"] == 0 and res[0]["ids_mask"] == 4 and res[0]["tc_cells_sampled"] > 0, res
        bad_id = [[1, 4], [1, 3, 4], [3], [2], [2, 7], []]                      # ID 7 is not in the track file
        write_masks(tmp, bad_id)
        ok, res = chk.check_masks(str(tmp), "ETC_test_tracks_x_y", new, ["202001"])
        assert not ok and res[0]["in_mask_not_track"] == 1 and res[0]["ids_mask_not_track"] == 1
        wrong_time = [[1, 4], [1, 3, 4], [3], [2], [1], []]                      # storm 1 has no point at t24
        write_masks(tmp, wrong_time)
        ok, res = chk.check_masks(str(tmp), "ETC_test_tracks_x_y", new, ["202001"])
        assert not ok and res[0]["in_mask_not_track"] == 1 and res[0]["ids_mask_not_track"] == 0
        ok, res = chk.check_masks(str(tmp), "ETC_test_tracks_x_y", new, ["202002"])
        assert not ok and res[0]["error"] == "file not found"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_command_line_and_exit_status():
    tmp = Path(tempfile.mkdtemp(prefix="chk_"))
    try:
        files(tmp)
        write_masks(tmp, [[1, 4], [1, 3, 4], [3], [2], [2], []])
        args = [sys.executable, SCRIPT, "--source", "icon", "--new-track-file", str(tmp / "new.txt"), "--old-track-file", str(tmp / "old.txt"),
                "--mask-dir", str(tmp), "--mask-basename", "ETC_test_tracks_x_y", "--months", "202001", "--json", str(tmp / "out.json")]
        r = subprocess.run(args, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        assert "== icon_d3hp003 ==" in r.stdout and "environment route: REUSE" in r.stdout and "RESULT: CONSISTENT" in r.stdout, r.stdout
        assert "new points not in old: 0" in r.stdout and "old storms removed entirely 1" in r.stdout
        import json
        assert json.load(open(tmp / "out.json"))["route"] == "reuse"
        write_masks(tmp, [[1, 4], [9], [3], [2], [2], []])
        r = subprocess.run(args, capture_output=True, text=True)
        assert r.returncode == 1 and "MISMATCH" in r.stdout
        r = subprocess.run([sys.executable, SCRIPT, "--source", "nosuchsource", "--new-track-file", str(tmp / "new.txt")], capture_output=True, text=True)
        assert r.returncode == 2
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_relation_of_the_new_file_to_the_old_one(); test_route_link_reuse_or_extract(); test_masks_against_the_track_file(); test_command_line_and_exit_status()
    print("test_check_tracking_inputs: all checks passed")
