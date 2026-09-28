"""Tests of scripts/combine_etc_cof_data.py: ETC track points are matched to Step 3's overlap-tracking parquet by (time, lon,
lat), never by storm_id/etc_track (a mask's own numbering, e.g. ERA5's ETC_int_tag/TC_int_tag from the full multi-decade
tracking, need not match this source's own track file's storm_id - see docs/procedures/rerun_after_tracking_update.md).

  - the ERA5 scenario: the COF parquet's etc_track uses a different (larger) numbering than the track file's storm_id; an
    id_reference_file resolves each etc_track to a position, and the points still match correctly by position
  - without id_reference_file, the same mismatched etc_track values fail to resolve (informative "N unresolved" warning,
    NaN overlap_flag) - this documents the bug the fix addresses, not a desired outcome
  - the ordinary case (etc_track already equals this source's own storm_id, no id_reference_file needed) still matches
    correctly, and the output's storm_id column is always the track file's own (never the mask's / cof_df's etc_track)

Run:  python tests/test_combine_etc_cof_data.py      (or pytest tests/)
"""
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from combine_etc_cof_data import combine_etc_cof_data      # noqa: E402

t = lambda h: (2020, 1, 1 + h // 24, h % 24)                # noqa: E731


def P(i, hour):
    return (100 + i, 10.0 + i, 20.0 + i, t(hour))            # grid_id, lon, lat, time


def write_track_file(path, storms):
    """storms: [[(grid_id, lon, lat, (y,m,d,h)), ...], ...], one list of points per storm, in file order."""
    with open(path, "w") as f:
        for st in storms:
            y, m, d, h = st[0][3]
            f.write(f"start\t{len(st)}\t{y}\t{m}\t{d}\t{h}\n")
            for gid, lon, lat, (y, m, d, h) in st:
                f.write(f"\t{gid}\t{lon}\t{lat}\t1.0\t2.0\t3.0\t4.0\t{y}\t{m}\t{d}\t{h}\n")


def write_cof_parquet(path, rows):
    """rows: [(etc_track, pd.Timestamp, overlap_flag), ...]."""
    pd.DataFrame({
        "etc_track": [r[0] for r in rows], "time": [r[1] for r in rows], "overlap_flag": [r[2] for r in rows],
        "ar_tracks": [[] for _ in rows], "mcs_tracks": [[] for _ in rows],
    }).to_parquet(path, index=False)


def timestamp(pt):
    y, m, d, h = pt[3]
    return pd.Timestamp(year=y, month=m, day=d, hour=h)


def run(tmp, storms, cof_rows, id_reference_storms=None):
    etc_file = tmp / "src.etc_stitched_nodes.filtered_out_tcs.txt"      # no "era5" in the name: parsed as unstructured
    write_track_file(etc_file, storms)
    cof_file = tmp / "src_etc_overlap_tracking.parquet"
    write_cof_parquet(cof_file, cof_rows)
    ref_file = None
    if id_reference_storms is not None:
        ref_file = tmp / "full_record.txt"
        write_track_file(ref_file, id_reference_storms)
    out_file = tmp / "out.parquet"
    df = combine_etc_cof_data(str(etc_file), str(cof_file), str(out_file), str(ref_file) if ref_file else None)
    return df, pd.read_parquet(out_file)


def test_id_reference_file_resolves_a_mismatched_numbering():
    """The ERA5 scenario: cof_df's etc_track (53101, 53102, ...) is not this file's own storm_id (1, 2) at all, but IS the
    storm_id of a reference file covering the same points - the merge must still find every point by position."""
    tmp = Path(tempfile.mkdtemp())
    try:
        storms = [[P(0, 0), P(0, 6)], [P(1, 0), P(1, 6), P(1, 12)]]                 # storm_id 1, 2 in the source's own file
        ref_storms = [[P(0, 0), P(0, 6)], [P(1, 0), P(1, 6), P(1, 12)]]             # identical points, storm_id 1, 2 too by
        # default (parse_etc_track_file numbers storms by file order) - relabel the reference file's storms to large IDs by
        # writing extra leading storms so numbering starts far from 1, mimicking "the mask's IDs are much larger"
        ref_storms = [[P(900 + i, 0)] for i in range(9998)] + ref_storms            # storms 1..9998 padding, then 9999, 10000
        df, saved = run(tmp, storms, cof_rows=[
            (9999, timestamp(P(0, 0)), 0.0), (9999, timestamp(P(0, 6)), 1.0),
            (10000, timestamp(P(1, 0)), 2.0), (10000, timestamp(P(1, 6)), 3.0), (10000, timestamp(P(1, 12)), 0.0),
        ], id_reference_storms=ref_storms)
        assert len(saved) == 5
        assert saved["overlap_flag"].notna().all(), "every point should resolve through the reference file"
        assert sorted(saved["storm_id"].unique().tolist()) == [1, 2], "storm_id is always the source's own numbering, never etc_track"
        # the overlap_flag values landed on the correct (time-matched) point, not just filled in file order
        row = saved[(saved["storm_id"] == 1) & (saved["hour"] == 6)]
        assert row["overlap_flag"].iloc[0] == 1.0
        row = saved[(saved["storm_id"] == 2) & (saved["hour"] == 12)]
        assert row["overlap_flag"].iloc[0] == 0.0
    finally:
        shutil.rmtree(tmp)


def test_without_id_reference_file_a_mismatched_numbering_fails_to_resolve():
    """Documents the bug: passing the same mismatched etc_track values without id_reference_file cannot resolve them (no
    reference file has storm 9999/10000), so overlap_flag stays NaN and the count is reported, not silently wrong data."""
    tmp = Path(tempfile.mkdtemp())
    try:
        storms = [[P(0, 0), P(0, 6)]]
        df, saved = run(tmp, storms, cof_rows=[(9999, timestamp(P(0, 0)), 0.0), (9999, timestamp(P(0, 6)), 1.0)])
        assert saved["overlap_flag"].isna().all()
    finally:
        shutil.rmtree(tmp)


def test_ordinary_case_etc_track_already_matches_this_files_own_storm_id():
    """The five GSRM sources' situation: the mask's storm ID already equals this file's own storm_id, no id_reference_file
    needed (the default) - matching by position still finds every point, and recovers points an ID-based join would miss
    whenever the mask's storms were split/renumbered slightly (see the ICON regression check in the fix's commit message)."""
    tmp = Path(tempfile.mkdtemp())
    try:
        storms = [[P(0, 0), P(0, 6)], [P(1, 0)]]                                    # storm_id 1, 2
        df, saved = run(tmp, storms, cof_rows=[
            (1, timestamp(P(0, 0)), 3.0), (1, timestamp(P(0, 6)), 2.0), (2, timestamp(P(1, 0)), 1.0),
        ])
        assert len(saved) == 3 and saved["overlap_flag"].notna().all()
        assert saved.set_index("hour")["overlap_flag"].to_dict()[6] == 2.0
    finally:
        shutil.rmtree(tmp)


def test_unmatched_points_are_kept_with_nan_and_counted():
    """A left join: every ETC track point is kept even without a COF match (fewer COF records than track points, e.g. Step
    3's own latitude limit), and the gap is printed, not silently dropped."""
    tmp = Path(tempfile.mkdtemp())
    try:
        storms = [[P(0, 0), P(0, 6), P(0, 12)]]
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            df, saved = run(tmp, storms, cof_rows=[(1, timestamp(P(0, 0)), 5.0)])   # only the first point has COF data
        assert len(saved) == 3
        assert saved["overlap_flag"].notna().sum() == 1
        assert "2 ETC points missing COF data" in buf.getvalue()
    finally:
        shutil.rmtree(tmp)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok: {name}")
    print("all tests passed")
