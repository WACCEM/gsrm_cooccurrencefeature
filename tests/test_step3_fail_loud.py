"""Tests that Step 3 (make_cooccurrence_masks.py) fails loudly: frames without a result are reported and give exit status 1.

  1. all_time_steps_written: complete / count short / frames listed
  2. stream_process_to_zarr(return_missing=True) on a toy store (Dask, threads) where one frame raises in the worker function: the frame is reported,
     stays NaN, the others are written; without return_missing the old two-value return is unchanged
  3. the script exits with status 1 (it used to exit 0) when its input store cannot be read

Run:  python tests/test_step3_fail_loud.py      (or pytest tests/)      Needs the pipeline environment (xarray, zarr 2, dask).
"""
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr
import zarr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))          # make_cooccurrence_masks (zarr_tools imports it by name in its sequential path)
sys.path.insert(0, str(REPO))
import make_cooccurrence_masks as step3            # noqa: E402
from src import zarr_tools                         # noqa: E402


def test_helper():
    log = logging.getLogger("t"); log.disabled = True
    assert step3.all_time_steps_written(10, 10, [], log) is True
    assert step3.all_time_steps_written(9, 10, [], log) is False, "a short count without identified frames must fail"
    assert step3.all_time_steps_written(9, 10, ["2020-01-01T06:00:00.000000000"], log) is False
    assert step3.all_time_steps_written(10, 10, ["x"], log) is False


def test_stream_reports_the_failed_frame():
    """The real parallel path (Dask, threads only so that the patched worker function is the one that runs)."""
    from dask.distributed import Client, LocalCluster
    tmp = Path(tempfile.mkdtemp(prefix="s3ff_"))
    cluster = client = None
    try:
        n_t, n_cell, bad = 7, 12, 4
        time = np.datetime64("2020-01-01T00:00", "ns") + np.arange(n_t) * np.timedelta64(6, "h")
        ds = xr.Dataset({"idx": (("time",), np.arange(n_t))}, coords={"time": time, "cell": np.arange(n_cell)})
        in_zarr = str(tmp / "in.zarr")
        ds.to_zarr(in_zarr, mode="w")                                 # the workers read their frame from this store
        variables = ["mcs_mask", "tot_pr"]
        out = str(tmp / "out.zarr")
        zarr_tools.initialize_zarr_store(out, time, variables, ds.coords, {}, chunk_size_time=3)

        def fake(_ds, verbose=False):
            i = int(_ds["idx"])
            if i == bad:
                raise ValueError("boom")                     # a worker error: caught per frame, the frame gets no result
            res = {v: np.full(n_cell, i + 1, dtype="float32") for v in variables}
            res["etc_overlap_records"] = [{"frame": i}]
            return res

        cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=2, dashboard_address=None)
        client = Client(cluster)
        real = step3.process_single_timestep_overlaps
        step3.process_single_timestep_overlaps = fake
        try:
            logging.disable(logging.CRITICAL)
            n_ok, records, missing = zarr_tools.stream_process_to_zarr(
                ds, time, variables, out, ds.coords, {}, client=client, parallel=True, chunk_size_time=3,
                input_zarr_path=in_zarr, return_missing=True)
            zarr_tools.initialize_zarr_store(out, time, variables, ds.coords, {}, chunk_size_time=3)
            old_style = zarr_tools.stream_process_to_zarr(
                ds, time, variables, out, ds.coords, {}, client=client, parallel=True, chunk_size_time=3, input_zarr_path=in_zarr)
        finally:
            step3.process_single_timestep_overlaps = real
            logging.disable(logging.NOTSET)
        assert n_ok == n_t - 1 and len(records) == n_t - 1, (n_ok, len(records))
        assert missing == [str(time[bad])], missing
        assert len(old_style) == 2, "without return_missing the return value is unchanged"
        z = zarr.open(out, mode="r")
        data = z["tot_pr"][:]
        assert np.isnan(data[bad]).all(), "the failed frame stays NaN"
        for i in range(n_t):
            if i != bad:
                assert (data[i] == i + 1).all(), f"frame {i} must be written"
    finally:
        if client is not None:
            client.close()
        if cluster is not None:
            cluster.close()
        shutil.rmtree(tmp, ignore_errors=True)


def test_exit_status_when_input_missing():
    tmp = tempfile.mkdtemp(prefix="s3ff_root_")
    try:
        env = dict(os.environ, COF_DATA_ROOT=tmp + "/", OMP_NUM_THREADS="1")
        r = subprocess.run([sys.executable, str(REPO / "scripts" / "make_cooccurrence_masks.py"), "-c", str(REPO / "config" / "config_sources.yaml"),
                            "--source", "scream_ne120", "--no-parallel"], capture_output=True, text=True, env=env, timeout=600)
        assert "Error loading dataset" in r.stdout, r.stdout[-500:] + r.stderr[-500:]
        assert r.returncode == 1, f"exit status {r.returncode}: an unreadable input must not exit with 0"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_helper(); test_stream_reports_the_failed_frame(); test_exit_status_when_input_missing()
    print("test_step3_fail_loud: all checks passed")
