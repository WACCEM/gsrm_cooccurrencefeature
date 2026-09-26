"""Tests of the path options that let the IMERG/ERA5 chain run inside a test data root (a re-run must not touch the production tree):

  - remap_era5_masks_healpix.py: --in-dir / --out-dir; the defaults are the production paths, unchanged
  - make_mcs_swath_masks.py: --tc-source-zarr replaces the config's tc_source_zarr (default: the config's)

Run:  python tests/test_era5_path_options.py      (or pytest tests/)
"""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import remap_era5_masks_healpix as remap        # noqa: E402


def test_remap_defaults_are_the_production_paths_and_can_be_overridden():
    a = remap.parse_args([])
    assert a.in_dir == "/pscratch/sd/b/beharrop/kmscale_hackathon/hackathon_pre/era5_tracking_etc_nocoldcoreonly/"
    assert a.out_dir == "/pscratch/sd/w/wcmca1/hackathon/all_masks/"
    b = remap.parse_args(["--in-dir", "/in/x", "--out-dir", "/test/root/all_masks"])
    assert (b.in_dir, b.out_dir) == ("/in/x", "/test/root/all_masks")
    src = (REPO / "scripts" / "remap_era5_masks_healpix.py").read_text()
    assert 'out_dir = args.out_dir' in src and 'in_dir = args.in_dir' in src, "main() must use the options"


def test_step1_takes_a_tc_source_zarr_override():
    out = subprocess.run([sys.executable, str(REPO / "scripts" / "make_mcs_swath_masks.py"), "--help"], capture_output=True, text=True).stdout
    assert "--tc-source-zarr" in out
    src = (REPO / "scripts" / "make_mcs_swath_masks.py").read_text()
    assert "'tc_source_zarr': args.tc_source_zarr or config.get('tc_source_zarr')" in src, "the override must win over the config, the config remain the default"


if __name__ == "__main__":
    test_remap_defaults_are_the_production_paths_and_can_be_overridden(); test_step1_takes_a_tc_source_zarr_override()
    print("test_era5_path_options: all checks passed")
