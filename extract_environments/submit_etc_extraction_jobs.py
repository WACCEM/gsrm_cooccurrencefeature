#!/usr/bin/env python3
"""
Submit SLURM job arrays for ETC 2D variable extraction across multiple models.

Each model's variables are grouped by catalog settings (model, params, pressure level).
Each group becomes one sbatch array job where every task processes a single variable.

Usage:
    python submit_etc_extraction_jobs.py --source ERA5
    python submit_etc_extraction_jobs.py --source SCREAM
    python submit_etc_extraction_jobs.py --source NICAM
    python submit_etc_extraction_jobs.py --source CASESM2
    python submit_etc_extraction_jobs.py --source ICON
    python submit_etc_extraction_jobs.py --source UM
    python submit_etc_extraction_jobs.py --source ERA5 --dry-run   # Print commands without submitting
    python submit_etc_extraction_jobs.py --source ERA5 --group 2   # Submit only group index 2

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import argparse
import subprocess
import os
import sys
import json
import tempfile
from pathlib import Path

# ============================================================
# SHARED CONSTANTS
# ============================================================

SCRIPT_DIR = Path(__file__).parent

# Path to the extraction Python script (relative to this file)
EXTRACT_SCRIPT = SCRIPT_DIR / "extract_etc_2d_vars.py"

# SLURM account and email
SLURM_ACCOUNT = "m1867"
MAIL_USER = "zhe.feng@pnnl.gov"

# Common extraction parameters (same for all models)
RADIUS = "20.0"
LON_RES = "0.25"
LAT_RES = "0.25"
CHUNK_SIZE = "1000"
PROGRESS_FREQ = "1000"

# COF mask variables (same for all sources)
COF_MASK_VARS = [
    "mcs_ar_etc_overlap_mask",
    "ar_mcs_etc_overlap_mask",
    "etc_mcs_ar_overlap_mask",
    "etc_ar_overlap_mask",
    "etc_mcs_overlap_mask",
    "mcs_etc_overlap_mask",
    "ar_etc_overlap_mask",
]

# ============================================================
# MODEL CONFIGURATIONS
# Each model has:
#   track_file        : path to ETC track file
#   output_dir        : where extracted zarr files go
#   catalog_url       : intake catalog URL or local path
#   current_location  : catalog location key (empty string = None)
#   structured_mesh   : True for ERA5 lat/lon grid, False for HEALPix
#   walltime          : SLURM wallclock limit
#   job_name          : short SLURM job name
#   log_prefix        : prefix for SLURM log filenames
#   job_groups        : list of variable group dicts (see below)
#
# Each job_group dict has:
#   label             : human-readable description (for logs)
#   catalog_model     : model key in the intake catalog
#   catalog_params    : dict of catalog call parameters
#   variables         : list of variable names (one task per variable)
#   pressure_levels   : pressure level string, e.g. "850" (None = 2D)
#   cof_mask          : True = pass --cof_mask flag
#   convert_wa_to_omega : True = pass --convert_wa_to_omega flag
#   catalog_url_override     : optional per-group catalog URL override
#   current_location_override: optional per-group catalog location override
# ============================================================

MODEL_CONFIGS = {

    # ----------------------------------------------------------
    "ERA5": {
        "track_file": "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/era5.etc_stitched_nodes.txt",
        "output_dir": "/pscratch/sd/w/wcmca1/hackathon/etc_data/era5/single_vars/",
        "catalog_url": "/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml",
        "current_location": "",          # empty → no --current_location argument
        "structured_mesh": True,         # ERA5 is on lat/lon grid
        "walltime": "00:45:00",
        "job_name": "era5",
        "log_prefix": "extract_etc_era5",
        "job_groups": [
            {
                "label": "2D variables",
                "catalog_model": "era5_3h",
                "catalog_params": {"zoom": 8},
                "variables": ["pr", "ps", "psl", "tas", "d2m", "u10", "v10", "sstk", "prw"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "3D variables at 850 hPa",
                "catalog_model": "era5_3h",
                "catalog_params": {"zoom": 8},
                "variables": ["ua", "va", "omega", "hus", "hur", "zg"],
                "pressure_levels": "850",
                "cof_mask": False,
            },
            {
                "label": "3D variables at 500 hPa",
                "catalog_model": "era5_3h",
                "catalog_params": {"zoom": 8},
                "variables": ["ua", "va", "omega", "hus", "hur", "zg"],
                "pressure_levels": "500",
                "cof_mask": False,
            },
            {
                "label": "COF masks",
                "catalog_model": "era5_3h",
                "catalog_params": {"zoom": 8},
                "variables": COF_MASK_VARS,
                "pressure_levels": None,
                "cof_mask": True,
            },
        ],
    },

    # ----------------------------------------------------------
    "SCREAM": {
        "track_file": "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/screamv2_ne120_hp8.etc_stitched_nodes.txt",
        "output_dir": "/pscratch/sd/w/wcmca1/hackathon/etc_data/scream/single_vars/",
        "catalog_url": "/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml",
        "current_location": "NERSC",
        "structured_mesh": False,
        "walltime": "00:15:00",
        "job_name": "scream",
        "log_prefix": "extract_etc_SCREAM",
        "job_groups": [
            {
                "label": "pr (hourly)",
                "catalog_model": "scream2D_hrly",
                "catalog_params": {"zoom": 8},
                "variables": ["pr"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "2D instantaneous surface variables",
                "catalog_model": "scream_ne120_inst",
                "catalog_params": {"zoom": 8},
                "variables": ["huss", "tas", "uas", "vas", "psl", "ps"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "2D derived pressure-level variables",
                "catalog_model": "scream_ne120",
                "catalog_params": {"zoom": 8},
                "variables": ["ua850", "va850", "ua500", "va500", "rh850", "uivt", "vivt", "zg500"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "3D variables at 850 hPa",
                "catalog_model": "scream_ne120",
                "catalog_params": {"zoom": 8},
                "variables": ["ua", "va", "omega", "hus"],
                "pressure_levels": "850",
                "cof_mask": False,
            },
            {
                "label": "3D variables at 500 hPa",
                "catalog_model": "scream_ne120",
                "catalog_params": {"zoom": 8},
                "variables": ["ua", "va", "omega", "hus"],
                "pressure_levels": "500",
                "cof_mask": False,
            },
            {
                "label": "COF masks",
                "catalog_model": "scream_ne120",
                "catalog_params": {"zoom": 8},
                "variables": COF_MASK_VARS,
                "pressure_levels": None,
                "cof_mask": True,
            },
        ],
    },

    # ----------------------------------------------------------
    "NICAM": {
        "track_file": "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/nicam_gl11_hp8.etc_stitched_nodes.txt",
        "output_dir": "/pscratch/sd/w/wcmca1/hackathon/etc_data/nicam_gl11/single_vars/",
        "catalog_url": "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
        "current_location": "NERSC",
        "structured_mesh": False,
        "walltime": "00:15:00",
        "job_name": "nicam",
        "log_prefix": "extract_etc_NICAM",
        "job_groups": [
            {
                "label": "pr (time-shifted, local catalog)",
                # pr uses a local catalog with time-shifted NICAM data
                "catalog_url_override": "/global/homes/f/feng045/program/hackathon/catalog/NERSC/main.yaml",
                "current_location_override": "",
                "catalog_model": "nicam_gl11_shifted",
                "catalog_params": {"zoom": 8},
                "variables": ["pr"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "2D variables",
                "catalog_model": "nicam_gl11",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["tas", "huss", "ps", "psl", "uas", "vas", "prw"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "3D variables at 850 hPa (excl. wa)",
                "catalog_model": "nicam_gl11",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "850",
                "cof_mask": False,
            },
            {
                "label": "3D variables at 500 hPa (excl. wa)",
                "catalog_model": "nicam_gl11",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "500",
                "cof_mask": False,
            },
            {
                "label": "wa at 850 hPa (→ omega)",
                "catalog_model": "nicam_gl11",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["wa"],
                "pressure_levels": "850",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "wa at 500 hPa (→ omega)",
                "catalog_model": "nicam_gl11",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["wa"],
                "pressure_levels": "500",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "COF masks",
                "catalog_model": "nicam_gl11",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": COF_MASK_VARS,
                "pressure_levels": None,
                "cof_mask": True,
            },
        ],
    },

    # ----------------------------------------------------------
    "CASESM2": {
        "track_file": "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/casesm2_10km_nocumulus_hp8.etc_stitched_nodes.txt",
        "output_dir": "/pscratch/sd/w/wcmca1/hackathon/etc_data/casesm2_10km_nocumulus/single_vars/",
        "catalog_url": "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
        "current_location": "online",
        "structured_mesh": False,
        "walltime": "02:00:00",   # CASESM2 data is from online, needs more time
        "job_name": "casesm2_10km_nocumulus",
        "log_prefix": "extract_etc_CASESM2",
        "job_groups": [
            {
                "label": "pr (hourly)",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8},
                "variables": ["pr"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "2D variables",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["tas", "huss", "ps", "psl", "uas", "vas", "prw"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "3D variables at 850 hPa (excl. wa)",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "850",
                "cof_mask": False,
            },
            {
                "label": "3D variables at 500 hPa (excl. wa)",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "500",
                "cof_mask": False,
            },
            {
                "label": "wa at 850 hPa (→ omega)",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["wa"],
                "pressure_levels": "850",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "wa at 500 hPa (→ omega)",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8, "time": "PT6H"},
                "variables": ["wa"],
                "pressure_levels": "500",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "COF masks",
                "catalog_model": "casesm2_10km_nocumulus",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": COF_MASK_VARS,
                "pressure_levels": None,
                "cof_mask": True,
            },
        ],
    },

    # ----------------------------------------------------------
    "ICON": {
        "track_file": "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/icon_d3hp003_hp8.etc_stitched_nodes.txt",
        "output_dir": "/pscratch/sd/w/wcmca1/hackathon/etc_data/icon_d3hp003/single_vars/",
        "catalog_url": "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
        "current_location": "NERSC",
        "structured_mesh": False,
        "walltime": "00:15:00",
        "job_name": "icon",
        "log_prefix": "extract_etc_ICON",
        "job_groups": [
            {
                "label": "pr (hourly inst)",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8},
                "variables": ["pr"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "2D variables (3h mean)",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8, "time": "PT3H", "time_method": "mean"},
                "variables": ["tas", "huss", "ps", "psl", "uas", "vas", "prw"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "3D variables at 850 hPa (excl. wa)",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8, "time": "PT6H", "time_method": "inst"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "850",
                "cof_mask": False,
            },
            {
                "label": "3D variables at 500 hPa (excl. wa)",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8, "time": "PT6H", "time_method": "inst"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "500",
                "cof_mask": False,
            },
            {
                "label": "wa at 850 hPa (→ omega)",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8, "time": "PT6H", "time_method": "inst"},
                "variables": ["wa"],
                "pressure_levels": "850",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "wa at 500 hPa (→ omega)",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8, "time": "PT6H", "time_method": "inst"},
                "variables": ["wa"],
                "pressure_levels": "500",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "COF masks",
                "catalog_model": "icon_d3hp003",
                "catalog_params": {"zoom": 8, "time": "PT3H", "time_method": "mean"},
                "variables": COF_MASK_VARS,
                "pressure_levels": None,
                "cof_mask": True,
            },
        ],
    },

    # ----------------------------------------------------------
    "UM": {
        "track_file": "/pscratch/sd/w/wcmca1/hackathon/etc_tracks/um_glm_n2560_RAL3p3_hp8.etc_stitched_nodes.txt",
        "output_dir": "/pscratch/sd/w/wcmca1/hackathon/etc_data/um_glm_n2560_RAL3p3/single_vars/",
        "catalog_url": "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml",
        "current_location": "online",
        "structured_mesh": False,
        "walltime": "01:00:00",   # UM data is from online (remote), needs more time
        "job_name": "um_glm_n2560_RAL3p3",
        "log_prefix": "extract_etc_UM",
        "job_groups": [
            {
                "label": "pr (hourly)",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8},
                "variables": ["pr"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "2D variables",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["tas", "huss", "ps", "psl", "uas", "vas", "prw"],
                "pressure_levels": None,
                "cof_mask": False,
            },
            {
                "label": "3D variables at 850 hPa (excl. wa)",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "850",
                "cof_mask": False,
            },
            {
                "label": "3D variables at 500 hPa (excl. wa)",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["ua", "va", "hus", "hur", "zg"],
                "pressure_levels": "500",
                "cof_mask": False,
            },
            {
                "label": "wa at 850 hPa (→ omega)",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["wa"],
                "pressure_levels": "850",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "wa at 500 hPa (→ omega)",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": ["wa"],
                "pressure_levels": "500",
                "convert_wa_to_omega": True,
                "cof_mask": False,
            },
            {
                "label": "COF masks",
                "catalog_model": "um_glm_n2560_RAL3p3",
                "catalog_params": {"zoom": 8, "time": "PT3H"},
                "variables": COF_MASK_VARS,
                "pressure_levels": None,
                "cof_mask": True,
            },
        ],
    },
}


# ============================================================
# SLURM SCRIPT GENERATION & SUBMISSION
# ============================================================

def build_python_cmd(model_cfg, group):
    """
    Build the python extract_etc_2d_vars.py command string for a given group.
    The VARIABLE placeholder will be filled in by the array task via $1.
    """
    # Resolve per-group catalog overrides
    catalog_url = group.get("catalog_url_override", model_cfg["catalog_url"])
    current_location = group.get("current_location_override", model_cfg["current_location"])

    catalog_params_str = json.dumps(group["catalog_params"])

    parts = [
        f'python {EXTRACT_SCRIPT}',
        f'  --catalog_url "{catalog_url}"',
        f'  --catalog_model "{group["catalog_model"]}"',
        f'  --catalog_params \'{catalog_params_str}\'',
        f'  --trackfile "{model_cfg["track_file"]}"',
        f'  --output_dir "{model_cfg["output_dir"]}"',
        f'  --variables "$CURRENT_VAR"',
        f'  --radius {RADIUS} --lon_res {LON_RES} --lat_res {LAT_RES}',
        f'  --chunk_size {CHUNK_SIZE} --progress_freq {PROGRESS_FREQ}',
    ]

    # current_location (only add if non-empty)
    if current_location:
        parts.append(f'  --current_location "{current_location}"')

    # Pressure levels
    if group.get("pressure_levels"):
        parts.append(f'  --pressure_levels {group["pressure_levels"]}')

    # Vertical velocity conversion
    if group.get("convert_wa_to_omega"):
        parts.append("  --convert_wa_to_omega")

    # COF mask flag
    if group.get("cof_mask"):
        parts.append("  --cof_mask")

    # Structured mesh (ERA5)
    if model_cfg.get("structured_mesh"):
        parts.append("  --structured_mesh")

    return " \\\n".join(parts)


def generate_sbatch_script(model_cfg, group, group_idx):
    """
    Generate the text of a self-contained sbatch script for one variable group.
    Each array task picks one variable from the group's variable list.
    """
    n_vars = len(group["variables"])
    array_spec = f"0-{n_vars - 1}%{min(n_vars, 10)}"

    # Build the VARIABLES bash array literal
    vars_line = " ".join(f'"{v}"' for v in group["variables"])

    # Build a concise description for the job name (max 15 chars for SLURM)
    plevel = group.get("pressure_levels") or "2D"
    cof_tag = "_cof" if group.get("cof_mask") else ""
    job_name = f'{model_cfg["job_name"]}_{plevel}{cof_tag}'[:15]

    log_name = f'{model_cfg["log_prefix"]}_g{group_idx}_%A_%a.log'
    python_cmd = build_python_cmd(model_cfg, group)

    script = f"""#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t {model_cfg["walltime"]}
#SBATCH -J {job_name}
#SBATCH -A {SLURM_ACCOUNT}
#SBATCH --array={array_spec}
#SBATCH --output=logs/{log_name}
#SBATCH --mail-user={MAIL_USER}
#SBATCH --mail-type=FAIL,END

# Auto-generated by submit_etc_extraction_jobs.py
# Source : {model_cfg.get("_source_key", "")}
# Group  : {group_idx} — {group.get("label", "")}
# Vars   : {", ".join(group["variables"])}
# Model  : {group["catalog_model"]}
# Params : {json.dumps(group["catalog_params"])}
# Plevel : {group.get("pressure_levels") or "N/A (2D)"}
# COF    : {group.get("cof_mask", False)}

source activate /global/common/software/m1867/python/hackathon

mkdir -p logs
mkdir -p "{model_cfg["output_dir"]}"

# Variable list for this group (each array task picks one)
VARIABLES=({vars_line})
CURRENT_VAR="${{VARIABLES[$SLURM_ARRAY_TASK_ID]}}"

echo "================================================"
echo "SLURM Array Job : $SLURM_ARRAY_JOB_ID"
echo "Array Task ID   : $SLURM_ARRAY_TASK_ID"
echo "Processing var  : $CURRENT_VAR"
echo "Group           : {group_idx} — {group.get("label", "")}"
echo "Catalog model   : {group["catalog_model"]}"
echo "================================================"

{python_cmd}

EXIT_CODE=$?
echo "================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: $CURRENT_VAR completed at $(date)"
else
    echo "FAILED:  $CURRENT_VAR exit code $EXIT_CODE at $(date)"
fi
echo "================================================"
exit $EXIT_CODE
"""
    return script


def submit_group(model_cfg, group, group_idx, dry_run=False):
    """Write sbatch script to a temp file and submit (or print in dry-run mode)."""
    script_text = generate_sbatch_script(model_cfg, group, group_idx)
    n_vars = len(group["variables"])
    label = group.get("label", f"group {group_idx}")

    print(f"\n{'='*60}")
    print(f"Group {group_idx}: {label}")
    print(f"  Variables ({n_vars}): {', '.join(group['variables'])}")
    print(f"  Catalog model : {group['catalog_model']}")
    print(f"  Params        : {json.dumps(group['catalog_params'])}")
    plevel = group.get("pressure_levels")
    if plevel:
        print(f"  Pressure level: {plevel} hPa")
    if group.get("cof_mask"):
        print("  COF mask      : yes")
    if group.get("convert_wa_to_omega"):
        print("  Convert       : wa → omega")

    if dry_run:
        print("\n--- SBATCH SCRIPT (dry run) ---")
        print(script_text)
        print("--- END SCRIPT ---")
        return None

    # Write to a named temp file and submit
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sh", prefix=f"sbatch_{model_cfg['job_name']}_g{group_idx}_",
        dir=SCRIPT_DIR, delete=False
    ) as f:
        f.write(script_text)
        tmp_path = f.name

    try:
        result = subprocess.run(
            ["sbatch", tmp_path],
            capture_output=True, text=True, check=True
        )
        job_id = result.stdout.strip()
        print(f"  Submitted: {job_id}  (script: {os.path.basename(tmp_path)})")
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"  ERROR submitting group {group_idx}: {e.stderr}", file=sys.stderr)
        return None
    finally:
        # Keep the script for reference (useful for resubmitting failed tasks)
        pass


def main():
    parser = argparse.ArgumentParser(
        description="Submit SLURM job arrays for ETC 2D variable extraction.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python submit_etc_extraction_jobs.py --source ERA5
  python submit_etc_extraction_jobs.py --source SCREAM --dry-run
  python submit_etc_extraction_jobs.py --source NICAM --group 0
  python submit_etc_extraction_jobs.py --source ERA5 --list-groups
        """
    )
    parser.add_argument(
        "--source", required=True,
        choices=list(MODEL_CONFIGS.keys()),
        help="Source model to process"
    )
    parser.add_argument(
        "--group", type=int, default=None,
        help="Submit only this group index (0-based). Default: submit all groups."
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print sbatch scripts without submitting."
    )
    parser.add_argument(
        "--list-groups", action="store_true",
        help="List variable groups for the source and exit."
    )
    args = parser.parse_args()

    model_cfg = MODEL_CONFIGS[args.source].copy()
    model_cfg["_source_key"] = args.source
    job_groups = model_cfg["job_groups"]

    if args.list_groups:
        print(f"\nVariable groups for {args.source}:")
        for i, g in enumerate(job_groups):
            plevel = g.get("pressure_levels") or "2D"
            cof = " [COF]" if g.get("cof_mask") else ""
            print(f"  [{i}] {g.get('label', '')}  |  {len(g['variables'])} vars  |  {plevel}{cof}")
            print(f"       catalog_model = {g['catalog_model']}")
            print(f"       vars: {', '.join(g['variables'])}")
        return

    os.makedirs(SCRIPT_DIR / "logs", exist_ok=True)

    groups_to_submit = (
        [job_groups[args.group]] if args.group is not None else job_groups
    )
    group_indices = (
        [args.group] if args.group is not None else list(range(len(job_groups)))
    )

    mode = "DRY RUN" if args.dry_run else "SUBMITTING"
    print(f"\n{mode}: {len(groups_to_submit)} job group(s) for {args.source}")
    print(f"Track file : {model_cfg['track_file']}")
    print(f"Output dir : {model_cfg['output_dir']}")
    print(f"Walltime   : {model_cfg['walltime']}")

    submitted = []
    for gidx, grp in zip(group_indices, groups_to_submit):
        job_id = submit_group(model_cfg, grp, gidx, dry_run=args.dry_run)
        if job_id:
            submitted.append((gidx, grp.get("label", ""), job_id))

    if submitted and not args.dry_run:
        print(f"\n{'='*60}")
        print(f"Submitted {len(submitted)} job(s) for {args.source}:")
        for gidx, label, jid in submitted:
            print(f"  Group {gidx} ({label}): {jid}")


if __name__ == "__main__":
    main()
