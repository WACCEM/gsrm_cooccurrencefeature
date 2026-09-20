"""
Where the COF pipeline reads and writes its data products.

Every script of the pipeline (Steps 1-3, the monthly rain map, the extreme-precipitation thresholds and attribution) keeps its
products under one root: mcs_masks/, all_masks/, cof_masks/, cof_masks/stats/monthly/ and extreme_precip/. The default is the
production tree. Set the environment variable COF_DATA_ROOT (scripts/run_cof_pipeline.py does this from --data-root) to work in
another tree, for example a test area, without touching production. Inputs that do not come from the pipeline itself (the hourly
MCS masks, the tracking files, the catalog data) are not under this root.

Author: Zhe Feng | zhe.feng@pnnl.gov
"""

import os

DEFAULT_DATA_ROOT = "/pscratch/sd/w/wcmca1/hackathon/"
ENV_VAR = "COF_DATA_ROOT"


def data_root(logger=None):
    """
    Root directory of the pipeline products, with a trailing slash.

    Parameters:
    -----------
    logger : logging.Logger, optional
        When given, logs the root once per call so that every run states where it reads and writes. A root other than the
        default is logged as a warning, because it comes from the environment and not from the command line.
    """
    root = os.environ.get(ENV_VAR) or DEFAULT_DATA_ROOT
    root = root.rstrip("/") + "/"
    if logger is not None:
        if root == DEFAULT_DATA_ROOT:
            logger.info(f"Data root: {root} (default)")
        else:
            logger.warning(f"Data root: {root} (from {ENV_VAR}, not the production default {DEFAULT_DATA_ROOT})")
    return root
