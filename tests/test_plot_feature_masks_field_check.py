"""Tests of resolve_field_layers in scripts/plot_feature_masks_with_field.py (2026-09-24):

  - the requested shading and contour variables are checked against the variables of the field dataset; a missing one is
    skipped with a warning instead of failing every frame, and a layer that was not requested stays off without a warning
  - IVT counts as present when uivt and vivt are (SCREAM has them, ICON has neither)

Run:  python tests/test_plot_feature_masks_field_check.py      (or pytest tests/)
"""
import contextlib
import io
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import plot_feature_masks_with_field as pf  # noqa: E402

SCREAM_VARS = ["psl", "uivt", "vivt", "ps"]
ICON_VARS = ["psl", "ua", "va", "hus", "zg", "orog", "rlut"]


def check(*args):
    """(result, printed text) of resolve_field_layers."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        result = pf.resolve_field_layers(*args)
    return result, buf.getvalue()


def test_layers_are_kept_when_their_variables_exist():
    assert check(["IVT", "psl"], "IVT", "psl", True, True) == ((True, True), "")
    assert check(SCREAM_VARS, "IVT", "psl", True, True) == ((True, True), "")     # IVT from uivt and vivt
    assert check(ICON_VARS, "rlut", "psl", True, True) == ((True, True), "")      # any other shading variable by name


def test_missing_shading_variable_is_skipped_with_a_warning():
    result, out = check(ICON_VARS, "IVT", "psl", True, True)
    assert result == (False, True)                                              # the contours stay
    assert "'IVT'" in out and "skipping the background shading" in out and "contour" not in out
    assert check(["psl", "uivt"], "IVT", "psl", True, True)[0] == (False, True)  # one of uivt, vivt is not enough
    assert check(SCREAM_VARS, "rlut", "psl", True, True)[0] == (False, True)     # a named variable must exist itself


def test_missing_contour_variable_is_skipped_with_a_warning():
    result, out = check(["IVT"], "IVT", "psl", True, True)
    assert result == (True, False)                                              # the shading stays
    assert "'psl'" in out and "skipping the contours" in out and "shading" not in out


def test_both_missing_gives_two_warnings_and_no_layers():
    result, out = check([], "IVT", "psl", True, True)
    assert result == (False, False)
    assert out.count("Warning") == 2


def test_layers_that_were_not_requested_stay_off_and_do_not_warn():
    assert check(ICON_VARS, "IVT", "psl", False, False) == ((False, False), "")
    assert check([], "IVT", "psl", False, False) == ((False, False), "")
    assert check([], "IVT", "psl", False, True)[0] == (False, False)             # only the requested layer is checked


if __name__ == "__main__":
    test_layers_are_kept_when_their_variables_exist()
    test_missing_shading_variable_is_skipped_with_a_warning()
    test_missing_contour_variable_is_skipped_with_a_warning()
    test_both_missing_gives_two_warnings_and_no_layers()
    test_layers_that_were_not_requested_stay_off_and_do_not_warn()
    print("test_plot_feature_masks_field_check: all checks passed")
