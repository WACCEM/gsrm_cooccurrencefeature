"""Tests of the label and count changes of the ETC and MCS-COF notebooks (notebooks/*.ipynb).

  - labels: COF names joined by "+" (ETC+MCS, ETC + MCS + AR, MCS+AR+ETC, the table labels A+M and E+A) are now joined by "-" in the six notebooks
    plot_etc_composites, _diff, _multipanel, plot_etc_spatialmean_stats_1source, _multisource and plot_mcs_cof_trackstats_multisource
  - the six notebooks are valid nbformat and every code cell parses (the missing line break "/ 100track_stats.head()" of the 1-source notebook is fixed, so it runs top to bottom)
  - counts: the three notebooks with sample sizes in their figures (1source, multisource, mcs_cof_trackstats_multisource) define the same helpers
    estimate_n_years and format_per_year (years = span of the time stamps / 365.25, robust to a record that crosses a calendar year, to gaps and to NaT;
    an average count per year is formatted by format_per_year)
  - the plot functions of the notebooks show the average count per year in their legend and say so in the legend title (run on synthetic data)

Run:  python tests/test_notebook_labels_counts.py      (or pytest tests/)
"""
import ast
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
NB_DIR = REPO / "notebooks"
LABEL_NOTEBOOKS = ["plot_etc_composites_diff", "plot_etc_composites", "plot_etc_composites_multipanel", "plot_etc_spatialmean_stats_1source",
                   "plot_etc_spatialmean_stats_multisource", "plot_mcs_cof_trackstats_multisource"]
COUNT_NOTEBOOKS = ["plot_etc_spatialmean_stats_1source", "plot_etc_spatialmean_stats_multisource", "plot_mcs_cof_trackstats_multisource"]
PLUS_LABEL = re.compile(r"\b(?:ETC|MCS|AR)\s?\+\s?(?:ETC|MCS|AR)")
HELPER_MARKER = "# ── Years covered by each source"


def load(name):
    return json.loads((NB_DIR / f"{name}.ipynb").read_text(encoding="utf-8"))


def sources(nb):
    return ["".join(c["source"]) for c in nb["cells"]]


def helper_text(name):
    """The helper definitions (before the notebook-specific years code) of the cell that defines estimate_n_years."""
    cells = [s for s in sources(load(name)) if "def estimate_n_years" in s]
    assert len(cells) == 1, (name, len(cells))
    return cells[0].split(HELPER_MARKER)[0]


def test_no_plus_joined_cof_labels():
    for name in LABEL_NOTEBOOKS:
        for i, s in enumerate(sources(load(name))):
            assert not PLUS_LABEL.search(s), f"{name} cell {i}: {PLUS_LABEL.search(s).group(0)!r}"
    assert not any(re.search(r"\bA\+M\b", s) for s in sources(load("plot_etc_spatialmean_stats_multisource")))
    assert not any(re.search(r"\bE\+A\b", s) for s in sources(load("plot_mcs_cof_trackstats_multisource")))
    # the new labels are there
    multi = "\n".join(sources(load("plot_etc_spatialmean_stats_multisource")))
    assert "'ETC-AR-MCS'" in multi and "('ETC-MCS', 'M')" in multi and "'A-M'" in multi
    a4 = "\n".join(sources(load("plot_mcs_cof_trackstats_multisource")))
    assert "'MCS-AR-ETC'" in a4 and "('MCS-ETC', 'E')" in a4 and "'E-A'" in a4
    assert "'ETC-AR_ar'" in "\n".join(sources(load("plot_etc_spatialmean_stats_1source")))        # the internal keys were renamed together with their lookups


def test_notebooks_are_valid_and_parse():
    import nbformat
    for name in LABEL_NOTEBOOKS:
        nbformat.validate(nbformat.read(str(NB_DIR / f"{name}.ipynb"), as_version=4))
        nb = load(name)
        for i, c in enumerate(nb["cells"]):
            if c["cell_type"] != "code":
                continue
            src = "".join(c["source"])
            code = "\n".join(l for l in src.split("\n") if not l.lstrip().startswith(("%", "!")))
            try:
                ast.parse(code)
            except SyntaxError as e:
                raise AssertionError(f"{name} cell {i}: {e}")
        ids = [c["id"] for c in nb["cells"]]
        assert len(ids) == len(set(ids)), f"{name}: duplicate cell ids"
        errors = [i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code" for o in c.get("outputs", []) if o.get("output_type") == "error"]
        assert not errors, f"{name}: saved error output in cells {errors}"                           # the saved outputs are those of a run without errors


def test_helpers_are_identical_and_correct():
    texts = {name: helper_text(name) for name in COUNT_NOTEBOOKS}
    assert len(set(texts.values())) == 1, "the three notebooks must define identical helpers"
    ns = {"pd": pd, "np": np}
    exec(texts[COUNT_NOTEBOOKS[0]], ns)
    yrs, fmt = ns["estimate_n_years"], ns["format_per_year"]
    scream = pd.date_range("2019-08-01", "2020-08-31 18:00", freq="6h")             # crosses a calendar year, 13 months
    era5 = pd.date_range("2019-01-01", "2021-12-31 18:00", freq="6h")
    nicam = pd.date_range("2020-03-01 06:00", "2021-02-28 18:00", freq="6h")
    assert 1.08 < yrs(scream) < 1.10 and len(set(scream.year)) == 2                   # 1.09 years, not 2 calendar years
    assert 2.99 < yrs(era5) < 3.01
    assert 0.99 < yrs(nicam) <= 1.0 and len(set(nicam.year)) == 2
    assert 0.99 < yrs(pd.date_range("2020-07-01", "2021-06-30", freq="1D")) < 1.01
    gap = pd.DatetimeIndex(list(pd.date_range("2019-01-01", "2019-03-31", freq="1D")) + list(pd.date_range("2020-10-01", "2020-12-31", freq="1D")))
    assert 1.98 < yrs(gap) < 2.0                                                      # the span, not the number of months present
    assert abs(yrs(pd.Series(scream)) - yrs(scream)) < 1e-12 and abs(yrs(np.asarray(scream.values)) - yrs(scream)) < 1e-12
    assert abs(yrs(pd.Series([pd.NaT, pd.Timestamp("2019-01-01"), pd.Timestamp("2020-01-01"), pd.NaT])) - 365 / 365.25) < 1e-12      # NaT ignored
    for bad in ([], [pd.NaT], [pd.Timestamp("2020-01-01")], [pd.Timestamp("2020-01-01")] * 3):
        try:
            yrs(pd.Series(bad, dtype="datetime64[ns]")); raise AssertionError("must fail")
        except ValueError:
            pass
    assert fmt(630, 3.0) == "210" and fmt(5, 1.16) == "4.3" and fmt(0, 2.0) == "0" and fmt(29.9, 3.0) == "10" and fmt(1234, None) == "1,234" and fmt(3000, 1.0) == "3,000"
    assert fmt(6, 1.09) == "5.5" and fmt(124, 1.09) == "114"


def test_boxplot_legends_show_counts_per_year():
    ns = {"pd": pd, "np": np, "plt": plt, "sns": sns, "List": List, "Dict": Dict, "Tuple": Tuple, "Optional": Optional}
    exec(helper_text("plot_etc_spatialmean_stats_multisource"), ns)
    # --- ETC multisource: plot_boxplots_by_etc_type
    nb = load("plot_etc_spatialmean_stats_multisource")
    fn_cell = [s for s in sources(nb) if "def plot_boxplots_by_etc_type" in s][0]
    exec(fn_cell, ns)
    rows = []
    for src, counts in (("ERA5", (630, 598, 47, 435)), ("SCREAM", (124, 296, 6, 153))):
        for typ, n in zip(["Isolated", "ETC-MCS", "ETC-AR", "ETC-AR-MCS"], counts):
            rows += [{"source_display": src, "overlap_type": typ, "v": float(k % 7)} for k in range(n)]
    df = pd.DataFrame(rows)
    kw = dict(figsize=(6, 4), source_order=["ERA5", "SCREAM"], source_colors={"ERA5": "gray", "SCREAM": "blue"})
    fig = ns["plot_boxplots_by_etc_type"](df, [["v"]], {"v": "v"}, n_years={"ERA5": 3.0, "SCREAM": 1.09}, **kw)
    leg = fig.legends[0]
    assert [t.get_text() for t in leg.get_texts()] == ["ERA5 (210 | 199 | 16 | 145)", "SCREAM (114 | 272 | 5.5 | 140)"], [t.get_text() for t in leg.get_texts()]
    assert leg.get_title().get_text() == "Average tracks per year (Isolated | ETC-MCS | ETC-AR | ETC-AR-MCS)"
    plt.close(fig)
    fig = ns["plot_boxplots_by_etc_type"](df, [["v"]], {"v": "v"}, **kw)                   # without n_years: the total counts as before, no title
    assert [t.get_text() for t in fig.legends[0].get_texts()] == ["ERA5 (630 | 598 | 47 | 435)", "SCREAM (124 | 296 | 6 | 153)"]
    assert fig.legends[0].get_title().get_text() == ""
    plt.close(fig)
    # --- Analysis 4: plot_boxplots_by_mcs_type
    ns = {"pd": pd, "np": np, "plt": plt, "sns": sns, "List": List, "Dict": Dict, "Tuple": Tuple, "Optional": Optional}
    exec(helper_text("plot_mcs_cof_trackstats_multisource"), ns)
    exec([s for s in sources(load("plot_mcs_cof_trackstats_multisource")) if "def plot_boxplots_by_mcs_type" in s][0], ns)
    rows = []
    for ds, counts in (("OBS", (900, 600, 300, 150)), ("SCREAM", (100, 80, 40, 10))):
        for typ, n in zip(["Isolated", "MCS-AR", "MCS-ETC", "MCS-AR-ETC"], counts):
            rows += [{"dataset": ds, "cof_type": typ, "v": float(k % 5)} for k in range(n)]
    df = pd.DataFrame(rows)
    kw = dict(figsize=(6, 4), dataset_order=["OBS", "SCREAM"], source_colors={"OBS": "gray", "SCREAM": "blue"})
    fig = ns["plot_boxplots_by_mcs_type"](df, [["v"]], {"v": "v"}, n_years={"OBS": 3.0, "SCREAM": 1.09}, **kw)
    assert [t.get_text() for t in fig.legends[0].get_texts()] == ["OBS (300 | 200 | 100 | 50)", "SCREAM (92 | 73 | 37 | 9.2)"], [t.get_text() for t in fig.legends[0].get_texts()]
    assert fig.legends[0].get_title().get_text() == "Average tracks per year (Isolated | MCS-AR | MCS-ETC | MCS-AR-ETC)"
    plt.close(fig)


def test_calls_pass_the_years_and_lifecycle_code_uses_them():
    multi = sources(load("plot_etc_spatialmean_stats_multisource"))
    a4 = sources(load("plot_mcs_cof_trackstats_multisource"))
    assert sum(s.count("n_years=source_n_years,") for s in multi) == 3                       # two lifecycle calls and the box plot
    life = [s for s in multi if "def plot_lifecycle_evolution_by_etc_type" in s][0]
    assert "n_years: Optional[Dict[str, float]] = None," in life and "format_per_year(_type_counts_pivot.loc[lbl, t], (n_years or {}).get(lbl))" in life
    assert "title=legend_title," in life and "Average tracks per year" in life
    for fn in [c for c in multi if "def plot_boxplots_by_etc_type" in c] + [life] + [c for c in a4 if "def plot_boxplots_by_mcs_type" in c]:
        assert "title_shift = 1.5 * fontsize / (fig.get_figheight() * 72) if legend_title else 0.0" in fn      # the legend title does not move the layout
        assert "bbox_to_anchor=(0.5, legend_y + title_shift)," in fn and "y=title_y + title_shift)" in fn
    assert sum(s.count("n_years=n_years_by_dataset,") for s in a4) == 3                      # all, ocean, land
    one = "\n".join(sources(load("plot_etc_spatialmean_stats_1source")))
    assert one.count("format_per_year(") >= 9 and "n_years = estimate_n_years(df.loc[df['overlap_flag'].notna(), 'time'])" in one
    assert "f'n={format_per_year(count, n_years)}/yr'" in one and "tracks/yr" in one


if __name__ == "__main__":
    test_no_plus_joined_cof_labels(); test_notebooks_are_valid_and_parse(); test_helpers_are_identical_and_correct(); test_boxplot_legends_show_counts_per_year()
    test_calls_pass_the_years_and_lifecycle_code_uses_them()
    print("test_notebook_labels_counts: all checks passed")
