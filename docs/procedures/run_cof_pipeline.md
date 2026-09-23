# Running Analyses 1 and 2 in One Go (`run_cof_pipeline.py`)

**Reference scripts:**
- `scripts/run_cof_pipeline.py` — runs every step of Analyses 1 and 2 for any set of sources, in dependency order, on one node
- `config/config_pipeline.yaml` — the sources, commands, CPU slots and estimated durations
- `scripts/check_zarr_store.py` — chunk-by-chunk completeness check of Zarr stores (used by the preflight, also usable alone)
- `slurm/run_interactive_cof_pipeline.sh`, `slurm/slurm_run_cof_pipeline.sh`, `slurm/slurm_run_cof_pipeline_array.sh` — launchers

**Author:** Zhe Feng | zhe.feng@pnnl.gov

---

## Why one runner

Running each step as its own Slurm array means one queue wait per step and a manual hand-over between steps. The runner starts
every step as soon as the steps it needs have finished and runs the sources side by side on one node, so a full re-process needs one
allocation (none to wait for on the interactive QOS) and about 1.6 h (95 min measured for all six sources). The time and memory of that run are in "Resources and time" below.

## The dependency graph

For each source:

```
s1 (make_mcs_swath_masks) -> s2 (combine_tracking_masks) -> s3 (make_cooccurrence_masks) -> monthly            [Analysis 1]
s1 -> thresholds (calc_extreme_precip_thresholds)  \
                                                     +-> attribution (calc_stormtype_extreme_precip_spatial)  [Analysis 2]
s3 -------------------------------------------------/
```

| Step | Script | Output (under the data root) | Needs |
|------|--------|------------------------------|-------|
| `s1` | `make_mcs_swath_masks.py` | `mcs_masks/{source}_mcs_masks_hp8.zarr` | the hourly MCS masks, catalog `pr` and Tb |
| `s2` | `combine_tracking_masks.py` (IMERG: `combine_era5_imerg_tracking_masks.py`) | `all_masks/{source}_allmasks_hp8_v1.zarr` | `s1` |
| `s3` | `make_cooccurrence_masks.py` | `cof_masks/{source}_cofmasks_hp8_v1.zarr` | `s2` |
| `monthly` | `calc_monthly_rainmap_by_cof.py` | `cof_masks/stats/monthly/{source}_monthly_rainmap_cof_hp8_v1.nc` | `s3` |
| `thresholds` | `calc_extreme_precip_thresholds.py` | `extreme_precip/{source}_precip_percentiles_6h_hp8_v1.nc` | `s1` for the models (they read Step 1's `tot_pr`); nothing for IMERG (the non-IR 6-hourly store) |
| `attribution` | `calc_stormtype_extreme_precip_spatial.py` | `extreme_precip/{source}_stormtype_spatial_p90.nc`, `..._p95.nc` | `s3` and `thresholds` |

A step starts only when every step it needs has finished with exit status 0. A failed step skips the later steps of its own source
only; the other sources go on. Steps that do not depend on each other (the thresholds and Steps 2 and 3, several sources at once) run in
parallel, as many as fit in the CPU-slot budget.

## Running it

**Interactive node** (no queue wait, 4 h limit), the recommended way for a full re-process:

```bash
salloc -N 1 -C cpu -q interactive -t 04:00:00 -A m1867      # then log in to the node if salloc does not
cd /global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature
bash slurm/run_interactive_cof_pipeline.sh --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/round2 --dry-run   # graph, commands, estimate
bash slurm/run_interactive_cof_pipeline.sh --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/round2             # the run, in the background
```

The launcher checks the options and the inputs first, then starts the runner with `nohup`, so a lost connection does not stop it, and
prints how to follow it. `--foreground` keeps it in the terminal.

**One Slurm job on one node** (all sources, one queue wait): `sbatch --export=ALL,DATA_ROOT=DIR slurm/slurm_run_cof_pipeline.sh`
(3 h; `SOURCES="scream icon"` and `PIPELINE_ARGS="--analysis 2"` are optional). The job forwards the time-limit warning to the runner,
which stops its steps cleanly.

**One source per node** (fallback if one node is not enough): `sbatch --export=ALL,DATA_ROOT=DIR slurm/slurm_run_cof_pipeline_array.sh`
(`--array=1-6`, index order IMERGv7, scream, icon_d3hp003, nicam_gl11, um_glm_n2560_RAL3p3, casesm2_10km_nocumulus).

**Directly:** `python scripts/run_cof_pipeline.py --data-root DIR [options]` with the `hackathon` environment.

### Choosing what runs

| Option | Effect |
|--------|--------|
| `--sources scream icon` | only these sources (names or aliases such as `imerg`, `um`, `IR_IMERG`); default all |
| `--analysis 1` / `2` / `both` | Analysis 1 = s1 s2 s3 monthly; Analysis 2 = s1 s2 s3 thresholds attribution; default both |
| `--steps s3 monthly` | only these steps; the steps they need must already exist (or be selected) |
| `--from s3` | this step and the ones after it in the order s1 s2 s3 monthly thresholds attribution |
| `--resume` | skip steps whose marker says success (and whose upstream steps did not run again); re-run the rest |
| `--force` | overwrite outputs that the runner did not create |
| `--dry-run` | print the graph, the commands and an estimated schedule; run nothing |
| `--preflight-only` | run the input checks and stop |
| `--max-slots N`, `--stagger-sec S`, `--min-free-gb G` | CPU-slot budget (default 80% of the logical CPUs), seconds between step starts (20), memory that must be available to start a step (60 GB) |
| `--step-args STEP "ARGS"` | extra arguments for one step, for tests or a special run |

Step aliases (either form works everywhere a step name does): `1`/`step1` → `s1`, `2`/`step2` → `s2`, `3`/`step3` → `s3`, `thr`/`extreme` → `thresholds`, `attr` → `attribution`.

### Example: thresholds and attribution only, against a different threshold period

`--step-args` appends to the command the runner already builds for that step, and a later occurrence of an option on the
command line overrides an earlier one (plain argparse behavior) — so pointing `thresholds` at a different Zarr store needs
no `config_sources.yaml` edit, even for IMERG, whose registry entry (`config/config_pipeline.yaml:54-56`) already passes its
own `--input_zarr`/`--input_var`:

```bash
# Side data-root: symlink cof_masks/ to the real one so Step 3 is reused, not re-run. Thresholds+attribution land under
# the side root only; production is untouched.
mkdir -p /pscratch/sd/w/wcmca1/hackathon/tmp/imerg20yr
ln -s /pscratch/sd/w/wcmca1/hackathon/cof_masks /pscratch/sd/w/wcmca1/hackathon/tmp/imerg20yr/cof_masks

python scripts/run_cof_pipeline.py \
  --data-root /pscratch/sd/w/wcmca1/hackathon/tmp/imerg20yr \
  --sources imerg --steps thresholds attribution \
  --step-args thresholds "--input_zarr /path/to/IMERG_20yr_store.zarr --input_var precipitation" \
  --dry-run          # drop --dry-run to run; expect one WARNING that s3's output has no runner marker (harmless: it's
                      # production's own Step 3 output, reused read-only through the symlink, not something this run made)
```

```bash
# In place: overwrites the current threshold+attribution files under production (needs --force: those files predate this
# run and have no runner marker).
python scripts/run_cof_pipeline.py --data-root /pscratch/sd/w/wcmca1/hackathon/ \
  --sources imerg --steps thresholds attribution --force \
  --step-args thresholds "--input_zarr /path/to/IMERG_20yr_store.zarr --input_var precipitation"
```

Both commands work for any source, not just IMERG — swap `--sources imerg` and the `--input_zarr` value.

**Trap:** `calc_stormtype_extreme_precip_spatial.py` always reads `extreme_precip/{source_name}_precip_percentiles_6h_hp8_v1.nc`
— there is no `--version`/`--threshold_file` flag on the attribution step to point it at a differently-versioned thresholds
file. A `thresholds` run with `--step-args thresholds "--version v20yr ..."` therefore produces a file the attribution step
will never read; leave `--version` at its default (as both commands above do) so attribution finds it. If you want the
20-year thresholds purely as a separate, differently-named artifact instead — not meant to feed the attribution at all —
call `scripts/calc_extreme_precip_thresholds.py` directly with `--version` and skip the runner for that step.

A source with several years of qualifying data will also get the per-calendar-year percentiles and their interannual IQR
(new `year` coordinate, `pr_annual_p*`/`pr_q25_p*`/`pr_q75_p*`/`pr_iqr_p*` variables) once at least 2 calendar years clear
`--min_year_coverage_days` (default 300 distinct days); add `--no_annual` to the `--step-args` string above to skip that.

## Where things go, and what is protected

`--data-root` is required. Every script reads `COF_DATA_ROOT` (see `src/cof_paths.py`, default the production tree), and the runner sets
it for each step from `--data-root`, so a test area gets its own `mcs_masks/`, `all_masks/`, `cof_masks/` and `extreme_precip/`
and production is never touched. The hourly MCS masks, the tracking files, the catalog data and the ERA5 masks are inputs and are not under the root.

- The production tree holds the round-2 outputs since 2026-09-20 (moved from `tmp/round2`, whose markers and logs stayed there). They have no markers under the production root, so a run with
  `--data-root /pscratch/sd/w/wcmca1/hackathon/` stops at once and lists them unless `--force` is given: that is the protection working as intended. The previous production is in `_prev_production_20260918/`.
- A **marker** `<data-root>/pipeline_state/<source>/<step>.json` records each step (`running`, `success`, `failed`, `interrupted`) with the command, times, host,
  code version and peak memory. It is written only by the runner.
- An output that exists **without a marker** is never overwritten unless `--force` is given. That is what protects a production tree: pointing
  `--data-root` at it stops at once, listing the files. Outputs the runner made itself can be re-made (a step with a `failed`, `interrupted` or
  `success` marker).
- `--resume` continues a stopped or timed-out run: steps with a success marker are skipped, but a step whose upstream step runs again
  is re-run too, so a later step never sits on a stale input.
- If a step needs an earlier step that is neither selected nor present, the runner says so and stops (an existing output without a marker is accepted, with a warning).

Logs go to `<data-root>/pipeline_logs/<run id>/`: one log per step, `status.json` (running, done, failed, pending; updated every few seconds), `resources.csv`
(node memory available, resident memory per step every 30 s), `summary.txt` and, from the launcher, `launch_<time>.out` with the runner's events.

## Preflight

Before anything starts, the runner checks (a few seconds, directory listings only): that the environment has zarr 2 (zarr 3 writes another layout); that the
hourly MCS mask store of each source, the IMERG IR store (`IR_IMERG_V7_1H_zoom8_...zarr`, the one Step 1 reads; the config's `zoom: 9` is overridden) and the
ERA5 mask store are complete, chunk file by chunk file, with `scripts/check_zarr_store.py`; and that the AR, TC and ETC files of each model exist
(they live in another account's scratch). Zarr returns NaN for a missing chunk without any error, so an incomplete input store would otherwise
give silently wrong results. `--skip-preflight` turns it off.

## Resources and time

`config/config_pipeline.yaml` gives each step a number of CPU slots (workers x threads) and an estimated duration; `--dry-run` schedules them. The runner starts the source
with the longest remaining chain first (IMERG).

**Measured in the full run of 2026-09-20** (all six sources, one 503 GB / 256-CPU node, 204 slots, the runner's defaults; 36 steps, all exit status 0): **94.6 min** wall time, at most 192 of 204 slots
busy. The six Step 1 runs overlapped for the first half hour; the summed resident memory of all steps peaked at **266 GB** (5.5 min after the start) and the node never had less than **221 GB available**.
IMERG's Step 1 is the critical path; after about 35 min only IMERG's steps run. Minutes per step, and the peak resident memory of the step in GB (sampled every 30 s, so short peaks can be missed):

| Source | s1 | s2 | s3 | monthly | thresholds | attribution |
|--------|----|----|----|---------|------------|-------------|
| IMERGv7 (48 workers in Step 1) | 55 (87) | 0.6 (2) | 21 (21) | 3.0 (94) | 3.0 (4) | 18 (31) |
| SCREAM (32 workers) | 35 (78) | 0.8 (27) | 9.7 (16) | 1.9 (39) | 1.7 (2) | 9.7 (33) |
| ICON | 34 (64) | 0.8 (7) | 9.2 (16) | 1.8 (36) | 1.5 (2) | 8.2 (32) |
| NICAM | 32 (59) | 0.8 (26) | 9.3 (15) | 1.7 (35) | 1.8 (2) | 8.8 (31) |
| UM | 34 (55) | 0.7 (22) | 9.8 (15) | 1.6 (41) | 1.6 (2) | 8.2 (31) |
| CASESM2 | 29 (44) | 1.1 (16) | 9.6 (15) | 1.3 (34) | 2.0 (2) | 8.9 (33) |

Running the sources together costs time per step: with five chains at once (the earlier test run) the models' Step 1 took 23-27 min, with all six 29-35 min and IMERG's 55 min instead of 46. The memory of Step 1 is larger
than a first snapshot suggested (about 2 GB per worker in steady state, not 0.6 GB). Each Dask cluster sizes its worker memory limit from the total node memory (80% / workers), so several
clusters on one node each assume they own it; the runner therefore holds back new steps below `--min-free-gb` and starts steps `--stagger-sec` apart, and every run records the memory of every step in `resources.csv`.
To shorten a run, give IMERG's Step 1 more workers (its `slots` in the registry): after the models finish the node is mostly idle.

**The results do not depend on how the steps are scheduled.** The five model sources from this run are bit-identical to the same sources run as manual chains at every stage (Steps 1-3, all variables; the monthly file;
the thresholds; the attribution). See `docs/AUDIT_FINDINGS.md` for the comparison with IMERG, which differs only through the removal of rain at pixels with missing Tb.

## Interrupting and troubleshooting

- `kill -TERM <pid>` or Ctrl-C stops the runner: it terminates the running steps (their whole process group), marks them `interrupted` and exits with status 1. Continue with `--resume`.
- Exit status: 0 all steps done; 1 a step failed, was skipped because of a failure, or the run was interrupted; 2 usage error, refused overwrite, missing prerequisite or failed preflight.
- A step that failed: read its log under `pipeline_logs/<run id>/<source>_<step>.log`, fix the cause, run again with `--resume`.
- Tiny test stores: a test with `--step-args s1 "--test-steps N"` gives a store of N/6 windows, fewer than the time chunk of the Step 2 writers (24 frames; 28 in the IMERG Step 2 script). Those writers used to
  divide by zero in that case; they now write one chunk. Checked for SCREAM with N = 24 through Steps 1-3; the later steps (monthly, thresholds, attribution) were run on tiny stores with N = 240 (40 windows).
- To check a store by hand: `python scripts/check_zarr_store.py STORE [STORE ...]` (exit status 1 when a chunk is missing or empty, when the store has no array metadata at all,
  as after a scratch purge that leaves empty directories, or when a requested `--arrays` name is not in it). Tests: `python tests/test_check_zarr_store.py`.

## Adding a source

Add an entry to `config/config_pipeline.yaml` (its `config_key` in `config_sources.yaml`, the Step 1 config, the Step 2 variant, where its thresholds come from, optional slots and
estimates), and the source name to the `--array` order of `slurm/slurm_run_cof_pipeline_array.sh`. No code change is needed.
