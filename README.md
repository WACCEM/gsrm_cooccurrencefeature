# Co-occurring Feature Analysis for Global Storm-Resolving Models

This repository identifies and analyzes co-occurring atmospheric features in global storm-resolving models (GSRMs), with emphasis on precipitation attribution, extreme precipitation, ETC-centered composites, and MCS lifecycle statistics.

**Author:** Zhe Feng (zhe.feng@pnnl.gov)

## Features Analyzed

- Tropical cyclones (TC)
- Extratropical cyclones (ETC)
- Atmospheric rivers (AR)
- Mesoscale convective systems (MCS)
- Non-MCS cloud type and precipitation categories

## Environment

On NERSC Perlmutter:

```bash
source activate /global/common/software/m1867/python/hackathon
cd /path/to/gsrm_cooccurrencefeature
```

## Pipeline Overview

The full workflow is described in the [analysis pipeline overview](docs/pipelines/overview.md). The main analysis paths are:

1. **Shared COF mask pipeline:** build MCS swath masks, combine tracking masks, and identify co-occurrence feature masks.
2. **Total precipitation by COF type:** compute monthly precipitation maps attributed to isolated and overlapping feature categories.
3. **Extreme precipitation by COF type:** compute precipitation thresholds and attribute extreme events to storm and COF categories.
4. **ETC composite analysis:** extract ETC-centered 2D fields and compute spatial statistics by COF overlap type.
5. **MCS COF track statistics:** merge MCS track statistics with COF overlap flags across sources.

## Where to Start

| Goal | Start Here | Main Scripts / Notebooks |
|------|------------|--------------------------|
| Understand the full workflow | [Analysis pipeline overview](docs/pipelines/overview.md) | Multiple processing scripts |
| Find docs for a specific script | [Documentation index](docs/index.md) | All documented scripts and notebooks |
| Build MCS swath masks | [MCS swath procedure](docs/procedures/mcs_swath_cloud_type.md) | `scripts/make_mcs_swath_masks.py` |
| Combine feature masks | [Combined mask procedure](docs/procedures/combine_tracking_masks.md) | `scripts/combine_tracking_masks.py` |
| Identify COF masks | [COF identification procedure](docs/procedures/cof_identification.md) | `scripts/make_cooccurrence_masks.py` |
| Compute monthly precipitation by COF | [Monthly precipitation procedure](docs/procedures/monthly_precip_by_cof.md) | `scripts/calc_monthly_rainmap_by_cof.py` |
| Analyze extreme precipitation | [Extreme precipitation procedure](docs/procedures/extreme_precip_by_stormtype.md) | `scripts/calc_extreme_precip_thresholds.py`, `scripts/calc_stormtype_extreme_precip_spatial.py` |
| Run ETC environment extraction | [Batched extraction notes](extract_environments/README_BATCHED_EXTRACTION.md) | `extract_environments/extract_etc_2d_vars.py` |

## Repository Layout

| Path | Contents |
|------|----------|
| `config/` | Source and processing configuration files |
| `docs/` | Documentation index, pipeline overview, procedure docs, and assets |
| `extract_environments/` | ETC-centered environmental variable extraction workflow |
| `notebooks/` | Analysis, diagnostic, and plotting notebooks |
| `scripts/` | Processing, aggregation, plotting, and batch helper scripts |
| `slurm/` | NERSC Slurm job scripts and task lists |
| `src/` | Shared Python utilities |

See [docs/index.md](docs/index.md) for the full documentation map, wrapper script table, and input/output file summary.
