# Spatial Storm Type Attribution for Extreme Precipitation

This module calculates the spatial distribution of extreme precipitation events and attributes them to different storm types using Dask for parallel processing. Storm type contributions are quantified as precipitation fractions: the sum of extreme precipitation from each storm type divided by the total extreme precipitation at each grid cell.

## Scripts

### 1. `calc_stormtype_extreme_precip_spatial.py`

Main Python script that processes precipitation and storm mask data to attribute extreme precipitation to specific storm types.

**Features:**
- Dask-based parallelization over the time dimension
- Each time step processed independently for maximum efficiency
- Supports multiple percentile thresholds (P90, P95, P99, etc.)
- Outputs per-cell precipitation fractions and occurrence counts for each storm type
- Calculates fractions based on actual precipitation amounts, not occurrence counts

**Usage:**
```bash
python calc_stormtype_extreme_precip_spatial.py \
    --catalog_source scream_ne120 \
    --percentiles P90 P95 P99 \
    --n_workers 16 \
    --compute_cloud_types
```

**Arguments:**
- `--catalog_source`: Data source name (e.g., scream_ne120, IR_IMERG)
- `--config_file`: Path to configuration YAML (default: config/config_sources.yaml)
- `--percentiles`: Percentiles to process (default: P90)
- `--start_date`: Start date in YYYY-MM-DD format (optional)
- `--end_date`: End date in YYYY-MM-DD format (optional)
- `--output_dir`: Output directory (default: /pscratch/sd/w/wcmca1/hackathon/extreme_precip)
- `--n_workers`: Number of Dask workers (default: 8)
- `--compute_cloud_types`: Include cloud type attribution (default: True)
- `--skip_cloud_types`: Skip cloud type computation

### 2. `run_stormtype_spatial.sh`

Bash wrapper script for running the Python script on multiple data sources.

**Usage:**
```bash
# Process a single source
bash run_stormtype_spatial.sh scream_ne120

# Process all sources (when no argument provided)
bash run_stormtype_spatial.sh
```

### 3. `slurm_stormtype_extreme_precip.sh `

SLURM batch script for running on NERSC compute nodes with optimized resources.

**Usage:**
```bash
# Submit job for SCREAM data with all percentiles
sbatch slurm_stormtype_extreme_precip.sh  scream_ne120 "P90 P95 P99"

# Submit job for IMERG data with P90 only
sbatch slurm_stormtype_extreme_precip.sh  IR_IMERG "P90"

# Submit with all default sources & arguments (DEFAULT_SOURCES, P90 P95)
sbatch slurm_stormtype_extreme_precip.sh 
```

**SLURM Configuration:**
- 1 node, 16 CPUs
- 20 min walltime (sufficient for computing two percentiles)
- Regular QOS
- CPU constraint
- Logs saved to `logs/stormtype_spatial_<jobid>.{out,err}`

## Storm Type Categories

The script attributes extreme precipitation to 13 storm type categories:

**Isolated Features (4 types):**
1. MCS (isolated)
2. AR (isolated)
3. ETC (isolated)
4. TC (Tropical Cyclone)

**Co-occurring Features (4 types):**
5. MCS+AR (2-way)
6. MCS+ETC (2-way)
7. AR+ETC (2-way)
8. MCS+AR+ETC (3-way)

**Cloud Types (4 types):**
9. DC (Deep Convection)
10. ND (Non-deep Convection)
11. ST (Stratiform)
12. DZ (Drizzle)

**Untracked:**
13. Unassigned

## Priority Assignment Logic

To avoid double-counting, cells are assigned using priority order:
1. 3-way co-occurrence (MCS+AR+ETC)
2. 2-way co-occurrences (MCS+AR, MCS+ETC, AR+ETC)
3. Isolated features (MCS, AR, ETC)
4. TC
5. Cloud types (DC, ND, ST, DZ)
6. Unassigned (remaining extreme precipitation)

## Precipitation Fraction Calculation

Storm type contributions are calculated as **precipitation fractions** rather than occurrence fractions:

```
pr_fraction = sum(stormtype_extreme_precip) / sum(all_extreme_precip)
```

For each grid cell:
- At each time step, precipitation values where `pr > threshold` are accumulated for each storm type
- The fraction represents the proportion of total extreme precipitation attributed to each storm type
- This provides a physically meaningful measure of storm type importance based on actual precipitation amounts

## Output Format

**NetCDF file with variables:**

*Totals:*
- `total_extreme_count`: Total count of extreme precipitation occurrences at each cell
- `total_extreme_precip`: Sum of all extreme precipitation amounts (mm/h) at each cell

*Per-storm-type variables (13 storm types):*
- `<stormtype>_count`: Count of extreme precipitation occurrences attributed to this storm type
- `<stormtype>_frac`: Precipitation fraction attributed to this storm type

*Coordinates:*
- `lat`, `lon`: HEALPix cell coordinates (from mask dataset)

**Variable details:**
- Counts: int32, compressed (number of time steps exceeding threshold)
- Fractions: float32, values between 0-1 (precipitation-weighted contributions)
- Total precipitation: float32 (mm/h accumulated over all extreme events)
- Shape: `(cell,)` for all variables

**Fraction calculation:**
Each `<stormtype>_frac` is calculated as:
```
stormtype_frac = sum(pr values where pr > threshold and stormtype mask) / 
                 sum(pr values where pr > threshold)
```

**Filename convention:**
```
# With specified date range
<source>_stormtype_spatial_<percentile>_<start_date>_<end_date>.nc

# Without date range (full dataset)
<source>_stormtype_spatial_<percentile>.nc
```

Examples:
```
scream_stormtype_spatial_p90_20200101_20201231.nc
IR_IMERG_stormtype_spatial_p95.nc
```

**Note:** Dominant storm type is **not** computed in this script. For ranking and dominant type analysis, see the analysis notebook `plot_cof_extreme_raintype_rank_map.ipynb`.

## Performance Optimization

**Dask Parallelization:**
- Each time step processed independently as a delayed task
- Both counts and precipitation amounts tracked separately
- Results aggregated after parallel computation
- Progress bar shows real-time status

**Memory Management:**
- Time steps processed in parallel, not loaded all at once
- Template arrays computed once to avoid dask memory issues
- Results accumulated incrementally as numpy arrays
- Optimized for NERSC Perlmutter nodes

**Scaling:**
- Full year (~1460-1588 time steps): ~5 minutes per percentile with 16 workers
- Two percentiles (P90, P95): ~10 minutes total processing time with 16 workers
- Multiple percentiles: processed sequentially within same job

**Observed Performance (16 CPUs, NERSC Perlmutter):**
| Source | Time Steps | Percentiles | Total Time |
|--------|-----------|-------------|------------|
| IR_IMERG | 1588 | P90, P95 | ~9 min |
| scream_ne120 | 1577 | P90, P95 | ~12 min |
| icon_d3hp003 | 1460 | P90, P95 | ~14 min |
| um_glm_n2560_RAL3p3 | 1576 | P90, P95 | ~9 min |
| nicam_gl11 | 1459 | P90, P95 | ~10 min |

Processing rate: ~400-450 time steps per minute with 16 workers

## Example Workflow

### Interactive Testing (Login Node)
```bash
# Test with a subset
python calc_stormtype_extreme_precip_spatial.py \
    --catalog_source scream_ne120 \
    --percentiles P90 \
    --start_date 2020-01-01 \
    --end_date 2020-01-31 \
    --n_workers 4
```

### Production Run (Compute Node)
```bash
# Create logs directory
mkdir -p logs

# Submit full year processing
sbatch slurm_stormtype_extreme_precip.sh  scream_ne120 "P90 P95 P99"

# Monitor job
squeue -u $USER
tail -f logs/stormtype_spatial_<jobid>.out
```

### Batch Processing (All Sources)
```bash
# Submit jobs for all data sources & default percentiles (P90 P95)
sbatch slurm_stormtype_extreme_precip.sh

# Submit jobs for all data sources for 1 percentile
sbatch slurm_stormtype_extreme_precip.sh "P95"

```

## Output Analysis

See analysis notebook `plot_cof_extreme_raintype_rank_map.ipynb`

## Requirements

**Python packages:**
- xarray
- numpy
- dask
- easygems.healpix
- intake
- pyyaml

**Data requirements:**
- Precipitation data (6-hourly, HEALPix zoom 8)
- Storm mask data (COF masks)
- Extreme precipitation thresholds (from `calc_extreme_precip_thresholds.py`)

## Troubleshooting

**Issue: Out of memory**
- Reduce `--n_workers`
- Process smaller time ranges with `--start_date` and `--end_date`

**Issue: Slow performance**
- Increase `--n_workers` (up to available CPUs)
- Use compute nodes instead of login nodes
- Check data locality (zarr stores should be on fast filesystem)

**Issue: Missing thresholds**
- Run `calc_extreme_precip_thresholds.py` first to generate percentile files
- Check that percentile names match (e.g., pr_p90 vs pr_P90)

**Issue: Unexpected fraction values**
- Fractions should sum to ≤1.0 at each cell (unassigned accounts for remainder)
- Check that storm masks have proper priority assignment
- Verify that cloud types are excluded from tracked storm regions

## Contact

Zhe Feng  
Email: zhe.feng@pnnl.gov  
Date: December 2025
