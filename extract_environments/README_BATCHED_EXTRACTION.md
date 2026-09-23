# Batched ETC 2D Variable Extraction

## Overview

This document explains the optimized batched extraction workflow used in `extract_etc_2d_vars.py` to efficiently extract 2D environmental variables around extratropical cyclone (ETC) tracks from HEALPix model output.

## Batched Time-Slice Loading

### Key Insight
**Many storms occur at the same time** (10-20 storms per 3-6 hour interval)
→ Load each unique time slice **once** and extract all storms from it

### Optimized Workflow

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: Group Storms by Timestamp                           │
└─────────────────────────────────────────────────────────────┘

Storms:  S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 ...
Times:   T1 T1 T2 T2 T2 T3 T3 T3 T3 T4  T4  T4  ...
              ↓
Group:  T1 → [S1, S2]
        T2 → [S3, S4, S5]
        T3 → [S6, S7, S8, S9]
        T4 → [S10, S11, S12]
        ...


┌─────────────────────────────────────────────────────────────┐
│ Step 2: For Each Unique Timestamp                           │
└─────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────┐
│ HEALPix Grid (786,432 cells)                       │
│                                                    │
│    ┌─────┐                                         │
│    │ S1  │  Storm 1 box (81×81)                    │
│    └─────┘                                         │
│                  ┌─────┐                           │
│                  │ S2  │  Storm 2 box              │
│                  └─────┘                           │
│                                                    │
│                          ┌─────┐                   │
│                          │ S3  │  Storm 3 box      │
│                          └─────┘                   │
│                                                    │
│  Union of all boxes → ~30,000 cells needed         │
└────────────────────────────────────────────────────┘

Instead of loading all 786k cells:
  1. Calculate which HEALPix cells are needed for ALL storms at time T
  2. Load ONLY those cells (typically ~30,000)
  3. Extract all storm boxes from this cached subset


┌─────────────────────────────────────────────────────────────┐
│ Step 3: Spatial Subsetting BEFORE Computation               │
└─────────────────────────────────────────────────────────────┘

❌ OLD (loads everything):
   variable_data.sel(time=T).compute().isel(cell=pixels)
   ↑ Loads 786k cells into memory first

✅ NEW (loads only needed cells):
   variable_data.sel(time=T).isel(cell=pixels).compute()
   ↑ Dask optimizes to load only 30k cells


┌─────────────────────────────────────────────────────────────┐
│ Step 4: Fast Extraction from Cached Subset                  │
└─────────────────────────────────────────────────────────────┘

For each storm at this timestamp:
  ┌──────────────────────────┐
  │ Cached time slice        │
  │ (30k cells in memory)    │
  │                          │
  │  cell_id → array_index   │
  │    157 → 0               │
  │    892 → 1               │
  │   1045 → 2               │
  │    ...                   │
  └──────────────────────────┘
           ↓
  Use numpy indexing (microseconds):
    subset_indices = [[cell_to_idx[c] for c in row] for row in storm.pixels]
    extracted_grid = data_subset.values[subset_indices]
```

---

## Implementation Details

### 1. Storm Grouping by Timestamp

```python
# Create mapping: timestamp → list of storms
time_to_storms = {}
for idx, row in storm_df.iterrows():
    storm_time = pd.Timestamp(row['base_time'])
    if storm_time not in time_to_storms:
        time_to_storms[storm_time] = []
    time_to_storms[storm_time].append((idx, row))

unique_times = sorted(time_to_storms.keys())
# Result: 18,225 storms → ~1,257 unique timestamps (~14.5 storms/timestamp)
```

### 2. Collect Pixel Union for All Storms

```python
# For all storms at this timestamp, collect all needed HEALPix pixels
all_pixels_needed = set()

for storm_idx, storm_row in time_to_storms[storm_time]:
    # Calculate 81×81 grid around storm center
    lon_grid = np.arange(lon_center - radius, lon_center + radius + lon_res, lon_res)
    lat_grid = np.arange(lat_center - radius, lat_center + radius + lat_res, lat_res)
    
    # Convert to HEALPix pixel indices
    pix = hp.ang2pix(nside, *np.meshgrid(lon_grid, lat_grid), nest=True, lonlat=True)
    
    # Add to union
    all_pixels_needed.update(pix.flatten())

# Typically ~20,000-40,000 pixels (vs 786,432 global)
```

### 3. Load Once with Spatial Subset

```python
# CRITICAL: Load only needed cells, defer .compute() until AFTER spatial selection
data_at_time_lazy = variable_data.sel(time=storm_time, method='nearest')
data_subset = data_at_time_lazy.isel(cell=list(all_pixels_needed)).compute()

# Create fast lookup: HEALPix cell ID → array index
cell_to_subset_idx = {cell_id: i for i, cell_id in enumerate(all_pixels_needed)}
```

### 4. Extract All Storms (Fast Numpy Indexing)

```python
for storm_idx, storm_row in time_to_storms[storm_time]:
    # Map storm's HEALPix pixels to indices in cached subset
    subset_indices = [[cell_to_subset_idx[cell] for cell in row] 
                      for row in storm.pix]
    
    # Extract using numpy array indexing (microseconds!)
    extracted_grid = data_subset.values[subset_indices]  # 81×81 array
    
    # Store in output array
    output_array[storm_idx, :, :] = extracted_grid
```

---


---

## Key Takeaways & Technical Notes

### What Makes This Fast

1. **Batching by time:** Load each timestamp once, not once per storm (e.g., 1,260 loads vs 18,225)
2. **Spatial subsetting:** Load only the union of needed cells for all storms at a time (e.g., ~30,000 vs 786,432)
3. **Deferred computation:** Use Dask lazy evaluation, compute after spatial selection (`.isel(cell=pixels).compute()`)
4. **Numpy indexing:** Extract each storm's grid from cached arrays using fast integer indexing
5. **Memory efficient:** ~500 MB peak, scales to millions of storms

### When to Use This Approach

- Multiple events at same times (MCS, ETC, fronts, TCs)
- Need 2D/3D spatial grids around events
- Large datasets (thousands to millions of events)
- HEALPix or other unstructured global grids

### Output Format

- **Format:** Zarr (chunked array storage)
- **Chunking:** (1000, 81, 81) for efficient time-slicing
- **Compression:** Blosc/zstd, clevel=3
- **Size:** ~0.24 GB per variable per year (for 18,225 storms)

### Polar Boundary Handling

- For storms at high latitudes (>80°), grid may extend beyond ±90°
- Invalid latitude rows are padded with NaN
- Maintains consistent 81×81 output shape for all storms

---

## Time Matching (changed 2026-09-20)

Each track point is matched to the source by its **exact** time (the track files are 6-hourly, 00/06/12/18 UTC). Earlier versions used `sel(time, method='nearest')`, which never fails: a time after the end of the record returns the last
frame and a time before the start the first frame. Against the current COF stores, whose periods are shorter than the track periods, that would put a copy of an edge frame at the ETC points outside the COF period (ICON 14.6% of the points, UM 3.2%,
SCREAM 0.5%; the March output has NaN `overlap_flag` at all of them, NaN masks for ICON and partly finite masks for UM and SCREAM). Now:

- a track time without a frame in the source gives a NaN slab; the point keeps its metadata (all variables share one point list) and is counted in the store attributes (`n_points`, `n_points_time_missing`, `n_points_failed`);
- `--max_missing_fraction` (default 0) fails the job *before it writes* when more points are missing than expected (`submit_etc_extraction_jobs.py` passes the expected share per model); a slice that cannot be loaded always fails the job, and so does a variable that cannot be extracted (exit status 1);
- whole time chunks of the source (up to 400 MB) are read once instead of once per storm time (a 48-frame chunk of the COF store was read 48 times): identical results, UM `pr` in 4.4 minutes;
- `--pr_source cof_tot_pr` (default) takes `pr` of the models from Step 1's `tot_pr` in the COF store (`<COF_DATA_ROOT>/cof_masks/<source>_cofmasks_hp8_v1.zarr`), the same [T, T+6 h) window as the COF masks; ERA5 keeps IMERG 6-hourly;
  `--pr_source legacy` reads the separate 6-hourly files and the catalog 6-hour mean as before (SCREAM's file is centred on T, ICON's catalog mean is a full window early; see [run_etc_pipeline.md](../docs/procedures/run_etc_pipeline.md)).

---

## Performance Results (NERSC Perlmutter, local scratch)

- **Old approach:** 75 min for 18,225 storms (serial, redundant loads)
- **Batched approach:** 4.6 min for 18,225 storms (16x speedup)
- **Peak memory:** ~500 MB (output + one time slice)
- **Output size:** 0.24 GB per variable per year

---

## Comparison: ETC vs MCS Extraction

| Aspect              | MCS (get_env_vars.py)         | ETC (extract_etc_2d_vars.py)         |
|---------------------|-------------------------------|--------------------------------------|
| Batching strategy   | Group by time, pixel union    | Group by time, pixel union           |
| Data loading        | .isel(cell=pixels).compute()  | .isel(cell=pixels).compute()         |
| Spatial extraction  | Circular areas (healpy)       | 2D rectangular grids (81×81)         |
| Output              | Statistics (mean, std, etc.)  | Full 2D spatial fields               |
| Speedup             | ~10-20x                       | ~16x                                 |

---

## Usage Example

```bash

```bash# Extract multiple 2D variables for full year

# Extract multiple variables for full yearpython extract_etc_2d_vars.py \

python extract_etc_2d_vars.py \

  --catalog_model scream_ne120 \

  --catalog_params '{"zoom": 8}' \

  --trackfile /path/to/etc_track.txt \

  --output_dir /path/to/output \

  --variable ua850 va850 ta850 rh850 ps \

  --start_date 2020-01-01 \

  --end_date 2020-12-31 \

  --radius 10.0 \

  --lon_res 0.25 \

  --lat_res 0.25

# Expected time: ~5 minutes per variable

# Expected time: ~5 variables × 4.6 min = 23 minutes# Total: 5 variables × 5 min = 25 minutes

``````



------


### 🚀 Potential Extensions

### Output Format

- **Format:** Zarr (chunked array storage)1. **Threading (2-3x speedup):** Add ThreadPoolExecutor for I/O-bound remote data

- **Chunking:** (1000, 81, 81) - efficient for time-slicing2. **Distributed (10x+ speedup):** Use Dask distributed for multi-node processing

- **Compression:** Blosc/zstd, clevel=33. **3D variables:** Extend to extract vertical profiles (pressure levels)

- **Size:** ~0.03 GB per variable per month4. **Multiple variables:** Load multiple variables per time slice simultaneously


---

## Acknowledgments

Optimization strategy inspired by Laura Paccini's `get_env_vars.py` implementation for MCS track extraction. The batched time-slice loading approach was adapted from her efficient circular area statistics extraction.

---

**Author:** Zhe Feng  
**Date:** November 11, 2025  
**Performance tested on:** NERSC Perlmutter (local scratch storage)
