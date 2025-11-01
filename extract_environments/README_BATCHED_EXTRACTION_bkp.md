# Batched ETC 2D Variable Extraction

## Overview

This document explains the optimized batched extraction workflow used in `extract_etc_2d_vars.py` to efficiently extract 2D environmental variables around extratropical cyclone (ETC) tracks from HEALPix model output.

## Performance Summary

**Test case:** 18,225 storm positions, 1 year of data (August 2019 - August 2020)

| Metric | Old Approach | Batched Approach | Speedup |
|--------|-------------|------------------|---------|
| **Total time** | 75.08 minutes | 4.60 minutes | **16.3x faster** |
| **Processing rate** | 4.0 storms/s | 66.0 storms/s | **16.5x faster** |
| **Data loads** | 18,225 time slices | 1,257 time slices | **14.5x fewer** |

**Hardware:** NERSC Perlmutter, local scratch storage (fast disk I/O)

---

## The Problem: Redundant Data Loading

### Old Approach (Inefficient)
```
For each storm (18,225 storms):
    1. Load ENTIRE global HEALPix grid at storm time (786,432 cells)
    2. Extract 81×81 box around storm (~6,500 cells)
    3. Repeat for next storm...
```

**Issues:**
- Loading 786k cells × 18,225 storms = massive redundant I/O
- Many storms share the same timestamp but each loads data independently
- Most loaded data (~99%) is immediately discarded

---

## The Solution: Batched Time-Slice Loading

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

## Performance Breakdown

### Time Reduction Analysis

**Test case:** 18,225 storm positions over 1 year

1. **Batching reduction:**
   - Old: 18,225 time selections
   - New: 1,257 unique time selections
   - **Speedup:** 14.5x fewer data loads

2. **Spatial subsetting:**
   - Old: Load 786,432 cells per time slice
   - New: Load ~30,000 cells per time slice
   - **Speedup:** ~26x less data per time slice

3. **Combined theoretical speedup:** 14.5 × 26 = **377x**

4. **Actual speedup:** **16.3x**
   - Why not 377x? Overhead from:
     - Pixel index calculation
     - Numpy array operations
     - File I/O latency
     - Memory allocation

### Memory Efficiency

**Old approach:**
- Peak memory: 786k cells × 4 bytes = 3.1 MB per time slice
- But: Loaded 18,225 times sequentially

**New approach:**
- Peak memory: ~30k cells × 4 bytes × 15 storms = ~1.8 MB per time batch
- Output array: 18,225 × 81 × 81 × 4 bytes = 480 MB (allocated once)
- **Total steady state:** ~500 MB

### Storage Efficiency

**Output file:**
- Format: Zarr with Blosc/zstd compression
- Size: 0.24 GB per variable (18,225 storm positions)
- Compression ratio: ~3x (vs uncompressed)
- Chunking: (1000, 81, 81) for efficient time-slicing

---

## Workflow Diagram

```
┌──────────────────────────────────────────────────────────────────┐
│                    INPUT: ETC Track File                         │
│  18,225 storm positions with (storm_id, lon, lat, time)          │
└──────────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────────┐
│              STEP 1: Parse and Group by Timestamp                │
│                                                                  │
│  Storm DataFrame  →  Time Groups                                 │
│  ┌──────────────┐    ┌─────────────────────────────┐             │
│  │ ID  Lon  Lat │    │ 2020-05-01 00:00 → 15 storms│             │
│  │ 1   180  45  │    │ 2020-05-01 03:00 → 14 storms│             │
│  │ 2   175  46  │    │ 2020-05-01 06:00 → 16 storms│             │
│  │ 3   170  47  │    │        ...                  │             │
│  │ ... ... ...  │    │ 1,257 unique timestamps     │             │
│  └──────────────┘    └─────────────────────────────┘             │
└──────────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────────┐
│        STEP 2: For Each Timestamp, Load Spatial Subset           │
│                                                                  │
│  ┌────────────────────────────────────────────────────┐          │
│  │  HEALPix Zarr Data (786,432 cells, 3,169 times)    │          │
│  │                                                    │          │
│  │  Time: 2020-05-01 00:00                            │          │
│  │  ┌──────────────────────────────────────────┐      │          │
│  │  │ Global grid (only load needed cells)     │      │          │
│  │  │                                          │      │          │
│  │  │   ◯ ◯ ◯   Storm locations at this time   │      │          │
│  │  │      ◯                                   │      │          │
│  │  │        ◯ ◯                               │      │          │
│  │  │                                          │      │          │
│  │  │  Load union of all storm boxes (~30k)    │      │          │
│  │  └──────────────────────────────────────────┘      │          │
│  └────────────────────────────────────────────────────┘          │
│                          ↓                                       │
│       variable_data.sel(time=T).isel(cell=pixels).compute()      │
└──────────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────────┐
│          STEP 3: Extract All Storms from Cached Subset           │
│                                                                  │
│  Cached Time Slice (30k cells in memory)                         │
│  ┌─────────────────────────────────┐                             │
│  │  cell_id: [157, 892, 1045, ...] │                             │
│  │  values:  [12.3, 8.7, 15.2, ...]│                             │
│  └─────────────────────────────────┘                             │
│           ↓           ↓           ↓                              │
│     ┌────────┐   ┌────────┐   ┌────────┐                         │
│     │Storm 1 │   │Storm 2 │   │Storm 3 │  ... (15 storms)        │
│     │ 81×81  │   │ 81×81  │   │ 81×81  │                         │
│     │  grid  │   │  grid  │   │  grid  │                         │
│     └────────┘   └────────┘   └────────┘                         │
│                                                                  │
│     Fast numpy indexing (microseconds per storm)                 │
└──────────────────────────────────────────────────────────────────┘
                              ↓
┌──────────────────────────────────────────────────────────────────┐
│              STEP 4: Save to Zarr with Chunking                  │
│                                                                  │
│  Output Array: [18,225 storms, 81 y, 81 x]                       │
│  ┌──────────────────────────────────────────┐                    │
│  │ Dimensions: (time, y, x)                 │                    │
│  │ Chunks: (1000, 81, 81)                   │                    │
│  │ Compression: Blosc/zstd, level 3         │                    │
│  │ Size: 0.24 GB per variable               │                    │
│  │                                          │                    │
│  │ Coordinates: x, y (relative to center)   │                    │
│  │ Metadata: storm_id, grid_id, lat, lon    │                    │
│  └──────────────────────────────────────────┘                    │
└──────────────────────────────────────────────────────────────────┘
```

---

## Usage Example

```bash
# Extract multiple 2D variables for full year
python extract_etc_2d_vars.py \
  --catalog_model scream_ne120 \
  --catalog_params '{"zoom": 8}' \
  --trackfile /path/to/etc_stitched_nodes.txt \
  --output_dir /path/to/output \
  --variable ua850 va850 ta850 rh850 ps \
  --start_date 2020-01-01 \
  --end_date 2020-12-31 \
  --radius 10.0 \
  --lon_res 0.25 \
  --lat_res 0.25

# Expected time: ~5 minutes per variable
# Total: 5 variables × 5 min = 25 minutes
```

---

## Key Takeaways

### ✅ What Makes This Fast

1. **Batching by time:** Load each timestamp once, not once per storm (14.5x reduction)
2. **Spatial subsetting:** Load only needed cells, not entire globe (26x reduction)
3. **Deferred computation:** Use Dask lazy evaluation, compute after spatial selection
4. **Numpy indexing:** Extract from cached arrays using integer indexing (microseconds)
5. **Memory efficient:** ~500 MB peak, scales to millions of storms

### 📊 When to Use This Approach

- ✅ **Multiple events at same times:** MCS, fronts, tropical cyclones
- ✅ **Spatial extraction:** Need 2D/3D grids around events
- ✅ **Large datasets:** Thousands to millions of events
- ✅ **HEALPix grids:** Unstructured meshes benefit from spatial indexing

### 🚀 Potential Extensions

1. **Threading (2-3x speedup):** Add ThreadPoolExecutor for I/O-bound remote data
2. **Distributed (10x+ speedup):** Use Dask distributed for multi-node processing
3. **3D variables:** Extend to extract vertical profiles (pressure levels)
4. **Multiple variables:** Load multiple variables per time slice simultaneously

---

## References

- **Similar implementation:** `get_env_vars.py` (MCS circular area extraction)
- **HEALPix library:** `healpy` for pixel indexing
- **Grid tools:** `easygems.healpix` for coordinate transformations
- **Performance test:** NERSC Perlmutter, 18,225 storms, 1 year, local storage

---

**Author:** Zhe Feng  
**Last updated:** October 31, 2025  
**Optimization:** 16.3x speedup via batched time-slice loading
