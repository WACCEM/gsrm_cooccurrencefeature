# MCS Swath Mask Processing - Resume Guide

## Performance Summary

**Verified Performance (32 workers, 100 batch-size):**
- **IMERG**: 427 seconds (7.1 minutes)
- **ICON**: 526 seconds (8.8 minutes)

**Quick Start for Production:**
```bash
# Maximum speed (dedicated node)
python make_mcs_swath_masks.py -c config.yml --workers 32 --batch-size 100

# Balanced (most cases)
python make_mcs_swath_masks.py -c config.yml --workers 16 --batch-size 100

# Conservative (shared node)
python make_mcs_swath_masks.py -c config.yml --workers 8 --batch-size 50
```

---

## Quick Commands

### 1. Check for Missing Chunks
```bash
python make_mcs_swath_masks.py -c ../config/config_mcs_tbpf_scream_healpix9.yml --check-missing
```

### 2. Resume Processing Missing Chunks
```bash
python make_mcs_swath_masks.py -c ../config/config_mcs_tbpf_scream_healpix9.yml --resume --workers 8 --batch-size 50
```

### 3. Full Processing (from scratch)
```bash
python make_mcs_swath_masks.py -c ../config/config_mcs_tbpf_scream_healpix9.yml --workers 8 --batch-size 50
```

## Memory Management

### Understanding Your Hardware
- **Your Node**: 2x AMD EPYC 7763 (128 cores, 512 GB RAM)
- **Optimal Config (verified)**: 32 workers × ~16 GB/worker = 512 GB (full utilization)
- **Safe Config**: 16 workers × ~32 GB/worker = comfortable headroom

### Memory Optimization Achieved
After implementing multiple optimizations, memory usage is now predictable and manageable:

1. **Task graph**: Reduced from 55 GB → <1 MB (using integer indices)
2. **Worker overhead**: Minimized by opening zarr with `chunks=None`
3. **Memory cleanup**: Explicit `gc.collect()` after each batch
4. **Removed unnecessary computation**: Only processing mcs_mask (not ccs_mask)

### Memory Usage Guidelines

**32 workers (maximum performance):**
- Memory per worker: ~16 GB (512 GB / 32)
- Total allocation: ~512 GB (full node)
- Best for: Dedicated nodes, production runs
- Verified stable with 100 batch-size

**16 workers (balanced):**
- Memory per worker: ~32 GB (512 GB / 16)
- Total allocation: ~512 GB with headroom
- Best for: Standard production, some shared usage
- Very stable configuration

**8 workers (conservative):**
- Memory per worker: ~64 GB (512 GB / 8)
- Total allocation: ~512 GB with large safety margin
- Best for: Shared nodes, memory-intensive datasets
- Guaranteed stability

### Recommended Configuration

#### For Maximum Speed (Dedicated Node)
```bash
python make_mcs_swath_masks.py \
  -c ../config/config_mcs_tbpf_scream_healpix9.yml \
  --workers 32 \
  --batch-size 100
```

**Benefits:**
- **Fastest processing**: 7-9 minutes for full dataset (verified on IMERG/ICON)
- Maximum CPU utilization on dedicated node
- Best throughput: ~2-3 chunks/second
- **Use when**: You have dedicated access to the node

**Watch for:**
- Memory usage should stay below 85% per worker
- If workers restart, reduce to 16 workers

#### Option 1: Balanced Performance (Recommended for Most Cases)
```bash
python make_mcs_swath_masks.py \
  -c ../config/config_mcs_tbpf_scream_healpix9.yml \
  --workers 16 \
  --batch-size 100 \
  --resume  # if resuming
```

**Benefits:**
- Fast processing: ~12-15 minutes for full dataset
- 16 workers × ~32 GB/worker = safe memory allocation
- Good balance between speed and stability
- **Use when**: Standard production runs

#### Option 2: Memory-Safe (For Shared Nodes or Large Datasets)
```bash
python make_mcs_swath_masks.py \
  -c ../config/config_mcs_tbpf_scream_healpix9.yml \
  --workers 8 \
  --batch-size 50 \
  --resume  # if resuming
```

**Benefits:**
- 8 workers × ~64 GB/worker = plenty of memory headroom
- Stable processing even on shared nodes
- Runtime: ~20-30 minutes for full dataset
- **Use when**: Sharing compute resources or memory-intensive datasets

## What Changed to Improve Performance

### 1. Removed CCS Mask Computation
**Impact**: ~40% faster processing

The script originally computed both MCS and CCS swath masks. Since only MCS masks are needed:
- **Before**: Processing both mcs_mask and ccs_mask
- **After**: Only processing mcs_mask
- **Result**: IMERG in 7.1 min, ICON in 8.8 min (32 workers, 100 batch-size)

### 2. Fixed Zarr Writing Bug
**Impact**: Eliminated write errors, enabled reliable processing

Fixed mismatch between aggregation window and zarr writing:
- **Before**: Passing 6 input times to zarr writer (causing shape mismatch errors)
- **After**: Passing 1 output time per aggregated swath
- **Result**: No more "expected array with shape (0, 786432), got (6, 786432)" errors

### 3. Reduced Task Graph Size
**Impact**: Reduced memory overhead from 55 GB to <1 MB

Changed worker communication to use integer indices:
- **Before**: 55 GB task graph (passing time arrays to workers)
- **After**: <1 MB task graph (passing only start_idx, end_idx)
- **Result**: Workers can be scheduled instantly without memory bottleneck

### 4. Added Explicit Garbage Collection
**Impact**: Better memory management, fewer worker restarts

```python
import gc
# ... after processing ...
del chunk_results
gc.collect()  # Force Python to free memory
```

### 5. Optimized Worker Data Loading
**Impact**: Lower memory overhead per worker

```python
# Workers open zarr with chunks=None to avoid overhead
ds = xr.open_dataset(zarr_path, engine='zarr', chunks=None)
```

### 6. Implemented Batched Submission
**Impact**: Scheduler stability, ability to process large datasets

- **Before**: Submit all 1576 chunks at once → scheduler overwhelmed
- **After**: Submit 100 chunks at a time → stable processing
- **Result**: Reliable completion without scheduler crashes

### 7. Separated Aggregation and Storage Chunking
**Impact**: Clearer code, optimized storage access

- **`aggregation_window`**: Processing parameter (e.g., 6 hours → 1 swath)
- **`zarr_chunk_size_time`**: Storage optimization (28 swaths = 1 week)
- **Result**: User can change aggregation window without affecting storage layout

## Resume Workflow Example

### Step 1: Check What's Missing
```bash
python make_mcs_swath_masks.py -c config.yml --check-missing
```

**Example Output:**
```
⚠️  Found 128 missing chunks out of 1576 total
Missing chunk indices: [1449, 1450, 1451, ..., 1576]

To resume processing, run:
python make_mcs_swath_masks.py -c config.yml --resume --workers 8
```

### Step 2: Resume Processing
```bash
python make_mcs_swath_masks.py -c config.yml --resume --workers 8 --batch-size 50
```

**What Happens:**
- Script checks zarr file for missing/incomplete chunks
- Only processes those specific chunks (e.g., 1449-1576)
- Writes directly to existing zarr file
- Preserves all previously completed chunks

### Step 3: Verify Completion
```bash
python make_mcs_swath_masks.py -c config.yml --check-missing
```

**Expected Output:**
```
✅ All 1576 chunks are complete!
```

## Understanding the Processing

### Chunk Organization
- **Input**: 9456 hourly time steps
- **Processing**: Aggregate every 6 hours → 1 swath
- **Output**: 1576 swath masks (one per 6-hour period)

### Batch Processing
With `--batch-size 50`:
- **Batch 1**: Submit chunks 1-50 to workers → Wait for completion → Write results
- **Batch 2**: Submit chunks 51-100 to workers → Wait for completion → Write results
- ...
- **Batch 32**: Submit chunks 1551-1576 to workers → Wait for completion → Write results

### Why Batching Helps
- **Without batching**: Submit all 1576 chunks at once → Scheduler overwhelmed
- **With batching**: Submit 50 chunks at a time → Scheduler can handle it

## Monitoring Progress

### Watch Dask Dashboard
The script prints the dashboard URL when starting:
```
Dask dashboard: http://127.0.0.1:8787/status
```

**What to Monitor:**
- **Task Stream**: Should show continuous activity across all workers
- **Memory**: Should stay below 85% per worker
- **Progress**: Should show steady completion rate

### Log Messages to Watch For

**Good Signs:**
```
BATCH 1/32: Processing chunks 1-50
Chunk 25/1576 complete: 1 swath mask(s) written
Batch 1/32 complete: 50/1576 total chunks processed
```

**Warning Signs:**
```
WARNING - Worker is at 90% memory usage. Pausing worker.
WARNING - Worker exceeded 95% memory budget. Restarting...
ERROR - Task marked as failed because 4 workers died
```

**If You See Warnings:**
- Kill the job
- Reduce number of workers
- Reduce batch size
- Resume with new settings

## Performance Benchmarks

### Actual Performance (32 workers, 100 batch-size)
After removing CCS mask computation (not needed), processing is significantly faster:

- **IMERG**: 427 seconds (7.1 minutes) for full dataset
- **ICON**: 526 seconds (8.8 minutes) for full dataset

**Configuration used:**
```bash
python make_mcs_swath_masks.py \
  -c config.yml \
  --workers 32 \
  --batch-size 100
```

### Estimated Runtime for Different Configurations

#### 32 Workers (Maximum Parallelism)
- **Full dataset**: ~7-9 minutes (verified)
- **Throughput**: ~2-3 chunks/second
- **Memory per worker**: ~16 GB (512 GB / 32)
- **Best for**: Quick production runs, when node is dedicated

#### 16 Workers (Balanced)
- **Full dataset**: ~12-15 minutes (estimated)
- **Throughput**: ~1.5-2 chunks/second
- **Memory per worker**: ~32 GB (512 GB / 16)
- **Best for**: Balance between speed and memory safety

#### 8 Workers (Conservative)
- **Full dataset**: ~20-30 minutes (estimated)
- **Throughput**: ~0.8-1.3 chunks/second
- **Memory per worker**: ~64 GB (512 GB / 8)
- **Best for**: Shared nodes, memory-constrained environments

## Troubleshooting

### Problem: "Large graph" warning still appears
**Solution**: You're using an old version of the script. Make sure you have the latest changes with integer indices.

**Check**: The warning should be <1 MB, not 55 GB.

### Problem: "ValueError: expected array with shape (0, 786432), got (6, 786432)"
**Solution**: This is fixed in the latest version. Make sure you're using the updated code that passes `output_time = np.array([time_coords[start_idx]])` instead of `chunk_times = time_coords[start_idx:end_idx]`.

**Why this happened**: We aggregate 6 hourly times into 1 swath, but were passing 6 times to zarr writer instead of 1.

### Problem: Workers still dying despite memory optimizations
**Solution**: 
1. Check system memory: `free -h`
2. Check if other processes are using memory: `ps aux --sort=-%mem | head -20`
3. Reduce workers: try 16 → 8 → 4
4. Reduce batch size: try 100 → 50 → 25

### Problem: Processing is slower than expected
**Solution**: 
1. Verify you're using 32 workers for maximum speed
2. Increase batch size: `--batch-size 200` (if memory permits)
3. Check Dask dashboard for idle workers (indicates bottleneck)
4. Ensure you're not running on a shared node with other jobs

### Problem: "All chunks complete" but data looks wrong
**Solution**: The check looks for non-zero data. If actual tracks are missing, inspect manually:
```python
import xarray as xr
ds = xr.open_zarr('output.zarr')
print(ds.mcs_mask.isel(time=0).values.sum())  # Should be > 0 if data exists
# Check a few time steps
for i in range(0, 100, 10):
    print(f"Time {i}: {ds.mcs_mask.isel(time=i).values.sum()}")
```

### Problem: Resume mode processes already-complete chunks
**Solution**: The check uses a simple heuristic (all zeros = missing). If your data legitimately has all-zero chunks, you'll need to manually specify chunk indices.

## Performance Optimization Tips

### 1. Use Maximum Workers on Dedicated Nodes
If you have exclusive access to a compute node:
```bash
--workers 32 --batch-size 100
```
This achieves the fastest processing times (7-9 minutes for full dataset).

### 2. Optimize Batch Size
The batch size controls how many chunks are submitted to the Dask scheduler at once:

- **Large batch (100-200)**: Better for faster processing, requires stable memory
- **Medium batch (50-100)**: Balanced approach, good for most cases
- **Small batch (20-50)**: Use if experiencing scheduler issues or memory pressure

**Rule of thumb**: `batch_size ≈ n_workers × 3` gives good overlap without overwhelming the scheduler.

### 3. Process Only What You Need
Use `--aggregation-window` to control the temporal resolution:
```bash
# 6-hourly swaths (default, 4 per day)
--aggregation-window 6

# 12-hourly swaths (2 per day, faster processing)
--aggregation-window 12

# Daily swaths (1 per day, fastest)
--aggregation-window 24
```

Larger aggregation windows = fewer output time steps = faster processing.

### 4. Resume from Failures
Always check for missing chunks before reprocessing:
```bash
# Check what's missing
python make_mcs_swath_masks.py -c config.yml --check-missing

# Resume only missing chunks
python make_mcs_swath_masks.py -c config.yml --resume --workers 32 --batch-size 100
```

### 5. Adjust Zarr Chunking for Your Access Pattern
The script uses `zarr_chunk_size_time = 28` (1 week) by default. If you typically access data at different intervals, you can modify this in the script:

- **Daily analysis**: 28 chunks (1 week) ✓ Default
- **Monthly analysis**: 120 chunks (1 month)
- **Seasonal analysis**: 360 chunks (3 months)

Larger zarr chunks = fewer chunks to read = faster for large time ranges, but slower for small queries.

### 6. Test with Small Subset First
Before processing the full dataset, test with a subset:
```bash
python make_mcs_swath_masks.py \
  -c config.yml \
  --workers 32 \
  --batch-size 100 \
  --test-steps 240  # Test with 240 time steps (10 days)
```

This helps verify configuration and estimate full runtime.

## Support Information

If you continue having issues:
1. Check Dask documentation: https://distributed.dask.org/en/latest/worker-memory.html
2. Monitor system resources: `htop` or `top` during processing
3. Check for other processes consuming memory
4. Consider processing in multiple stages (e.g., chunks 0-500, then 501-1000, etc.)
