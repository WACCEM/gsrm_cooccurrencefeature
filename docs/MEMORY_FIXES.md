# Memory Leak Fixes for Zarr Writing

## Changes Made to Address Memory Leak at Chunk 94/131

### 1. Enhanced Memory Management in zarr_tools.py

**Improvements to `write_zarr_chunked()` function:**

- **Aggressive Garbage Collection**: Added `gc.collect()` after each chunk processing
- **Memory Debugging**: Optional detailed memory logging with `debug_memory=True`
- **Memory Monitoring**: Real-time memory usage reporting during processing
- **Better Error Handling**: Non-blocking zarr consolidation to prevent hangs
- **Explicit Cleanup**: More thorough deletion of intermediate variables

**Key changes:**
```python
# Force aggressive garbage collection
import gc
gc.collect()

# Optional memory debugging
def log_memory_usage(stage=""):
    if debug_memory:
        import psutil
        memory = psutil.virtual_memory()
        logger.info(f"Memory {stage}: {memory.percent:.1f}% used")
```

### 2. Command Line Options in make_cooccurrence_masks.py

**New debug flag:**
```bash
--debug-memory    Enable detailed memory usage logging
```

**Updated function call:**
```python
successful_times = write_zarr_chunked(
    all_results=all_results,
    mask_variables=mask_variables,
    output_path=output_path,
    coords=ds.coords,
    chunks=chunks,           # Properly defined chunk sizes
    time_chunks=24,          # Process 24 time steps at a time
    logger=logger,
    debug_memory=args.debug_memory  # Enable memory debugging
)
```

### 3. Conservative Memory Test Script

**Created `test_memory_conservative.py`:**
- 16 processing workers (half of Option 1)
- 4 zarr workers (half of Option 1)
- 1 thread per worker
- Memory debugging enabled

**Usage:**
```bash
python test_memory_conservative.py
```

## Memory Leak Prevention Strategy

### Root Cause Analysis
The memory leak at chunk 94/131 was likely caused by:
1. **Incomplete cleanup** of intermediate xarray datasets
2. **Memory fragmentation** from repeated zarr append operations
3. **Zarr consolidation issues** consuming extra memory
4. **Worker memory accumulation** without forced garbage collection

### Solutions Implemented
1. **Explicit garbage collection** after each chunk
2. **Memory usage monitoring** to detect leaks early
3. **Conservative worker scaling** to reduce memory pressure
4. **Non-blocking zarr operations** to prevent hangs
5. **Real-time memory reporting** for debugging

### Expected Results
- **Reduced memory accumulation** between chunks
- **Earlier detection** of memory pressure
- **Better recovery** from temporary memory spikes
- **Completion** of full 131-chunk dataset

## Testing Protocol

1. **Run conservative test first:**
   ```bash
   python test_memory_conservative.py
   ```

2. **Monitor memory usage patterns:**
   - Watch for gradual memory increases
   - Check if memory is released after each chunk
   - Verify completion of all chunks

3. **Scale up gradually:**
   - If conservative test passes, try Option 1 with debug enabled
   - Adjust worker counts based on memory patterns

## Next Steps If Issues Persist

1. **Reduce chunk size:** From 24 to 12 time steps
2. **Further reduce workers:** From 4 to 2 zarr workers  
3. **Implement memory pressure monitoring:** Auto-scale workers down if memory > 90%
4. **Split zarr writing:** Write multiple smaller zarr files and merge later