# Variable Scaling Configuration

This document explains the variable scaling system used to standardize units across different model sources.

## Overview

Different climate models and reanalysis products use different units for the same physical quantities. To enable meaningful comparisons, we apply source-specific scaling factors to convert all variables to standardized units.

## Configuration File

The scaling configuration is defined in `variable_scaling_config.py`, which contains:

- **VARIABLE_SCALING**: Main configuration dictionary
- **Helper functions**: `get_scaling_info()`, `get_all_scaled_variables()`, `get_all_sources()`
- **Pattern matching**: Supports wildcard patterns for variable groups

## Configuration Structure

```python
VARIABLE_SCALING = {
    'source_name': {
        'variable_pattern': (scale_factor, target_units, description),
        ...
    },
    ...
}
```

Where:
- **variable_pattern**: Exact variable name OR pattern with wildcard (e.g., `"zg*"`)
- **scale_factor**: Multiply the variable by this value
- **target_units**: Standardized units after scaling
- **description**: Brief explanation of the conversion

## Pattern Matching

The configuration supports **wildcard patterns** using `*` to match multiple variables:

- **Exact match**: `"pr"` matches only `pr`
- **Pattern match**: `"zg*"` matches `zg850`, `zg500`, `zg_500hPa`, `zg_850hPa`, etc.
- **Precedence**: Exact matches take precedence over pattern matches

### Pattern Matching Rules

1. Patterns must end with `*` (e.g., `"zg*"`, `"ta*"`)
2. The pattern matches any variable **starting with** the prefix
3. If both exact and pattern matches exist, exact match is used
4. Pattern matching is safe - only intended variables are matched

### Example

```python
VARIABLE_SCALING = {
    'scream': {
        'pr': (3600000.0, 'mm h-1', 'Convert precipitation'),
        'zg*': (1.0 / 9.81, 'm', 'Convert geopotential'),  # Matches zg850, zg500, etc.
    }
}
```

This will match:
- `pr` → uses exact match
- `zg850`, `zg500`, `zg_500hPa` → all use the `zg*` pattern

## Currently Configured Variables

### Precipitation (pr)

| Source | Original Units | Scale Factor | Target Units | Notes |
|--------|---------------|--------------|--------------|-------|
| SCREAM | kg m⁻² s⁻¹ | 3,600,000 | mm h⁻¹ | ×1000 (kg→mm) × 3600 (s→h) |
| ERA5   | mm h⁻¹ | 1.0 | mm h⁻¹ | Already in target units |
| NICAM  | kg m⁻² s⁻¹ | 3,600 | mm h⁻¹ | ×1000 (kg→mm) × 3.6 (s→h) |
| ICON   | kg m⁻² s⁻¹ | 3,600 | mm h⁻¹ | ×1000 (kg→mm) × 3.6 (s→h) |
| CESM2  | kg m⁻² s⁻¹ | 3,600 | mm h⁻¹ | ×1000 (kg→mm) × 3.6 (s→h) |
| UM     | kg m⁻² s⁻¹ | 3,600 | mm h⁻¹ | ×1000 (kg→mm) × 3.6 (s→h) |

### Geopotential Height (zg*)

**Pattern**: `zg*` matches all variables starting with `zg` (e.g., `zg850`, `zg500`, `zg_500hPa`, `zg_850hPa`)

| Source | Original Units | Scale Factor | Target Units | Notes |
|--------|---------------|--------------|--------------|-------|
| SCREAM | m² s⁻² (geopotential) | 1/9.81 | m | Divide by gravity |
| ERA5   | m | 1.0 | m | Already geopotential height |
| NICAM  | m² s⁻² (geopotential) | 1/9.81 | m | Divide by gravity |
| ICON   | m² s⁻² (geopotential) | 1/9.81 | m | Divide by gravity |
| CESM2  | m | 1.0 | m | Already geopotential height |
| UM     | m | 1.0 | m | Already geopotential height |

## Usage

### In Python Scripts

```python
from variable_scaling_config import VARIABLE_SCALING

def standardize_variable_units(ds, source):
    """Apply scaling to dataset variables."""
    if source not in VARIABLE_SCALING:
        return ds
    
    scaling_config = VARIABLE_SCALING[source]
    
    for var_name, (scale_factor, target_units, description) in scaling_config.items():
        if var_name in ds.data_vars:
            # Apply scaling
            ds[var_name] = ds[var_name] * scale_factor
            
            # Update attributes
            ds[var_name].attrs['units'] = target_units
            ds[var_name].attrs['scaling_applied'] = description
    
    return ds
```

### View Configuration

Run the config file directly to see all configured scaling:

```bash
python variable_scaling_config.py
```

Output:
```
Variable Scaling Configuration
================================================================================

Source: SCREAM
--------------------------------------------------------------------------------
  pr              :    ×3.6e+06  →  mm h-1      (Convert kg m-2 s-1 to mm/h...)
  zg_500hPa       :   ×0.101937  →  m           (Convert geopotential...)
  ...
```

## Adding New Variables or Sources

### Add a new variable to existing source:

```python
VARIABLE_SCALING = {
    'scream': {
        'pr': (3600000.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'),
        'zg*': (1.0 / 9.81, 'm', 'Convert geopotential'),
        # Add new exact variable
        'tas': (1.0, 'K', 'Already in Kelvin'),
        # Add new pattern for temperature at pressure levels
        'ta*': (1.0, 'K', 'Temperature already in Kelvin'),
    },
    ...
}
```

### Add a new source:

```python
VARIABLE_SCALING = {
    ...
    'new_model': {
        'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h'),
        'zg*': (1.0, 'm', 'Already in geopotential height'),
    },
}
```

### Pattern Matching Best Practices

**Safe patterns** (recommended):
- `"zg*"` - geopotential at all levels (zg850, zg500, zg_500hPa)
- `"ta*"` - temperature at all levels (ta850, ta500, tas)
- `"ua*"` / `"va*"` - wind components at all levels

**Potentially risky** (use with caution):
- `"pr*"` - could match `pr`, `prc` (convective), `prl` (large-scale) if different scaling needed
- Short prefixes like `"t*"` - too broad, could match unintended variables

**Best practice**: Use specific prefixes and always test with `python variable_scaling_config.py`

## Integration with create_etc_composites.py

The `create_etc_composites.py` script automatically applies scaling:

```python
# Load data
ds = load_etc_zarr_data(zarr_file)

# Standardize variable units (automatic)
ds = standardize_variable_units(ds, source)

# Continue with compositing...
```

The scaling is applied **before** creating composites, ensuring all output files have standardized units.

## Benefits

1. **Consistency**: All model outputs use the same units
2. **Traceability**: Original units and scale factors stored in attributes
3. **Maintainability**: Single configuration file to update
4. **Reusability**: Same config can be used by multiple scripts
5. **Documentation**: Clear description of each conversion

## Verification

After scaling, check variable attributes:

```python
import xarray as xr

ds = xr.open_dataset('composite_file.nc')
print(ds['pr'].attrs)

# Output:
# {
#     'units': 'mm h-1',
#     'scaling_applied': 'Convert kg m-2 s-1 to mm/h',
#     'original_units': 'kg m-2 s-1',
#     'scale_factor': 3600000.0
# }
```

## Common Conversions

### Precipitation

- **kg m⁻² s⁻¹ → mm h⁻¹**
  - 1 kg m⁻² = 1 mm (water density = 1000 kg/m³)
  - 1 s = 1/3600 h
  - Factor: 1000 × 3600 = 3,600,000

### Geopotential → Geopotential Height

- **m² s⁻² → m**
  - Φ = g × Z (geopotential = gravity × height)
  - Z = Φ / g
  - g = 9.81 m/s²
  - Factor: 1/9.81 ≈ 0.10194

## Notes

- If a variable is not in the configuration, no scaling is applied
- Scaling is applied in-place to the dataset
- Original units are preserved in variable attributes
- The config file has fallback handling if import fails
