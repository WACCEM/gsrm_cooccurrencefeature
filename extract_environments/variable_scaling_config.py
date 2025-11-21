"""
Variable scaling configuration for standardizing units across different model sources.

This module defines scaling factors to convert variables from their native units
to standardized units for cross-model comparison.

Structure: {source: {variable_pattern: (scale_factor, target_units, description)}}

where:
- variable_pattern: exact variable name OR pattern with wildcard (e.g., "zg*" matches zg850, zg500, etc.)
- scale_factor: multiply the variable by this value
- target_units: standardized units after scaling
- description: brief explanation of the conversion

Pattern Matching Rules:
- Exact names (no wildcard) match only that specific variable
- Patterns ending with * match any variable starting with that prefix (e.g., "zg*" matches "zg850", "zg500")
- Exact matches take precedence over pattern matches

Author: Zhe Feng | zhe.feng@pnnl.gov
Date: November 2025
"""

# =============================================================================
# VARIABLE RENAMING CONFIGURATION
# =============================================================================

VARIABLE_RENAMING = {
    'scream': {
        # Example: SCREAM uses different naming conventions
        # 'original_name': 'standard_name'
        'zg500': 'zg_500hPa',
        'rh850': 'rh_850hPa',
    },
    
    'era5': {
        # ERA5 variable renaming
    },
    
    'nicam': {
        # NICAM variable renaming
    },
    
    'icon': {
        # ICON variable renaming
    },
    
    'cesm2': {
        # CESM2 variable renaming
    },
    
    'um': {
        # UM variable renaming
    },
}


# =============================================================================
# VARIABLE SCALING CONFIGURATION
# =============================================================================

VARIABLE_SCALING = {
    'scream': {
        # Precipitation: kg m-2 s-1 → mm/h
        'pr': (3600000.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h (×1000 kg→mm, ×3600 s→h)'),
        
        # Geopotential: m2 s-2 → geopotential height in m
        # Using wildcard pattern to match all zg variables (zg850, zg500, zg_850hPa, etc.)
        'zg*': (1.0, 'm', 'Convert geopotential (m2 s-2) to geopotential height (m)'),
        # Specific Humidity: kg/kg → g/kg
        'hus*': (1000.0, 'g kg-1', 'Convert specific humidity (kg/kg) to g/kg'),
        # Relative Humidity: 1 → %
        'rh*': (100.0, '%', 'Convert relative humidity (1) to percentage (%)'),
    },
    
    'era5': {
        # ERA5 variables are typically already in standard units
        'pr': (1.0, 'mm h-1', 'Already in mm/h'),
        # Geopotential: m2 s-2 → geopotential height in m
        'zg*': (1.0 / 9.81, 'm', 'Convert geopotential (m2 s-2) to geopotential height (m)'),
        # Specific Humidity: kg/kg → g/kg
        'hus*': (1000.0, 'g kg-1', 'Convert specific humidity (kg/kg) to g/kg'),
    },
    
    'nicam': {
        # Precipitation: kg m-2 s-1 → mm/h
        'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h (×1000 kg→mm, ×3.6 s→h)'),
        
        # Geopotential: m2 s-2 → geopotential height in m
        'zg*': (1.0 / 9.81, 'm', 'Convert geopotential (m2 s-2) to geopotential height (m)'),
        # Specific Humidity: kg/kg → g/kg
        'hus*': (1000.0, 'g kg-1', 'Convert specific humidity (kg/kg) to g/kg'),
    },
    
    'icon': {
        # Precipitation: kg m-2 s-1 → mm/h
        'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h (×1000 kg→mm, ×3.6 s→h)'),
        
        # Geopotential: m2 s-2 → geopotential height in m
        'zg*': (1.0 / 9.81, 'm', 'Convert geopotential (m2 s-2) to geopotential height (m)'),
        # Specific Humidity: kg/kg → g/kg
        'hus*': (1000.0, 'g kg-1', 'Convert specific humidity (kg/kg) to g/kg'),
    },
    
    'cesm2': {
        # Precipitation: kg m-2 s-1 → mm/h
        'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h (×1000 kg→mm, ×3.6 s→h)'),
        
        # Geopotential height: already in meters
        'zg*': (1.0, 'm', 'Already in geopotential height'),
        # Specific Humidity: kg/kg → g/kg
        'hus*': (1000.0, 'g kg-1', 'Convert specific humidity (kg/kg) to g/kg'),
    },
    
    'um': {
        # Precipitation: kg m-2 s-1 → mm/h
        'pr': (3600.0, 'mm h-1', 'Convert kg m-2 s-1 to mm/h (×1000 kg→mm, ×3.6 s→h)'),
        
        # Geopotential height: already in meters
        'zg*': (1.0, 'm', 'Already in geopotential height'),
        # Specific Humidity: kg/kg → g/kg
        'hus*': (1000.0, 'g kg-1', 'Convert specific humidity (kg/kg) to g/kg'),
    },
}


def get_scaling_info(source, variable):
    """
    Get scaling information for a specific variable and source.
    
    Supports both exact matching and wildcard patterns (e.g., "zg*" matches "zg850", "zg500").
    Exact matches take precedence over pattern matches.
    
    Parameters:
    -----------
    source : str
        Source identifier (e.g., 'scream', 'era5')
    variable : str
        Variable name (e.g., 'pr', 'zg_500hPa')
    
    Returns:
    --------
    tuple or None : (scale_factor, target_units, description) or None if not configured
    """
    if source not in VARIABLE_SCALING:
        return None
    
    config = VARIABLE_SCALING[source]
    
    # First, try exact match
    if variable in config:
        return config[variable]
    
    # Second, try pattern matching (patterns ending with *)
    for pattern, scaling_info in config.items():
        if pattern.endswith('*'):
            prefix = pattern[:-1]  # Remove the * at the end
            if variable.startswith(prefix):
                return scaling_info
    
    return None


def get_all_scaled_variables(source):
    """
    Get all variables that have scaling configured for a source.
    
    Parameters:
    -----------
    source : str
        Source identifier
    
    Returns:
    --------
    list : Variable names that have scaling configured
    """
    if source in VARIABLE_SCALING:
        return list(VARIABLE_SCALING[source].keys())
    return []


def get_all_sources():
    """
    Get all configured sources.
    
    Returns:
    --------
    list : Source identifiers
    """
    return list(VARIABLE_SCALING.keys())


def rename_variables(ds, source):
    """
    Rename variables in a dataset according to source-specific naming conventions.
    
    This function standardizes variable names across different model sources.
    For example, converting model-specific names to standardized names like:
    - 'precipitation' → 'pr'
    - 'Temperature' → 'ta'
    - 'geopotential_height' → 'zg'
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing variables to rename
    source : str
        Source identifier (e.g., 'scream', 'era5', 'nicam')
    
    Returns:
    --------
    ds : xarray.Dataset
        Dataset with renamed variables
    renamed_vars : dict
        Dictionary mapping old names to new names {old_name: new_name}
    """
    if source not in VARIABLE_RENAMING:
        print(f"  No renaming configuration for source '{source}' - skipping variable renaming")
        return ds, {}
    
    rename_map = VARIABLE_RENAMING[source]
    
    if not rename_map:
        print(f"  No variables to rename for source '{source}'")
        return ds, {}
    
    # Find which variables in the dataset need renaming
    vars_to_rename = {}
    for old_name, new_name in rename_map.items():
        if old_name in ds.data_vars:
            vars_to_rename[old_name] = new_name
    
    if vars_to_rename:
        print(f"\n  Renaming variables for source: {source}")
        for old_name, new_name in vars_to_rename.items():
            print(f"    {old_name} → {new_name}")
        
        # Apply renaming
        ds = ds.rename(vars_to_rename)
    else:
        print(f"  No matching variables found to rename for source '{source}'")
    
    return ds, vars_to_rename


def get_renamed_variable(source, original_name):
    """
    Get the renamed variable name for a specific source.
    
    Parameters:
    -----------
    source : str
        Source identifier
    original_name : str
        Original variable name
    
    Returns:
    --------
    str or None : Renamed variable name, or None if no renaming configured
    """
    if source in VARIABLE_RENAMING:
        return VARIABLE_RENAMING[source].get(original_name, None)
    return None


def get_all_renamings(source):
    """
    Get all variable renamings configured for a source.
    
    Parameters:
    -----------
    source : str
        Source identifier
    
    Returns:
    --------
    dict : Dictionary mapping original names to new names
    """
    if source in VARIABLE_RENAMING:
        return VARIABLE_RENAMING[source].copy()
    return {}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

if __name__ == "__main__":
    print("Variable Scaling Configuration")
    print("=" * 80)
    
    for source in get_all_sources():
        print(f"\nSource: {source.upper()}")
        print("-" * 80)
        
        variables = get_all_scaled_variables(source)
        for var in variables:
            scale_factor, target_units, description = get_scaling_info(source, var)
            if var.endswith('*'):
                print(f"  {var:15s} : ×{scale_factor:10.6g}  →  {target_units:10s}  ({description})")
                print(f"                   [Pattern: matches all variables starting with '{var[:-1]}']")
            else:
                print(f"  {var:15s} : ×{scale_factor:10.6g}  →  {target_units:10s}  ({description})")
    
    print("\n" + "=" * 80)
    print("\nPattern Matching Examples:")
    print("-" * 80)
    
    # Demonstrate pattern matching
    test_vars = ['pr', 'zg850', 'zg500', 'zg_500hPa', 'zg_850hPa', 'ta850', 'unknown']
    test_source = 'scream'
    
    print(f"Testing with source: {test_source}")
    print(f"Variables in config: {get_all_scaled_variables(test_source)}")
    print(f"\nMatching results:")
    
    for var in test_vars:
        info = get_scaling_info(test_source, var)
        if info:
            scale_factor, target_units, description = info
            print(f"  {var:15s} → MATCH   (×{scale_factor:.6g} → {target_units})")
        else:
            print(f"  {var:15s} → NO MATCH")
    
    print("\n" + "=" * 80)
    print("\nVariable Renaming Examples:")
    print("-" * 80)
    
    print(f"Testing renaming configuration:")
    for source in get_all_sources():
        renamings = get_all_renamings(source)
        if renamings:
            print(f"\n  {source.upper()}:")
            for old_name, new_name in renamings.items():
                print(f"    {old_name} → {new_name}")
        else:
            print(f"\n  {source.upper()}: No renamings configured")
    
    print("\n" + "=" * 80)
