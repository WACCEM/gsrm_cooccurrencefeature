#!/usr/bin/env python
"""
Test script for variable scaling configuration with pattern matching.

This script tests that:
1. Pattern matching works correctly (zg* matches zg850, zg500, etc.)
2. Exact matches take precedence over patterns
3. No unintended variables are matched
4. All expected variables are matched

Usage:
    python test_variable_scaling.py
"""

import sys
from variable_scaling_config import get_scaling_info, VARIABLE_SCALING, get_all_sources


def test_pattern_matching():
    """Test pattern matching functionality."""
    print("="*80)
    print("Testing Pattern Matching")
    print("="*80)
    
    test_cases = {
        'scream': {
            'should_match': ['pr', 'zg850', 'zg500', 'zg_500hPa', 'zg_850hPa', 'zg250'],
            'should_not_match': ['ta850', 'ua850', 'va500', 'hus850', 'unknown'],
        },
        'era5': {
            'should_match': ['pr', 'zg850', 'zg500', 'zg_500hPa'],
            'should_not_match': ['pr_ar', 'ta850', 'unknown'],
        },
    }
    
    all_passed = True
    
    for source, tests in test_cases.items():
        print(f"\n{source.upper()}:")
        print("-"*80)
        
        # Test variables that should match
        for var in tests['should_match']:
            info = get_scaling_info(source, var)
            if info:
                scale_factor, units, desc = info
                print(f"  ✓ {var:20s} matched (×{scale_factor:.6g} → {units})")
            else:
                print(f"  ✗ {var:20s} FAILED - should have matched!")
                all_passed = False
        
        # Test variables that should NOT match
        for var in tests['should_not_match']:
            info = get_scaling_info(source, var)
            if info:
                scale_factor, units, desc = info
                print(f"  ✗ {var:20s} FAILED - should NOT have matched (got ×{scale_factor:.6g})")
                all_passed = False
            else:
                print(f"  ✓ {var:20s} correctly NOT matched")
    
    return all_passed


def test_exact_match_precedence():
    """Test that exact matches take precedence over patterns."""
    print("\n" + "="*80)
    print("Testing Exact Match Precedence")
    print("="*80)
    
    # Create a test config with both exact and pattern matches
    test_config = {
        'test_source': {
            'pr': (1000.0, 'mm/day', 'Exact match for pr'),
            'pr*': (1.0, 'mm/h', 'Pattern match for pr*'),
            'zg850': (10.0, 'dam', 'Exact match for zg850'),
            'zg*': (1.0 / 9.81, 'm', 'Pattern match for zg*'),
        }
    }
    
    # Temporarily add to global config
    original_config = VARIABLE_SCALING.copy()
    VARIABLE_SCALING['test_source'] = test_config['test_source']
    
    all_passed = True
    
    # Test exact match for pr
    info = get_scaling_info('test_source', 'pr')
    if info and info[0] == 1000.0:
        print(f"  ✓ 'pr' used exact match (×{info[0]:.6g})")
    else:
        print(f"  ✗ 'pr' FAILED - should use exact match")
        all_passed = False
    
    # Test pattern match for pr_something
    info = get_scaling_info('test_source', 'pr_ar')
    if info and info[0] == 1.0:
        print(f"  ✓ 'pr_ar' used pattern match (×{info[0]:.6g})")
    else:
        print(f"  ✗ 'pr_ar' FAILED - should use pattern match")
        all_passed = False
    
    # Test exact match for zg850
    info = get_scaling_info('test_source', 'zg850')
    if info and info[0] == 10.0:
        print(f"  ✓ 'zg850' used exact match (×{info[0]:.6g})")
    else:
        print(f"  ✗ 'zg850' FAILED - should use exact match")
        all_passed = False
    
    # Test pattern match for zg500
    info = get_scaling_info('test_source', 'zg500')
    expected_factor = 1.0 / 9.81
    if info and abs(info[0] - expected_factor) < 0.0001:
        print(f"  ✓ 'zg500' used pattern match (×{info[0]:.6g})")
    else:
        print(f"  ✗ 'zg500' FAILED - should use pattern match")
        all_passed = False
    
    # Restore original config
    del VARIABLE_SCALING['test_source']
    
    return all_passed


def test_all_sources():
    """Test that all sources have valid configuration."""
    print("\n" + "="*80)
    print("Testing All Sources Configuration")
    print("="*80)
    
    all_passed = True
    
    for source in get_all_sources():
        print(f"\n{source.upper()}:")
        
        config = VARIABLE_SCALING[source]
        
        # Check that config has required variables
        has_pr = False
        has_zg = False
        
        for pattern in config.keys():
            if pattern == 'pr' or pattern.startswith('pr'):
                has_pr = True
            if pattern.startswith('zg'):
                has_zg = True
        
        if has_pr:
            print(f"  ✓ Has precipitation scaling")
        else:
            print(f"  ✗ Missing precipitation scaling")
            all_passed = False
        
        if has_zg:
            print(f"  ✓ Has geopotential scaling")
        else:
            print(f"  ✗ Missing geopotential scaling")
            all_passed = False
    
    return all_passed


def main():
    """Run all tests."""
    print("\nVariable Scaling Configuration Tests")
    print("="*80)
    
    results = []
    
    # Run tests
    results.append(("Pattern Matching", test_pattern_matching()))
    results.append(("Exact Match Precedence", test_exact_match_precedence()))
    results.append(("All Sources Configuration", test_all_sources()))
    
    # Print summary
    print("\n" + "="*80)
    print("Test Summary")
    print("="*80)
    
    all_passed = True
    for test_name, passed in results:
        status = "PASSED" if passed else "FAILED"
        symbol = "✓" if passed else "✗"
        print(f"  {symbol} {test_name:30s}: {status}")
        if not passed:
            all_passed = False
    
    print("="*80)
    
    if all_passed:
        print("\n✓ All tests PASSED!")
        return 0
    else:
        print("\n✗ Some tests FAILED!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
