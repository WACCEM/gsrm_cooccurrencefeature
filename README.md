# Co-occurring Feature Analysis for Global Storm-Resolving Models

Identification and analysis of co-occurring atmospheric features and their associated precipitation statistics for global storm-resolving models (GSRMs). 

The atmospheric features include:
- Tropical cyclone
- Extra tropical cyclone
- Atmospheric river
- Mesoscale convective system

**Author**: Zhe Feng (zhe.feng@pnnl.gov)

## Analysis Steps

### 1. Environment Setup
```bash
# On NERSC Perlmutter
source activate /global/common/software/m1867/python/hackathon
cd /path/to/gsrm_cooccurrencefeature
```

### 2. Process Data

**Combine 4 feature masks (for GSRMs):**
```bash
cd /scripts
python combine_tracking_masks.py
```

**Remap ERA5 lat/lon masks to HEALPix:**
```bash
cd /scripts
python remap_era5_masks_healpix.py
```

**Combine 4 feature masks (for observations):**
```bash
cd /scripts
python combine_era5_imerg_tracking_masks.py
```

**Make co-occurrence feature masks:**
```bash
cd /scripts
python make_cooccurrence_masks.py --workers 64 --source SOURCE
```

### 3. Compute statistics

**Calculate monthly total precipitation for individual features:**
```bash
cd /scripts
python calc_monthly_rainmap_by_featuretypes.py --config CONFIG.yaml --source SOURCE
```

### 4. Visualization

**Make feature mask animations:**
```bash
cd /scripts
python make_animation.py
```