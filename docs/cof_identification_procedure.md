# Co-occurrence Feature (COF) Identification Procedure

**Reference script:** `scripts/make_cooccurrence_masks.py`  
**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

Co-occurrence features (COFs) are defined as spatially overlapping combinations of mesoscale convective systems (MCS), atmospheric rivers (AR), and extratropical cyclones (ETC). The identification procedure operates on gridded (HEALPix) feature masks at each time step and classifies every feature into a mutually exclusive category based on its co-occurrence relationships with other feature types. Tropical cyclones (TC) are treated as a filtering criterion rather than a COF category.

The procedure reflects an atmospheric **scale hierarchy** (ETC > AR > MCS), which is encoded in the overlap fraction thresholds: smaller-scale features (MCS) require a larger fractional overlap to be considered associated, while larger-scale features (ETC) appear in COFs even if only a small fraction of their area is involved.

The following steps are applied independently at each time step.

---

## Summary Flowchart

```
Input: MCS, AR, ETC, TC gridded feature masks (each time step)
         |
         v
[Step 1] TC Filtering
  └─ Remove MCS tracks with ≥10% TC pixel overlap
  └─ Remove cloud-type areas with ≥1% TC overlap
         |
         v
[Step 2] Binary Masks & Summation
  └─ B_MCS, B_AR, B_ETC → pairwise sums, 3-way sum
         |
         v
[Step 3] Three-Way Overlap Detection
  └─ Scan S_3way == 3 pixels
  └─ Validate each (MCS, AR, ETC) triplet with mutual overlap
  └─ Threshold: MCS ≥ 20%, AR ≥ 10%, ETC ≥ 0%
         |
         v
[Step 4] Two-Way Overlap Detection (excluding Step 3 tracks)
  └─ MCS-AR: MCS ≥ 20%, AR ≥ 0%
  └─ AR-ETC: AR ≥ 10%, ETC ≥ 1%
  └─ MCS-ETC: MCS ≥ 20%, ETC ≥ 0%
         |
         v
[Step 5] Dual Two-Way Promotion
  └─ ETC in both AR-ETC and MCS-ETC → promote to 3-way
  └─ Comprehensive expansion of connected tracks
  └─ Remove promoted tracks from 2-way lists
         |
         v
[Step 6] Mutually Exclusive Mask Generation
  └─ Isolated and overlap masks for each category
         |
         v
[Step 7] ETC Overlap Statistics (flags 0–3, track ID records)
         |
         v
Output: COF masks (zarr) + ETC statistics (CSV, Parquet)
```

---

## Key Steps at a Glance

- **Step 1 — TC Filtering:** MCS tracks that substantially overlap with tropical cyclones are removed before any co-occurrence analysis begins. This prevents TC-embedded convection from being misclassified as a non-TC COF system.

- **Step 2 — Binary Mask Construction:** Each feature type (MCS, AR, ETC) is converted to a binary presence mask, then pairwise and three-way summation masks are computed. Candidate co-occurrence regions are identified wherever the summation reaches 2 (pairwise) or 3 (three-way).

- **Step 3 — Three-Way Overlap Detection:** Locations where all three features simultaneously overlap are scanned, and each MCS–AR–ETC triplet is validated using feature-specific overlap fraction thresholds reflecting the scale hierarchy. This step identifies the highest-order COF category first, before any two-way analysis.

- **Step 4 — Two-Way Overlap Detection:** Among tracks not already classified as three-way, each feature pair (MCS-AR, AR-ETC, MCS-ETC) is independently checked for spatial co-occurrence using pairwise overlap fraction thresholds.

- **Step 5 — Dual Two-Way Promotion:** If an ETC independently overlaps both an AR and an MCS in separate two-way pairs, the three features are promoted to a three-way COF even if they never all converge at the same pixel. The promotion propagates transitively to all connected tracks to ensure physically coherent grouping.

- **Step 6 — Mutually Exclusive Mask Generation:** Each feature is assigned to exactly one category—isolated, one of three two-way pair types, or the three-way group—based on its highest-order co-occurrence relationship. Separate output masks are produced for each category to support feature-type-specific analysis.

- **Step 7 — ETC Overlap Statistics:** Each ETC track receives an overlap flag (0–3) indicating whether it co-occurs with MCS, AR, both, or neither. Associated MCS and AR track IDs are recorded and saved for subsequent statistical analysis.

---

## Step 1 — TC Filtering

Prior to any co-occurrence analysis, MCS features that substantially overlap with tropical cyclones are removed. This prevents TC-embedded convection from being misclassified as non-TC COF systems.

- A binary overlap mask is constructed by summing the MCS binary mask and the TC binary mask.
- For each MCS track, the fraction of its total pixels that overlap with any TC pixel is computed.
- MCS tracks with an overlap fraction **≥ 10%** are removed from subsequent processing.
- The same procedure is applied to the cloud-type mask using a more conservative threshold of **≥ 1%**, to exclude cloud areas even loosely associated with TC circulation.

**Result:** A TC-filtered MCS mask and a TC-filtered cloud-type mask are used for all subsequent steps.

---

## Step 2 — Binary Mask Construction and Summation

For each of the three feature types—MCS (TC-filtered), AR, and ETC—a binary presence mask is constructed:

$$B_f(x) = \begin{cases} 1 & \text{if feature } f \text{ is present at pixel } x \\ 0 & \text{otherwise} \end{cases}$$

Three **pairwise summation masks** are formed:
$$S_{\text{MCS-AR}} = B_{\text{MCS}} + B_{\text{AR}}, \quad S_{\text{AR-ETC}} = B_{\text{AR}} + B_{\text{ETC}}, \quad S_{\text{MCS-ETC}} = B_{\text{MCS}} + B_{\text{ETC}}$$

A **three-way summation mask** is also formed:
$$S_{\text{3-way}} = B_{\text{MCS}} + B_{\text{AR}} + B_{\text{ETC}}$$

Pixels with a summation value of 2 (pairwise) or 3 (three-way) identify candidate overlap regions.

---

## Step 3 — Three-Way Overlap Detection

Three-way COFs involve the simultaneous spatial co-occurrence of all three feature types (MCS, AR, ETC). The identification requires **mutual spatial verification** to ensure all three features genuinely co-locate, not merely that they each appear somewhere in the domain.

1. Candidate **three-way pixels** are identified where $S_{\text{3-way}} = 3$.
2. For each combination of candidate tracks $(i_{\text{MCS}},\; j_{\text{AR}},\; k_{\text{ETC}})$ present in the three-way region:
   - The **mutual overlap pixel count** $N_{ijk}$ is computed as the number of pixels where all three tracks simultaneously overlap within the three-way region.
   - **Overlap fractions** are calculated for each track:
     $$f_{\text{MCS}} = \frac{N_{ijk}}{N_{\text{MCS},i}}, \quad f_{\text{AR}} = \frac{N_{ijk}}{N_{\text{AR},j}}, \quad f_{\text{ETC}} = \frac{N_{ijk}}{N_{\text{ETC},k}}$$
     where $N_{\cdot,i}$ is the total pixel count of track $i$.
3. The triplet $(i, j, k)$ is accepted as a **validated three-way COF** if all three tracks meet their respective overlap fraction thresholds:

| Feature | Three-way threshold |
|---------|-------------------|
| MCS     | ≥ 20%             |
| AR      | ≥ 10%             |
| ETC     | ≥ 0%              |

The asymmetric thresholds follow the scale hierarchy: since an ETC typically spans a much larger area than an MCS, even a small fraction of ETC area overlapping with an MCS is considered physically meaningful.

---

## Step 4 — Two-Way Overlap Detection

Two-way COFs are defined for each feature pair (MCS-AR, AR-ETC, MCS-ETC) among tracks **not already classified as three-way**.

For each feature pair, overlap detection follows three sub-steps:

1. **Overlap region identification:** Pixels where the pairwise summation mask $S \geq 2$ define the candidate overlap region, excluding NaN values.
2. **Overlap fraction calculation:** For each track present in the overlap region:
   $$f_{\text{track}} = \frac{N_{\text{overlap}}}{N_{\text{total}}}$$
3. **Pair validation:** A track pair $(i, j)$ is accepted only if both tracks meet their overlap fraction thresholds AND they are spatially co-located (i.e., their qualifying pixels mutually intersect).

The thresholds applied for each pair type are:

| Pair type | Feature 1 threshold | Feature 2 threshold |
|-----------|--------------------|--------------------|
| MCS-AR    | MCS ≥ 20%          | AR ≥ 0%            |
| AR-ETC    | AR ≥ 10%           | ETC ≥ 1%           |
| MCS-ETC   | MCS ≥ 20%          | ETC ≥ 0%           |

Tracks validated as three-way (Step 3) are excluded before forming two-way pairs.

---

## Step 5 — Dual Two-Way Overlap Promotion to Three-Way

An ETC track may form two-way pairs with both an AR and an MCS independently, without those three features mutually co-locating at a single set of pixels. Such cases represent an implicit three-way association (the ETC mediates a link between the MCS and the AR) and are **promoted to the three-way category**.

The promotion procedure works as follows:

1. **Identify dual-overlap ETC tracks:** ETC tracks appearing in both the AR-ETC two-way pair list and the MCS-ETC two-way pair list are flagged.
2. **Form promoted triplets:** For each flagged ETC $k$, all combinations with its AR partners $j$ and MCS partners $i$ are formed as promoted three-way COF triplets $(i, j, k)$.
3. **Comprehensive expansion:** To ensure group coherence, the promotion is propagated transitively:
   - All MCS tracks that have a two-way overlap with any promoted AR track are added to the three-way category.
   - All AR tracks that overlap with any promoted MCS or ETC track are added.
   - All ETC tracks that overlap with any promoted AR or MCS track are added.
4. **Update two-way pair lists:** All tracks included in the expanded three-way group are removed from the two-way pair lists.

This step ensures that all COF features within a physically connected system are consistently classified at the highest-order overlap level.

---

## Step 6 — Mutually Exclusive Mask Generation

After completion of Steps 1–5, each feature is assigned to exactly one mutually exclusive category based on the highest-order co-occurrence relationship it participates in:

| Category | Description |
|----------|-------------|
| **Isolated MCS** | MCS not overlapping AR or ETC (after TC removal) |
| **Isolated AR** | AR not overlapping MCS or ETC |
| **Isolated ETC** | ETC not overlapping MCS or AR |
| **MCS-AR (2-way)** | MCS and AR with direct spatial overlap |
| **AR-ETC (2-way)** | AR and ETC with direct spatial overlap |
| **MCS-ETC (2-way)** | MCS and ETC with direct spatial overlap |
| **MCS-AR-ETC (3-way)** | All three features mutually co-located |

Separate output masks are produced for both perspectives of each pair (e.g., the MCS component of an MCS-AR pair and the AR component of the same pair) to facilitate feature-type-specific analysis.

---

## Step 7 — ETC Overlap Statistics

For each ETC track at each time step, an **overlap flag** is recorded to characterize its co-occurrence state:

| Flag | Description |
|------|-------------|
| 0 | Isolated (no overlap with MCS or AR) |
| 1 | Overlaps with MCS only |
| 2 | Overlaps with AR only |
| 3 | Overlaps with both MCS and AR |

For ETC tracks in three-way COFs, the associated MCS and AR track IDs are recorded using the transitive relationship: MCS tracks that overlap with an AR that overlaps with the ETC are also considered associated with that ETC. This information is saved as tabular data (CSV and Parquet) for subsequent statistical analysis.

---

## Overlap Fraction Thresholds Summary

| Analysis | Feature | Threshold | Rationale |
|----------|---------|-----------|-----------|
| TC filtering | MCS | ≥ 10% | Remove TC-embedded convection |
| TC filtering | Cloud types | ≥ 1% | Conservative removal of TC cloud areas |
| 3-way COF | MCS | ≥ 20% | High threshold for smallest feature |
| 3-way COF | AR | ≥ 10% | Moderate threshold for medium feature |
| 3-way COF | ETC | ≥ 0% | Any presence accepted for largest feature |
| 2-way MCS-AR | MCS | ≥ 20% | Consistent with 3-way |
| 2-way MCS-AR | AR | ≥ 0% | AR need only touch MCS |
| 2-way AR-ETC | AR | ≥ 10% | Consistent with 3-way |
| 2-way AR-ETC | ETC | ≥ 1% | Minimal ETC area overlap required |
| 2-way MCS-ETC | MCS | ≥ 20% | Consistent with 3-way |
| 2-way MCS-ETC | ETC | ≥ 0% | Any ETC presence accepted |
