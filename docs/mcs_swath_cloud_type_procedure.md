# MCS Swath Mask and Cloud Type Classification Procedure

**Reference script:** `scripts/make_mcs_swath_masks.py`  
**Author:** Zhe Feng | zhe.feng@pnnl.gov  
**Purpose:** Methodology notes for writing the JGR: Atmosphere paper

---

## Overview

This procedure generates temporally aggregated MCS swath masks and cloud type classifications from high-resolution global storm-resolving model (GSRM) output on a HEALPix grid. Starting from hourly tracked MCS pixel masks, brightness temperature (Tb), and precipitation (pr), each sub-daily aggregation window (e.g., 6 hours) is collapsed into a single time step describing:

1. **MCS swath mask** — the spatial footprint of each MCS track over the aggregation window, with overlapping tracks resolved by coverage priority.
2. **Cloud type classification** — a mutually exclusive, priority-based classification of non-MCS cloud areas outside the MCS swath.
3. **Frequency-weighted mean precipitation** — mean precipitation attributed to each cloud type, weighted by how frequently that type occurred at each grid cell during the window.

The procedure operates independently on each aggregation window and produces output on the same HEALPix grid as the input.

---

## Summary Flowchart

```
Input: hourly MCS pixel masks + Tb + Precipitation (HEALPix grid)
         |
         v
[Step 1] Temporal Aggregation Setup
  └─ Group hourly time steps into N-hour windows (default: 6h)
  └─ Align output times to standard hour boundaries (00, 06, 12, 18 UTC)
         |
         v
[Step 2] MCS Swath Construction
  └─ Per-track: union of pixels over all timesteps in window → swath
  └─ Per-track: count how many timesteps cover each pixel → coverage
  └─ Combine tracks: overlapping pixels assigned to track with highest coverage
         |
         v
[Step 3] Latitude-Dependent Tb Threshold
  └─ Tropics (|lat| ≤ 30°): 250 K
  └─ Mid-latitudes (30° < |lat| ≤ 60°): linearly decrease 250→230 K
  └─ High latitudes (|lat| > 60°): 230 K
         |
         v
[Step 4] Per-Timestep Cloud Type Classification
  └─ Applied only outside MCS areas (mcs_mask == 0)
  └─ Type 1: Deep convective   (Tb < Tb_thresh  &  pr ≥ 0.5 mm/h)
  └─ Type 2: Stratiform        (Tb < Tb_thresh  &  pr < 0.5 mm/h)
  └─ Type 3: Non-deep conv.    (Tb ≥ Tb_thresh  &  pr ≥ 0.5 mm/h)
  └─ Type 4: Drizzle           (Tb ≥ Tb_thresh  &  pr < 0.5 mm/h)
         |
         v
[Step 5] Priority-Based Cloud Type Aggregation
  └─ Aggregate over all timesteps in window using priority 1 > 2 > 3 > 4
  └─ Apply MCS priority: set cloud type = 0 where MCS swath exists
         |
         v
[Step 6] Frequency-Weighted Mean Precipitation by Cloud Type
  └─ For each type: conditional mean pr × frequency of that type
  └─ MCS swath pixels set to 0 for all cloud type precipitation
         |
         v
Output: zarr file with mcs_mask, cloud_types, dc_pr, st_pr, nd_pr, dz_pr
```

---

## Key Steps at a Glance

- **Step 1 — Temporal Aggregation Setup:** Hourly time steps are grouped into fixed-duration aggregation windows (default: 6 hours) aligned to standard UTC boundaries (00, 06, 12, 18). Each output time step thus represents a compact swath of activity across multiple model hours.

- **Step 2 — MCS Swath Construction:** For each MCS track, the union of all pixels it occupies during the window is computed to form a 2D swath footprint, along with a pixel-level count of how many timesteps it was present. Where multiple tracks overlap, the track with the highest pixel coverage count takes priority.

- **Step 3 — Latitude-Dependent Tb Threshold:** A spatially varying brightness temperature threshold is constructed that decreases poleward, reflecting the lower cloud top temperatures needed to detect deep convection at higher latitudes. This ensures physically consistent cloud detection across the globe.

- **Step 4 — Per-Timestep Cloud Type Classification:** At each hourly time step, every non-MCS grid cell is classified into one of four mutually exclusive cloud types using Tb relative to the threshold and whether precipitation exceeds 0.5 mm/h. MCS pixels are excluded from this classification entirely.

- **Step 5 — Priority-Based Aggregation:** Cloud types across all timesteps in the window are collapsed to a single value per cell using a strict priority order (deep convective > stratiform > non-deep convective > drizzle). The final cloud type map is then masked to zero wherever the MCS swath is present, enforcing mutual exclusivity between MCS and cloud type categories.

- **Step 6 — Frequency-Weighted Mean Precipitation:** For each cloud type, mean precipitation is computed as the product of the conditional mean (average precipitation when the cell exhibited that type) and the type frequency (fraction of timesteps). The four contributions sum to the simple time-mean precipitation over the window, and are set to zero within MCS swath pixels.

---

## Step 1 — Temporal Aggregation Setup

Hourly input time steps are grouped into aggregation windows of configurable length (default: **6 hours**). Output times are aligned to standard UTC hour boundaries so that all model data sources produce comparable output time coordinates:

- **Boundary hours:** 00, 06, 12, 18 UTC
- **Alignment method:** floor each input timestamp to the nearest aggregation-window boundary

Each output time step aggregates all input time steps that fall within the same window. The number of input steps per window is typically equal to the aggregation window length (e.g., 6 steps for 6-hour windows from hourly data), but may differ near the start and end of the dataset.

---

## Step 2 — MCS Swath Construction

The MCS swath for a given aggregation window is derived in two sub-steps.

### 2a — Per-Track Swath and Coverage

For each unique MCS track $i$ present during the aggregation window, two 2D arrays are computed from the 3D track number array $T(t, x)$:

**Swath** — the spatial union of all pixels occupied by track $i$ at any time step in the window:
$$\text{swath}_i(x) = \begin{cases} i & \text{if } \exists\, t : T(t, x) = i \\ 0 & \text{otherwise} \end{cases}$$

**Coverage count** — the number of time steps during which track $i$ occupies pixel $x$:
$$\text{cov}_i(x) = \sum_t \mathbf{1}[T(t, x) = i]$$

### 2b — Priority-Based Combination

All per-track swaths are combined into a single 2D MCS swath mask. Where multiple tracks overlap at pixel $x$, the track with the highest coverage count is assigned:
$$\text{MCS\_swath}(x) = \underset{i\,:\;\text{swath}_i(x) > 0}{\arg\max}\;\text{cov}_i(x)$$

If no track is present at pixel $x$, the combined swath value is 0.

---

## Step 3 — Latitude-Dependent Brightness Temperature Threshold

A spatially varying Tb threshold $T_{\text{thresh}}(\phi)$ is constructed from absolute latitude $|\phi|$:

$$T_{\text{thresh}}(\phi) = \begin{cases} 250\;\text{K} & |\phi| \leq 30° \\ 250 - 20 \times \dfrac{|\phi| - 30}{30}\;\text{K} & 30° < |\phi| \leq 60° \\ 230\;\text{K} & |\phi| > 60° \end{cases}$$

The threshold decreases linearly from 250 K in the tropics to 230 K at 60° latitude, then remains constant poleward. This accounts for the lower tropopause at higher latitudes, which reduces the minimum attainable cloud-top temperatures for deep convection.

---

## Step 4 — Per-Timestep Cloud Type Classification

At each time step $t$ within the aggregation window, grid cells are classified into four mutually exclusive cloud types based on two criteria:

| Cloud type | Value | Condition |
|-----------|-------|-----------|
| Deep convective | 1 | $T_b < T_{\text{thresh}}$ **and** $\text{pr} \geq 0.5\;\text{mm\,h}^{-1}$ |
| Stratiform | 2 | $T_b < T_{\text{thresh}}$ **and** $\text{pr} < 0.5\;\text{mm\,h}^{-1}$ |
| Non-deep convective | 3 | $T_b \geq T_{\text{thresh}}$ **and** $\text{pr} \geq 0.5\;\text{mm\,h}^{-1}$ |
| Drizzle | 4 | $T_b \geq T_{\text{thresh}}$ **and** $\text{pr} < 0.5\;\text{mm\,h}^{-1}$ |
| Unclassified | 0 | Cell is inside an MCS region ($\text{mcs\_mask} > 0$) |

Classification is performed only at cells where the MCS mask equals zero. All MCS pixels receive cloud type 0 at this stage, ensuring that MCS precipitation and cloud-top properties are not double-counted.

---

## Step 5 — Priority-Based Cloud Type Aggregation

### 5a — Temporal Priority Aggregation

After per-timestep classification, cloud types across all $N_t$ time steps in the window are aggregated into a single value per cell using a **strict priority order** (1 > 2 > 3 > 4):

$$\text{cloud\_type\_agg}(x) = \min\bigl\{c \in \{1,2,3,4\} : \exists\,t\; \text{ s.t. } \text{cloud\_type}(t, x) = c\bigr\}$$

In practice: if a cell ever exhibits deep convection (type 1) during the window, it is classified as type 1 regardless of other types. Only if no type 1 occurs does type 2 take effect, and so on. Cells with no occurrence of any type remain 0.

### 5b — MCS Priority Masking

After temporal aggregation, an additional MCS priority mask is applied:

$$\text{cloud\_type\_final}(x) = \begin{cases} 0 & \text{if } \text{MCS\_swath}(x) > 0 \\ \text{cloud\_type\_agg}(x) & \text{otherwise} \end{cases}$$

This enforces the classification hierarchy: **MCS > cloud types > unclassified**, so any cell within the MCS swath footprint is excluded from cloud type statistics.

---

## Step 6 — Frequency-Weighted Mean Precipitation by Cloud Type

For each cloud type $c \in \{1, 2, 3, 4\}$, the **frequency-weighted mean precipitation** is defined as:

$$\overline{P}_c(x) = \bar{P}_{c \mid c}(x) \cdot f_c(x)$$

where:
- $\bar{P}_{c \mid c}(x) = \dfrac{1}{N_{c}(x)} \displaystyle\sum_{t:\,\text{type}(t,x)=c} P(t,x)$ is the **conditional mean precipitation** (average pr when cell $x$ has type $c$)
- $f_c(x) = \dfrac{N_c(x)}{N_t}$ is the **frequency** (fraction of time steps cell $x$ had type $c$)
- $N_c(x)$ is the number of time steps with type $c$ at cell $x$
- $N_t$ is the total number of time steps in the window

This decomposition satisfies the additive property:

$$\sum_{c=1}^{4} \overline{P}_c(x) = \bar{P}(x) \quad \text{(simple time-mean precipitation)}$$

After computing the four components, all are set to zero within MCS swath pixels, so that precipitation within MCS footprints is counted separately via the MCS swath mask.

---

## Output Variables

All output variables are saved on the HEALPix grid (1D cell dimension) at the aggregated 6-hourly time coordinate:

| Variable | Units | Description |
|----------|-------|-------------|
| `mcs_mask` | 1 | MCS swath mask; pixel value = MCS track number (0 = no MCS) |
| `cloud_types` | 1 | Priority-based aggregated cloud type (0–4); 0 inside MCS swath |
| `dc_pr` | mm h⁻¹ | Frequency-weighted mean precipitation for deep convective clouds |
| `st_pr` | mm h⁻¹ | Frequency-weighted mean precipitation for stratiform clouds |
| `nd_pr` | mm h⁻¹ | Frequency-weighted mean precipitation for non-deep convective clouds |
| `dz_pr` | mm h⁻¹ | Frequency-weighted mean precipitation for drizzle |

`dc_pr + st_pr + nd_pr + dz_pr` equals the time-mean precipitation at all non-MCS cells.

---

## Threshold and Parameter Summary

| Parameter | Value | Description |
|-----------|-------|-------------|
| Aggregation window | 6 h (default) | Number of hourly steps collapsed into one output swath |
| Output time alignment | 00, 06, 12, 18 UTC | Standard hour boundaries for aligned output times |
| Tb threshold (tropics) | 250 K | |lat| ≤ 30° |
| Tb threshold (mid-latitudes) | 250–230 K (linear) | 30° < |lat| ≤ 60° |
| Tb threshold (high latitudes) | 230 K | |lat| > 60° |
| Precipitation threshold | 0.5 mm h⁻¹ | Separates precipitating from non-precipitating cloud types |
| Cloud type priority order | 1 > 2 > 3 > 4 | Deep convective takes highest priority in temporal aggregation |
| MCS coverage tie-breaking | Max coverage count | Track with most timestep overlap takes priority at contested pixels |
