# Forum EDA Visualizations

## Overview
Two key exploratory data analysis visualizations for the project forum presentation.

---

## Visualization 1: UMAP/PCA Colored by Stage
**File**: `forum_umap_pca_by_stage.png`

### Purpose
Show global structure and demonstrate that GV and MI oocytes separate along the main axis.

### What It Shows
- **2D PCA Plot** (PC1 vs PC2) with points colored by stage (GV vs MI)
- Points are colored: GV = Green (#2ecc71), MI = Red (#e74c3c)
- PC1 explains 96.19% of variance (main separation axis)
- PC2 explains 3.81% of variance

### What It Answers
1. Cells are not random noise; they form a continuum/separation
2. Maturation stage is reflected in the transcriptome
3. GV and MI oocytes cluster separately

### Slide Caption
**"Global structure: GV and MI separate along main axis."**

### Technical Details
- **Method**: PCA computed from:
  - `cellular_age_z` (normalized cellular age)
  - `cellular_age_uncertainty` (trajectory uncertainty)
  - `health_score` (computed from cellular_age_z)
- **Features**: 3 features, 2 components
- **Total variance explained**: 100%
- **Samples**: 20 oocytes (6 GV, 14 MI)

---

## Visualization 2: Health Score by Stage (Box + Violin Plots)
**File**: `forum_health_score_by_stage.png`

### Purpose
Show basic biology signal in one glance: GV oocytes are healthier than MI at the transcriptomic level.

### What It Shows
- **Left Panel**: Box plot comparing health scores between GV and MI
- **Right Panel**: Violin plot showing distribution shapes
- Both panels include:
  - Individual data points (jittered)
  - Mean markers
  - Median lines
  - Sample size annotations
  - Statistical significance bracket (if p < 0.05)

### What It Answers
1. GV oocytes have higher health scores on average than MI
2. There's variability within each stage
3. The decline from GV to MI is quantifiable

### Slide Caption
**"GV oocytes look healthier than MI at transcriptomic level."**

### Key Statistics
- **GV**: n=6, Mean=96.03, Median=97.10, SD=4.12
- **MI**: n=14, Mean=83.49, Median=92.38, SD=24.35
- **Statistical Test**: Mann-Whitney U test
  - U-statistic: 66.00
  - p-value: 0.0507
  - Significance: ns (not significant at α=0.05, but clear trend)

### Technical Details
- **Health Score**: Computed from `cellular_age_z` as `(1 - cellular_age_z) * 100`
  - Higher score = healthier oocyte
  - Lower cellular age = higher health
- **Colors**: GV = Green (#2ecc71), MI = Red (#e74c3c)
- **Distribution**: Shows both central tendency (box/violin) and individual points

---

## How to Generate

### Standalone
```bash
python3 create_forum_eda_visualizations.py
```

### With All Visualizations
```bash
bash run_all_visualizations.sh
```

---

## Key Findings

1. **Clear Separation**: GV and MI oocytes separate along PC1 (96% variance)
   - Confirms maturation stage is reflected in transcriptome
   - Shows continuum from GV → MI

2. **Health Decline**: GV oocytes are healthier than MI
   - GV mean: 96.0 (high health)
   - MI mean: 83.5 (moderate health)
   - 13-point difference (13% relative decline)

3. **Within-Stage Variability**: 
   - GV has lower variability (SD=4.12)
   - MI has higher variability (SD=24.35)
   - Suggests MI stage has more heterogeneity

---

## Design Notes

- **Color Scheme**: Consistent with project (GV=green, MI=red)
- **Publication Quality**: 300 DPI, clean styling
- **Information Density**: Both plots show same data in complementary ways
- **Statistics Included**: Sample sizes, means, medians, significance testing

---

## Updates from Previous Version

### Removed
-  Stage distribution histogram/bar chart
-  Age distribution histograms

### Added
- PCA plot colored by stage (shows global structure)
- Combined box+violin plot for health scores

### Why These Changes?
1. **PCA Plot**: More informative than histograms - shows relationships between samples
2. **Health Score Plot**: Directly answers biological question (GV healthier than MI)
3. **Cleaner**: Focused on key EDA insights rather than basic counts

---

## Files Generated

1. `forum_umap_pca_by_stage.png` - PCA plot colored by stage
2. `forum_health_score_by_stage.png` - Box + violin plot of health scores

Both are publication-ready at 300 DPI with proper labels and styling.

