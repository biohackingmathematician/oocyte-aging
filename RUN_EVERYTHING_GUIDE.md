# How to Run Everything Again

## Quick Start - Run All Visualizations

### Option 1: Run Everything with One Command (Easiest)

```bash
bash run_all_visualizations.sh
```

Or if that doesn't work:

```bash
./run_all_visualizations.sh
```

This will:
1. Generate forum visualizations from CSV data
2. Generate improved/fixed versions
3. Fix all existing visualizations (risk, gplvm, complete summary)

---

### Option 2: Run Scripts Individually

Run these scripts in order:

```bash
# 1. Generate forum visualizations (stage overview & boxplot)
python3 run_visualizations_from_csv.py

# 2. Generate improved forum visualizations with better spacing
python3 review_and_fix_visualizations.py

# 3. Fix existing visualizations (risk stratification, gplvm, complete summary)
python3 fix_existing_visualizations.py
```

---

### Option 3: Run Individual Scripts (If Needed)

#### Forum Visualizations Only:
```bash
python3 run_visualizations_from_csv.py
```
**Generates:**
- `forum_stage_overview.png`
- `forum_cellular_age_z_boxplot.png`

#### Improved Forum Visualizations:
```bash
python3 review_and_fix_visualizations.py
```
**Generates:**
- `forum_stage_overview_fixed.png` (improved)
- `forum_comparison_boxplot_fixed.png` (improved)

#### Fix All Existing Visualizations:
```bash
python3 fix_existing_visualizations.py
```
**Generates:**
- `risk_stratification_fixed.png`
- `gplvm_trajectory_analysis_fixed.png`
- `complete_results_summary_fixed.png`

---

## What Each Script Does

### `run_visualizations_from_csv.py`
- Loads data from CSV files (no h5ad needed)
- Creates stage overview and boxplot visualizations
- Works with basic Python packages only

### `review_and_fix_visualizations.py`
- Creates improved versions of forum visualizations
- Better spacing and layout
- No blank spaces
- Publication-ready quality

### `fix_existing_visualizations.py`
- Fixes risk stratification visualization (3-panel)
- Fixes GPLVM trajectory analysis (3-panel)
- Fixes complete results summary (multi-panel GridSpec)
- Removes all blank spaces
- Adds placeholder text when data missing

---

## Generated Files

After running all scripts, you'll have:

### Forum Visualizations
- `forum_stage_overview.png` (original)
- `forum_stage_overview_fixed.png` (✅ **use this one**)
- `forum_cellular_age_z_boxplot.png` (original)
- `forum_comparison_boxplot_fixed.png` (✅ **use this one**)

### Main Analysis Visualizations  
- `risk_stratification_fixed.png` (✅ **use this one**)
- `gplvm_trajectory_analysis_fixed.png` (✅ **use this one**)
- `complete_results_summary_fixed.png` (✅ **use this one**)

**Recommendation:** Use the `*_fixed.png` versions as they have:
- No blank spaces
- Better spacing
- Consistent styling
- Ready for publication

---

## Requirements

- Python 3 (works with 3.14, but 3.11.9 is recommended)
- Required packages:
  - pandas
  - numpy
  - matplotlib
  - scipy

These are usually already installed. If not:
```bash
pip3 install pandas numpy matplotlib scipy
```

---

## Troubleshooting

### If you get "command not found" errors:
```bash
# Use python3 explicitly
python3 run_visualizations_from_csv.py
```

### If you get missing package errors:
```bash
pip3 install pandas numpy matplotlib scipy
```

### If CSV files are missing:
Make sure these files exist:
- `sample_metadata_with_age.csv`
- `clinical_decision_framework_final.csv`

---

## Quick Command Reference

```bash
# Run everything
bash run_all_visualizations.sh

# Or step by step:
python3 run_visualizations_from_csv.py
python3 review_and_fix_visualizations.py  
python3 fix_existing_visualizations.py
```

All scripts run locally without needing Jupyter or h5ad files! ✅

