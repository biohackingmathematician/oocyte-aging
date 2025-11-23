# Quick Run Guide - Generate All Visualizations

## ✅ Easiest Way: Run Everything at Once

```bash
bash run_all_visualizations.sh
```

**That's it!** This will generate all visualization files.

---

## What This Does

The script runs three Python scripts in sequence:

1. **`run_visualizations_from_csv.py`**
   - Creates forum visualizations from CSV data
   - Generates: `forum_stage_overview.png` and boxplot

2. **`review_and_fix_visualizations.py`**
   - Creates improved/fixed versions with better spacing
   - Generates: `forum_stage_overview_fixed.png` and `forum_comparison_boxplot_fixed.png`

3. **`fix_existing_visualizations.py`**
   - Fixes all existing visualizations (risk, gplvm, complete summary)
   - Generates: `risk_stratification_fixed.png`, `gplvm_trajectory_analysis_fixed.png`, `complete_results_summary_fixed.png`

---

## Output Files

After running, you'll have **9 visualization files**:

### Forum Visualizations
- `forum_stage_overview.png` (original)
- `forum_stage_overview_fixed.png` ✅ **USE THIS**
- `forum_cellular_age_z_boxplot.png` (original)  
- `forum_comparison_boxplot_fixed.png` ✅ **USE THIS**

### Main Analysis Visualizations
- `risk_stratification_fixed.png` ✅ **USE THIS**
- `gplvm_trajectory_analysis_fixed.png` ✅ **USE THIS**
- `complete_results_summary_fixed.png` ✅ **USE THIS**

**Recommendation:** Always use the `*_fixed.png` versions - they have no blank spaces and better layouts!

---

## Manual Option: Run Scripts Individually

If you prefer to run scripts one at a time:

```bash
# Step 1: Forum visualizations
python3 run_visualizations_from_csv.py

# Step 2: Improved forum visualizations
python3 review_and_fix_visualizations.py

# Step 3: Fix existing visualizations
python3 fix_existing_visualizations.py
```

---

## Requirements

- Python 3 (any version works, even 3.14)
- Required packages (usually already installed):
  - pandas
  - numpy
  - matplotlib
  - scipy

All scripts work locally - no Jupyter or h5ad files needed! ✅

---

## Quick Reference

**To regenerate everything:**
```bash
bash run_all_visualizations.sh
```

**To run individual scripts:**
```bash
python3 run_visualizations_from_csv.py          # Forum visuals
python3 review_and_fix_visualizations.py        # Fixed forum visuals
python3 fix_existing_visualizations.py          # Fix all existing visuals
```

That's it! 🎉

