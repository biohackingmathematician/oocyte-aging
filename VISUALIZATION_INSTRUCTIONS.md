# Project Forum Visualization Instructions

## Issue Encountered

The visualization codes have been created in the notebook (Cells 88-89), but running them from the command line encountered Python 3.14 compatibility issues with `numba` and other dependencies required by `scanpy`.

## Solution: Run in Jupyter Notebook

The visualization codes are ready to run in your Jupyter notebook environment where `scanpy` and other dependencies are properly configured.

### Steps:

1. **Open the notebook** `ADSPROJECT_new.ipynb` in Jupyter
2. **Navigate to the "PROJECT FORUM: DATA EDA VISUALIZATIONS" section** (Cells 88-89)
3. **Ensure your data file exists** (one of):
   - `adata_trajectory_complete.h5ad`
   - `adata_with_pathway_scores.h5ad`
   - `adata_final_with_intervention.h5ad`
4. **Run Cell 88** - This generates `forum_umap_by_stage.png`
5. **Run Cell 89** - This generates `forum_health_score_boxplot.png`

### What the Visualizations Show:

1. **UMAP Plot (Cell 88)**:
   - 2D UMAP embedding of oocytes
   - Colored by developmental stage (GV = green, MI = red)
   - Shows separation between stages

2. **Health Score Boxplot (Cell 89)**:
   - Side-by-side comparison of health scores
   - GV vs MI stages
   - Includes statistical test (Mann-Whitney U)
   - Shows significance bracket if p < 0.05

### If Data Files Don't Exist:

The visualization cells will automatically:
- Try multiple file paths
- Compute UMAP if missing (from PCA)
- Compute health scores if missing (from pathway scores)

### Alternative: Use Existing Images

If you need to proceed without running the code, check these existing files:
- `complete_results_summary.png` - May contain UMAP visualizations
- `risk_stratification.png` - Contains risk group visualizations

### Code Location:

The visualization code is in:
- **Cell 88**: UMAP plot colored by stage
- **Cell 89**: Health score boxplot GV vs MI

Both cells are fully commented and ready to execute.

