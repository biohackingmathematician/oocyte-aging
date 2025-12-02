# scVI Setup Complete 

## Summary

The notebook environment has been successfully configured and tested with scVI. The data has been processed and scVI batch correction has been completed.

## What Was Done

1. **Environment Deployed**: Python 3.11.9 virtual environment with all required packages
   - scVI-tools: 1.4.0.post1
   - scanpy: 1.11.5
   - TensorFlow: 2.15.1
   - GPflow: 2.10.0

2. **Data Processed**: 
   - Loaded 20 kallisto samples from `zenodo_data/final_code/kallisto/`
   - Created AnnData object with 20 samples × 126,966 genes
   - Applied scVI batch correction
   - Extracted 10-dimensional latent representation
   - Computed UMAP coordinates

3. **Output Files Created**:
   - `adata_with_scvi.h5ad` - Complete AnnData object with scVI latent space

## Using the Notebook

### Option 1: Load Pre-processed Data (Recommended)

If you just want to use the scVI-processed data, you can load it directly:

```python
import scanpy as sc
adata = sc.read_h5ad('adata_with_scvi.h5ad')

# scVI latent space is available at:
print(adata.obsm['X_scvi'].shape)  # (20, 10)

# UMAP coordinates:
print(adata.obsm['X_umap'].shape)  # (20, 2)
```

### Option 2: Run Notebook Cells

The notebook will automatically detect scVI and use it instead of PCA fallback. The relevant cells are:

- **Section 2: scVI Batch Correction** (around cell 60-65)
  - Automatically detects scVI availability
  - Uses scVI for batch correction if available
  - Falls back to PCA only if scVI is not installed

### Option 3: Re-run Complete Pipeline

To reprocess the data from scratch:

```bash
source venv/bin/activate
python run_notebook_complete.py
```

This will:
1. Load kallisto abundance files
2. Create AnnData object
3. Train scVI model (takes ~2-5 minutes)
4. Extract latent representation
5. Compute UMAP
6. Save to `adata_with_scvi.h5ad`

## Verification

To verify everything is working:

```bash
source venv/bin/activate
python test_scvi_in_notebook.py
```

This will confirm:
- scVI is installed and importable
- Data file exists and contains scVI latent space
- Environment is ready for notebook execution

## Next Steps

1. **Continue with Section 3 (GPLVM)**: The notebook can now proceed to compute cellular age trajectories using the scVI latent space
2. **Visualize Results**: Use the scVI latent space for UMAP visualization
3. **Risk Stratification**: Continue with pathway scoring and risk group assignment

## Troubleshooting

If the notebook still tries to use PCA fallback:

1. **Check scVI is available in the kernel**:
   ```python
   import scvi
   print(scvi.__version__)
   ```

2. **Ensure you're using the correct kernel**:
   - The notebook should use the `venv` environment
   - In Jupyter: Kernel → Change Kernel → Select Python environment with scVI

3. **Force scVI usage**: The notebook cell checks `HAS_SCVI` flag. If scVI is installed, it will automatically be set to `True`.

## Files Created

- `venv/` - Virtual environment (already in .gitignore)
- `adata_with_scvi.h5ad` - Processed data with scVI latent space
- `run_notebook_complete.py` - Script to run complete pipeline
- `test_scvi_in_notebook.py` - Verification script
- `SCVI_SETUP_COMPLETE.md` - This file

## Notes

- The scVI model was trained with 10 latent dimensions
- Training stopped early due to validation metric, which is normal for small datasets
- The latent representation is saved and ready for downstream analysis
- All original data layers are preserved (counts, log1p_norm, etc.)

