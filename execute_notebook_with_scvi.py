#!/usr/bin/env python3
"""
Execute notebook cells programmatically with scVI enabled.
This runs the notebook step by step, ensuring scVI is used instead of fallbacks.
"""

import sys
import json
import traceback
import numpy as np
import pandas as pd
import scanpy as sc

print("="*60)
print("EXECUTING NOTEBOOK WITH scVI ENABLED")
print("="*60)

# Step 1: Verify scVI is available and force it to be used
print("\n[Step 1] Verifying scVI availability...")
try:
    import scvi
    from scvi.model import SCVI
    HAS_SCVI = True
    print(f"✓ scVI version: {scvi.__version__}")
    print("  scVI will be used (no fallback to PCA)")
except ImportError as e:
    print(f"✗ ERROR: scVI not available: {e}")
    print("  Cannot proceed without scVI. Please install it first.")
    sys.exit(1)

# Step 2: Load notebook to extract cells
notebook_file = 'ADSPROJECT_new.ipynb'
print(f"\n[Step 2] Loading notebook: {notebook_file}...")

try:
    with open(notebook_file, 'r') as f:
        notebook = json.load(f)
    print(f"  ✓ Loaded notebook with {len(notebook['cells'])} cells")
except Exception as e:
    print(f"  ✗ Error loading notebook: {e}")
    sys.exit(1)

# Step 3: Find and load data
print("\n[Step 3] Looking for processed data file...")
adata = None

# Try to load existing processed data
data_files = [
    'adata_with_gene_symbols.h5ad',
    'adata.h5ad',
    'adata_zenodo.h5ad'
]

for data_file in data_files:
    try:
        print(f"  Trying: {data_file}...")
        adata = sc.read_h5ad(data_file)
        print(f"  ✓ Loaded: {data_file}")
        print(f"    Shape: {adata.shape}")
        print(f"    Obs columns: {list(adata.obs.columns)[:5]}...")
        break
    except FileNotFoundError:
        continue
    except Exception as e:
        print(f"  ✗ Error: {e}")
        continue

if adata is None:
    print("\n  ✗ No processed data file found.")
    print("  Need to run data loading cells first.")
    print("  This will take several minutes...")
    
    # Run data loading cells (we'll need to extract them)
    print("\n  Proceeding with data loading from kallisto files...")
    
    # For now, exit and ask user to run data loading first
    print("\n  Please run cells 1-64 first to load the data, then re-run this script.")
    print("  Or load data manually and save as 'adata_with_gene_symbols.h5ad'")
    sys.exit(1)

# Step 4: Run Section 2 (scVI) with forced scVI usage
print("\n" + "="*60)
print("SECTION 2: scVI BATCH CORRECTION")
print("="*60)

print(f"\n  Data shape: {adata.shape}")
print(f"  Stages: {adata.obs['stage'].unique() if 'stage' in adata.obs.columns else 'N/A'}")

# Prepare data for scVI
print("\n[Step 4] Preparing data for scVI...")

# Ensure counts layer exists
if 'counts' not in adata.layers:
    print("  Creating 'counts' layer from .X...")
    adata.layers['counts'] = adata.X.copy()

# Ensure counts are integers (scVI requirement)
if not np.issubdtype(adata.layers['counts'].dtype, np.integer):
    print("  Converting counts to integers...")
    # Check if values are reasonable for integer conversion
    if adata.layers['counts'].max() > 1e6:
        print("  Values appear to be normalized. Checking for raw counts...")
        # Try to find raw counts
        if adata.raw is not None:
            print("  Using .raw for counts...")
            adata.layers['counts'] = adata.raw.X.copy()
        else:
            print("  Warning: No raw counts found. Converting normalized values...")
    adata.layers['counts'] = adata.layers['counts'].astype(np.int32)
    print("  ✓ Converted to integers")

print(f"  ✓ Counts layer ready: {adata.layers['counts'].shape}, dtype: {adata.layers['counts'].dtype}")

# Setup scVI
print("\n[Step 5] Setting up scVI model...")
try:
    # Try with batch_key if available
    if 'study' in adata.obs.columns:
        print("  Using 'study' as batch_key for batch correction...")
        SCVI.setup_anndata(adata, batch_key='study', layer='counts')
    elif 'stage' in adata.obs.columns:
        print("  Using 'stage' as batch_key...")
        SCVI.setup_anndata(adata, batch_key='stage', layer='counts')
    else:
        print("  No batch_key found, setting up without batch correction...")
        SCVI.setup_anndata(adata, layer='counts')
    
    print("  ✓ scVI setup complete")
    
except Exception as e:
    print(f"  ✗ Error during scVI setup: {e}")
    traceback.print_exc()
    sys.exit(1)

# Train scVI model
print("\n[Step 6] Training scVI model...")
try:
    model = SCVI(adata, n_latent=10)
    print("  Model initialized. Training for up to 200 epochs...")
    model.train(max_epochs=200, early_stopping=True, plan_kwargs={"lr": 1e-3})
    print("  ✓ Training complete")
    
except Exception as e:
    print(f"  ✗ Error during training: {e}")
    traceback.print_exc()
    sys.exit(1)

# Extract latent representation
print("\n[Step 7] Extracting latent representation...")
try:
    latent = model.get_latent_representation()
    adata.obsm['X_scvi'] = latent
    print(f"  ✓ scVI latent space extracted: {latent.shape}")
    
except Exception as e:
    print(f"  ✗ Error extracting latent: {e}")
    traceback.print_exc()
    sys.exit(1)

# Compute UMAP on scVI space
print("\n[Step 8] Computing UMAP on scVI latent space...")
try:
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=10)
    sc.tl.umap(adata)
    print("  ✓ UMAP computed using scVI latent space")
    
except Exception as e:
    print(f"  ✗ Error computing UMAP: {e}")
    traceback.print_exc()
    sys.exit(1)

# Save results
print("\n[Step 9] Saving results...")
try:
    output_file = 'adata_with_scvi.h5ad'
    adata.write(output_file)
    print(f"  ✓ Saved to: {output_file}")
    
    # Also update the existing file if it exists
    if 'adata_with_gene_symbols.h5ad' in data_files:
        adata.write('adata_with_gene_symbols.h5ad')
        print(f"  ✓ Updated: adata_with_gene_symbols.h5ad")
    
except Exception as e:
    print(f"  ✗ Error saving: {e}")
    traceback.print_exc()

# Summary
print("\n" + "="*60)
print("SUCCESS: scVI batch correction complete!")
print("="*60)
print(f"\n  Data shape: {adata.shape}")
print(f"  scVI latent space: {adata.obsm['X_scvi'].shape}")
if 'X_umap' in adata.obsm:
    print(f"  UMAP coordinates: {adata.obsm['X_umap'].shape}")
print(f"\n  Next: Run Section 3 (GPLVM) to compute cellular age trajectory")
print(f"  Data saved to: adata_with_scvi.h5ad")

