#!/usr/bin/env python3
"""
Run notebook cells programmatically with scVI enabled.
This script executes the key notebook cells and ensures scVI is used.
"""

import sys
import traceback
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scvi.model import SCVI

print("="*60)
print("RUNNING NOTEBOOK WITH scVI")
print("="*60)

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Step 1: Verify scVI is available
print("\n[Step 1] Verifying scVI availability...")
try:
    import scvi
    print(f"✓ scVI version: {scvi.__version__}")
    HAS_SCVI = True
except ImportError as e:
    print(f"✗ scVI not available: {e}")
    HAS_SCVI = False
    sys.exit(1)

# Step 2: Check for existing processed data
print("\n[Step 2] Checking for existing data files...")
data_files = [
    'adata_with_gene_symbols.h5ad',
    'adata.h5ad',
    'adata_zenodo.h5ad'
]

adata = None
for data_file in data_files:
    try:
        print(f"  Trying to load: {data_file}...")
        adata = sc.read_h5ad(data_file)
        print(f"  ✓ Loaded: {data_file}")
        print(f"    Shape: {adata.shape}")
        break
    except FileNotFoundError:
        continue
    except Exception as e:
        print(f"  ✗ Error loading {data_file}: {e}")
        continue

if adata is None:
    print("\n[Step 3] Data file not found. Need to load raw data first.")
    print("  Please run the data loading cells from the notebook first.")
    print("  Looking for kallisto files...")
    
    import os
    kallisto_dir = 'zenodo_data/final_code/kallisto'
    if os.path.exists(kallisto_dir):
        print(f"  Found kallisto directory: {kallisto_dir}")
        print("  Data loading will be needed from notebook cells 1-64")
    else:
        print("  Kallisto directory not found.")
    sys.exit(1)

# Step 3: Ensure counts layer exists
print("\n[Step 3] Checking data preparation...")
if 'counts' not in adata.layers:
    print("  Creating 'counts' layer from .X...")
    adata.layers['counts'] = adata.X.copy()
    print("  ✓ Created counts layer")

# Ensure counts are integers (scVI requirement)
if not np.issubdtype(adata.layers['counts'].dtype, np.integer):
    print("  Converting counts to integers...")
    adata.layers['counts'] = adata.layers['counts'].astype(np.int32)
    print("  ✓ Converted to integers")

# Step 4: Setup scVI
print("\n[Step 4] Setting up scVI...")
try:
    # Check if batch_key exists
    if 'study' in adata.obs.columns:
        print("  Using 'study' as batch_key...")
        SCVI.setup_anndata(adata, batch_key='study', layer='counts')
    else:
        print("  No batch_key found, setting up without batch correction...")
        SCVI.setup_anndata(adata, layer='counts')
    
    print("  ✓ scVI setup complete")
    
    # Step 5: Train scVI model
    print("\n[Step 5] Training scVI model...")
    model = SCVI(adata, n_latent=10)
    print("  Model initialized. Training...")
    model.train(max_epochs=200, early_stopping=True, plan_kwargs={"lr": 1e-3})
    print("  ✓ Training complete")
    
    # Step 6: Extract latent representation
    print("\n[Step 6] Extracting latent representation...")
    latent = model.get_latent_representation()
    adata.obsm['X_scvi'] = latent
    print(f"  ✓ scVI latent space shape: {latent.shape}")
    
    # Step 7: Compute UMAP on scVI space
    print("\n[Step 7] Computing UMAP on scVI latent space...")
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=10)
    sc.tl.umap(adata)
    print("  ✓ UMAP computed")
    
    # Step 8: Save results
    print("\n[Step 8] Saving results...")
    output_file = 'adata_with_scvi.h5ad'
    adata.write(output_file)
    print(f"  ✓ Saved to: {output_file}")
    
    print("\n" + "="*60)
    print("SUCCESS: scVI batch correction complete!")
    print("="*60)
    print(f"\nscVI latent space: {adata.obsm['X_scvi'].shape}")
    print(f"UMAP coordinates: {adata.obsm['X_umap'].shape}")
    print(f"\nData saved to: {output_file}")
    
except Exception as e:
    print(f"\n✗ Error during scVI processing:")
    traceback.print_exc()
    sys.exit(1)

