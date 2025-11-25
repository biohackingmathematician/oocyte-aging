#!/usr/bin/env python3
"""
Test that scVI is properly configured and can be used in the notebook.
"""

import sys
import json

print("Testing scVI availability for notebook...")

# Test import
try:
    import scvi
    from scvi.model import SCVI
    print(f"✓ scVI version: {scvi.__version__}")
    HAS_SCVI = True
except ImportError as e:
    print(f"✗ scVI not available: {e}")
    HAS_SCVI = False
    sys.exit(1)

# Check if we have the processed data
try:
    import scanpy as sc
    adata = sc.read_h5ad('adata_with_scvi.h5ad')
    print(f"✓ Data file found: {adata.shape}")
    if 'X_scvi' in adata.obsm:
        print(f"✓ scVI latent space present: {adata.obsm['X_scvi'].shape}")
    else:
        print("  Warning: X_scvi not found, but file exists")
except Exception as e:
    print(f"  Note: Could not load data file: {e}")

print("\n" + "="*60)
print("scVI is ready to use in the notebook!")
print("="*60)
print("\nTo ensure the notebook uses scVI:")
print("  1. The notebook cell will detect scVI automatically")
print("  2. It will skip PCA fallback if scVI is available")
print("  3. You can also load 'adata_with_scvi.h5ad' directly")
print("\nIf you want to rerun scVI training:")
print("  python run_notebook_complete.py")

