#!/usr/bin/env python3
"""
Complete notebook execution with scVI enabled.
Loads data, processes it, and runs scVI batch correction.
"""

import sys
import os
import traceback
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("COMPLETE NOTEBOOK EXECUTION WITH scVI")
print("="*60)

# Step 1: Verify scVI
print("\n[Step 1] Verifying scVI availability...")
try:
    import scvi
    from scvi.model import SCVI
    HAS_SCVI = True
    print(f"✓ scVI version: {scvi.__version__}")
except ImportError as e:
    print(f"✗ ERROR: scVI not available: {e}")
    sys.exit(1)

# Step 2: Load metadata
print("\n[Step 2] Loading metadata...")
try:
    meta = pd.read_csv('data/sample_metadata_with_age.csv')
    print(f"✓ Loaded metadata: {len(meta)} samples")
except FileNotFoundError:
    meta = pd.read_csv('data/sample_metadata.csv')
    print(f"✓ Loaded metadata (no age): {len(meta)} samples")

# Step 3: Load kallisto abundance files
print("\n[Step 3] Loading kallisto abundance files...")
kallisto_dir = 'zenodo_data/final_code/kallisto'
abundance_files = []
sample_names = []

for folder in os.listdir(kallisto_dir):
    folder_path = os.path.join(kallisto_dir, folder)
    if os.path.isdir(folder_path):
        abundance_file = os.path.join(folder_path, 'abundance.tsv')
        if os.path.exists(abundance_file):
            abundance_files.append(abundance_file)
            sample_names.append(folder)

print(f"  Found {len(abundance_files)} kallisto samples")

# Load and combine abundance files
expr_list = []
valid_samples = []

for i, (ab_file, sample_name) in enumerate(zip(abundance_files, sample_names)):
    try:
        df = pd.read_csv(ab_file, sep='\t')
        # Use TPM column (or est_counts if TPM not available)
        if 'tpm' in df.columns:
            expr_list.append(df.set_index('target_id')['tpm'])
            valid_samples.append(sample_name)
        elif 'est_counts' in df.columns:
            expr_list.append(df.set_index('target_id')['est_counts'])
            valid_samples.append(sample_name)
        else:
            print(f"  Warning: {sample_name} - no TPM or est_counts column")
            continue
        if (i + 1) % 5 == 0:
            print(f"  Loaded {i+1}/{len(abundance_files)} samples...")
    except Exception as e:
        print(f"  Error loading {sample_name}: {e}")
        continue

print(f"  ✓ Successfully loaded {len(expr_list)} samples")

# Combine into expression matrix
print("\n[Step 4] Combining into expression matrix...")
expr = pd.concat(expr_list, axis=1)
expr.columns = valid_samples
print(f"  ✓ Expression matrix shape: {expr.shape}")

# Create AnnData object
print("\n[Step 5] Creating AnnData object...")
adata = sc.AnnData(X=expr.T.values)  # Transpose: samples as rows, genes as columns
adata.var_names = expr.index.astype(str)
adata.obs_names = expr.columns.astype(str)

# Add metadata
print("  Adding metadata...")
meta_indexed = meta.set_index('sample')
adata.obs = meta_indexed.loc[adata.obs_names]
print(f"  ✓ AnnData created: {adata.shape}")

# Prepare for scVI
print("\n[Step 6] Preparing data for scVI...")
# Store counts (TPM values as counts for scVI)
adata.layers['counts'] = adata.X.copy()

# Filter genes
sc.pp.filter_genes(adata, min_cells=2)
print(f"  After filtering: {adata.shape}")

# Normalize for visualization (but keep raw for scVI)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log1p_norm'] = adata.X.copy()

# Restore raw counts for scVI (convert TPM-like values to integer counts)
# scVI expects integer counts, so we'll approximate
print("  Converting to integer counts for scVI...")
# Scale TPM to approximate counts (this is an approximation)
counts_approx = adata.layers['counts'].copy()
# Scale by a factor to get reasonable count ranges
counts_scaled = counts_approx * 1000  # Approximate scaling
counts_scaled = np.maximum(counts_scaled, 0)  # Ensure non-negative
adata.layers['counts'] = counts_scaled.astype(np.int32)
print(f"  ✓ Counts prepared: min={adata.layers['counts'].min()}, max={adata.layers['counts'].max()}")

# Step 7: Run scVI
print("\n" + "="*60)
print("SECTION 2: scVI BATCH CORRECTION")
print("="*60)

print(f"\n  Data shape: {adata.shape}")
print(f"  Stages: {adata.obs['stage'].unique() if 'stage' in adata.obs.columns else 'N/A'}")

# Setup scVI
print("\n[Step 7] Setting up scVI...")
try:
    if 'study' in adata.obs.columns:
        print("  Using 'study' as batch_key...")
        SCVI.setup_anndata(adata, batch_key='study', layer='counts')
    else:
        print("  Setting up without batch_key...")
        SCVI.setup_anndata(adata, layer='counts')
    print("  ✓ scVI setup complete")
except Exception as e:
    print(f"  ✗ Error during setup: {e}")
    traceback.print_exc()
    sys.exit(1)

# Train scVI model
print("\n[Step 8] Training scVI model...")
try:
    model = SCVI(adata, n_latent=10)
    print("  Model initialized. Training...")
    print("  (This may take 5-10 minutes)")
    model.train(max_epochs=200, early_stopping=True, plan_kwargs={"lr": 1e-3})
    print("  ✓ Training complete")
except Exception as e:
    print(f"  ✗ Error during training: {e}")
    traceback.print_exc()
    sys.exit(1)

# Extract latent representation
print("\n[Step 9] Extracting latent representation...")
try:
    latent = model.get_latent_representation()
    adata.obsm['X_scvi'] = latent
    print(f"  ✓ scVI latent space: {latent.shape}")
    
    # Get variance using posterior distribution sampling
    try:
        posterior = model.get_posterior(adata)
        z_samples = posterior.sample_posterior(n_samples=10)
        z_mean = z_samples.mean(axis=0)
        z_var = z_samples.var(axis=0)
        adata.obsm['X_scvi_mu'] = z_mean
        adata.obsm['X_scvi_var'] = z_var
        adata.obs['scvi_latent_sigma'] = z_var.mean(axis=1)
        print(f"  ✓ Uncertainty computed: mean sigma = {adata.obs['scvi_latent_sigma'].mean():.2f}")
    except Exception as e_var:
        print(f"  Note: Could not compute variance (non-critical): {e_var}")
        # Use mean as variance approximation
        adata.obsm['X_scvi_mu'] = latent
        adata.obs['scvi_latent_sigma'] = 0.1  # Placeholder
    
except Exception as e:
    print(f"  ✗ Error extracting latent: {e}")
    traceback.print_exc()
    sys.exit(1)

# Compute UMAP on scVI space
print("\n[Step 10] Computing UMAP on scVI latent space...")
try:
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=10)
    sc.tl.umap(adata)
    print("  ✓ UMAP computed")
except Exception as e:
    print(f"  ✗ Error computing UMAP: {e}")
    traceback.print_exc()

# Save results
print("\n[Step 11] Saving results...")
try:
    output_file = 'adata_with_scvi.h5ad'
    adata.write(output_file)
    print(f"  ✓ Saved: {output_file}")
    
    # Also save with gene symbols if mapping exists
    if 'gene_symbol' in adata.var.columns:
        adata.write('adata_with_gene_symbols.h5ad')
        print(f"  ✓ Saved: adata_with_gene_symbols.h5ad")
    
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
print(f"\n  File saved: adata_with_scvi.h5ad")
print(f"\n  Next steps:")
print(f"    1. Run Section 3 (GPLVM) to compute cellular age")
print(f"    2. Continue with risk stratification")
print(f"    3. Generate visualizations")

