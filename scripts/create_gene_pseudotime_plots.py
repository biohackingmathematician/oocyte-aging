#!/usr/bin/env python3
"""
Create gene-pseudotime correlation plots for post-hoc analysis.
Plots expression vs cellular age z for top trajectory-associated genes.
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("="*70)
print("GENE-PSEUDOTIME CORRELATION ANALYSIS")
print("="*70)

# Create output directory
os.makedirs('visualizations', exist_ok=True)

# ============================================================================
# Load data
# ============================================================================
print("\nLoading data...")

try:
    import scanpy as sc
    import anndata as ad
    
    # Try to load the complete annotated data
    if os.path.exists('adata_complete_scvi.h5ad'):
        adata = sc.read_h5ad('adata_complete_scvi.h5ad')
        print(f"✓ Loaded adata_complete_scvi.h5ad: {adata.shape}")
    elif os.path.exists('adata_with_scvi.h5ad'):
        adata = sc.read_h5ad('adata_with_scvi.h5ad')
        print(f"✓ Loaded adata_with_scvi.h5ad: {adata.shape}")
    else:
        print("✗ ERROR: No AnnData file found. Please run pipeline first.")
        sys.exit(1)
    
    # Check required fields
    if 'cellular_age_z' not in adata.obs.columns:
        print("✗ ERROR: 'cellular_age_z' not found in adata.obs")
        sys.exit(1)
    
    if 'log1p_norm' not in adata.layers:
        print("  Creating log1p_norm layer...")
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        adata.layers['log1p_norm'] = adata_norm.X.copy()
        print("  ✓ Created log1p_norm layer")
    
    print(f"  Observations: {adata.n_obs}")
    print(f"  Variables: {adata.n_vars}")
    
except Exception as e:
    print(f"✗ Error loading data: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# Compute correlations for all genes
# ============================================================================
print("\n" + "="*70)
print("Computing gene-trajectory correlations")
print("="*70)

# Get expression matrix (log-normalized)
if 'log1p_norm' in adata.layers:
    expr_matrix = adata.layers['log1p_norm']
else:
    expr_matrix = adata.X

# Get cellular age z
z = adata.obs['cellular_age_z'].values

# Compute Spearman correlations for all genes
print(f"  Computing Spearman correlations for {adata.n_vars} genes...")
correlations = []
p_values = []
gene_names = []

for i in range(adata.n_vars):
    gene_expr = expr_matrix[:, i]
    # Remove any NaN or infinite values
    valid_mask = np.isfinite(gene_expr) & np.isfinite(z)
    if valid_mask.sum() < 3:  # Need at least 3 valid points
        correlations.append(np.nan)
        p_values.append(1.0)
    else:
        rho, pval = spearmanr(gene_expr[valid_mask], z[valid_mask])
        correlations.append(rho)
        p_values.append(pval)
    
    gene_names.append(adata.var_names[i])

# Create DataFrame
gene_df = pd.DataFrame({
    'gene': gene_names,
    'correlation': correlations,
    'p_value': p_values
})

# Remove NaN correlations
gene_df = gene_df.dropna(subset=['correlation'])

# Apply FDR correction
gene_df['fdr'] = multipletests(gene_df['p_value'].values, method='fdr_bh')[1]

# Sort by absolute correlation
gene_df['abs_correlation'] = gene_df['correlation'].abs()
gene_df = gene_df.sort_values('abs_correlation', ascending=False)

print(f"  ✓ Computed correlations for {len(gene_df)} genes")
print(f"  Top correlation: {gene_df.iloc[0]['correlation']:.3f}")
print(f"  Genes with |ρ| > 0.7: {(gene_df['abs_correlation'] > 0.7).sum()}")
print(f"  Genes with |ρ| > 0.7 and FDR < 0.1: {((gene_df['abs_correlation'] > 0.7) & (gene_df['fdr'] < 0.1)).sum()}")

# Save full results
gene_df.to_csv('pipeline_results_scvi/tables/gene_trajectory_correlations.csv', index=False)
print("  ✓ Saved full correlation results to pipeline_results_scvi/tables/gene_trajectory_correlations.csv")

# ============================================================================
# Select top genes for plotting
# ============================================================================
print("\n" + "="*70)
print("Selecting top genes for visualization")
print("="*70)

# Select top decreasing and increasing genes
top_decreasing = gene_df[gene_df['correlation'] < 0].head(3)
top_increasing = gene_df[gene_df['correlation'] > 0].head(3)

print("\nTop Decreasing Genes (GV→MI):")
for idx, row in top_decreasing.iterrows():
    print(f"  {row['gene']}: ρ = {row['correlation']:.3f}, p = {row['p_value']:.4f}, FDR = {row['fdr']:.4f}")

print("\nTop Increasing Genes (GV→MI):")
for idx, row in top_increasing.iterrows():
    print(f"  {row['gene']}: ρ = {row['correlation']:.3f}, p = {row['p_value']:.4f}, FDR = {row['fdr']:.4f}")

# Combine for plotting
top_genes = pd.concat([top_decreasing, top_increasing])

# ============================================================================
# Create visualization
# ============================================================================
print("\n" + "="*70)
print("Creating gene-pseudotime plots")
print("="*70)

# Create figure with subplots
n_genes = len(top_genes)
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i, (idx, row) in enumerate(top_genes.iterrows()):
    gene = row['gene']
    rho = row['correlation']
    pval = row['p_value']
    fdr = row['fdr']
    
    # Get gene expression
    gene_idx = adata.var_names.get_loc(gene)
    gene_expr = expr_matrix[:, gene_idx]
    
    # Create scatter plot with smooth curve
    ax = axes[i]
    
    # Scatter plot
    scatter = ax.scatter(z, gene_expr, c=adata.obs['cellular_age_uncertainty'] if 'cellular_age_uncertainty' in adata.obs.columns else 'blue',
                        alpha=0.7, s=100, cmap='viridis' if 'cellular_age_uncertainty' in adata.obs.columns else None)
    
    # Add smooth curve (lowess/loess)
    try:
        from scipy.interpolate import interp1d
        # Sort by z for smooth curve
        sort_idx = np.argsort(z)
        z_sorted = z[sort_idx]
        expr_sorted = gene_expr[sort_idx]
        
        # Use moving average for smoothness
        window_size = max(3, len(z) // 5)
        if window_size % 2 == 0:
            window_size += 1
        
        from scipy.ndimage import uniform_filter1d
        expr_smooth = uniform_filter1d(expr_sorted, size=window_size, mode='nearest')
        
        ax.plot(z_sorted, expr_smooth, 'r-', linewidth=2, alpha=0.8, label='Trend')
        ax.legend(fontsize=8)
    except:
        # Fallback: simple linear fit
        z_valid = z[np.isfinite(gene_expr)]
        expr_valid = gene_expr[np.isfinite(gene_expr)]
        if len(z_valid) > 1:
            z_sorted_idx = np.argsort(z_valid)
            z_sorted = z_valid[z_sorted_idx]
            expr_sorted = expr_valid[z_sorted_idx]
            z_line = np.linspace(z_sorted.min(), z_sorted.max(), 100)
            expr_line = np.interp(z_line, z_sorted, expr_sorted)
            ax.plot(z_line, expr_line, 'r-', linewidth=2, alpha=0.8, label='Interpolated')
            ax.legend(fontsize=8)
    
    # Labels and title
    direction = "Decreasing" if rho < 0 else "Increasing"
    ax.set_xlabel('Cellular Age (z)', fontsize=10)
    ax.set_ylabel('Log-normalized Expression', fontsize=10)
    ax.set_title(f'{gene}\n{direction} (ρ = {rho:.3f}, FDR = {fdr:.3f})', fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Add colorbar for uncertainty if available
    if 'cellular_age_uncertainty' in adata.obs.columns and i == 0:
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Uncertainty (σ)', fontsize=8)

plt.tight_layout()
plt.savefig('visualizations/gene_pseudotime_correlations.png', dpi=300, bbox_inches='tight')
print("  ✓ Saved visualization to visualizations/gene_pseudotime_correlations.png")

# ============================================================================
# Summary statistics
# ============================================================================
print("\n" + "="*70)
print("Summary Statistics")
print("="*70)

n_significant = ((gene_df['abs_correlation'] > 0.7) & (gene_df['fdr'] < 0.1)).sum()
print(f"\nGenes with |ρ| > 0.7 and FDR < 0.1: {n_significant}")

print(f"\nTop 10 genes by absolute correlation:")
print(gene_df.head(10)[['gene', 'correlation', 'p_value', 'fdr']].to_string(index=False))

print("\n" + "="*70)
print("COMPLETE")
print("="*70)

