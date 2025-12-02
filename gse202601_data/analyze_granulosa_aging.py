#!/usr/bin/env python3
"""
GSE202601 Granulosa Cell Aging Analysis

This script analyzes granulosa cell transcriptomes from young vs aged ovaries
to identify aging signatures that could serve as clinically accessible biomarkers.

Clinical relevance: Granulosa cells are routinely collected and discarded during
IVF procedures, making them an accessible source for non-invasive aging assessment.
"""

import sys
import os
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Set base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GSE_DIR = Path(BASE_DIR) / "gse202601_data"
RESULTS_DIR = GSE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

def load_granulosa_data(h5ad_path):
    """Load and filter for granulosa cells only."""
    adata = sc.read_h5ad(h5ad_path)
    
    # Filter for granulosa cells (adjust based on actual cell type annotations)
    granulosa_markers = ['FOXL2', 'AMH', 'FSHR', 'CYP19A1', 'INHA']
    
    # Subset to granulosa cells if cell type annotation exists
    if 'cell_type' in adata.obs.columns:
        adata_gc = adata[adata.obs['cell_type'].str.contains('granulosa', case=False, na=False)]
    else:
        print("Warning: No cell_type annotation. Using all cells.")
        adata_gc = adata
    
    return adata_gc

def compute_aging_trajectory(adata, age_column='age_group'):
    """
    Compute aging trajectory for granulosa cells.
    Mirrors the approach used for oocytes in our main analysis.
    """
    # Preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    # Dimensionality reduction
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)
    
    # Compute pseudotime/trajectory using diffusion pseudotime
    # Anchor on young samples
    sc.tl.diffmap(adata)
    
    # If age groups are available, use young cells as root
    if age_column in adata.obs.columns:
        young_mask = adata.obs[age_column].str.contains('young', case=False, na=False)
        if young_mask.sum() > 0:
            # Use first young cell as root (or centroid)
            adata.uns['iroot'] = np.where(young_mask)[0][0]
            sc.tl.dpt(adata)
    
    return adata

def correlate_genes_with_trajectory(adata, trajectory_col='dpt_pseudotime'):
    """
    Find genes correlated with aging trajectory.
    Returns DataFrame with correlation statistics.
    """
    if trajectory_col not in adata.obs.columns:
        print(f"Warning: {trajectory_col} not found. Skipping correlation.")
        return None
    
    trajectory = adata.obs[trajectory_col].values
    
    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        expr_matrix = adata.X.toarray()
    else:
        expr_matrix = adata.X
    
    correlations = []
    for i, gene in enumerate(adata.var_names):
        expr = expr_matrix[:, i]
        
        # Skip low-expression genes
        if np.mean(expr > 0) < 0.1:
            continue
            
        rho, pval = stats.spearmanr(expr, trajectory, nan_policy='omit')
        
        correlations.append({
            'gene': gene,
            'spearman_rho': rho,
            'pvalue': pval,
            'mean_expr': np.mean(expr),
            'pct_expressed': np.mean(expr > 0) * 100
        })
    
    corr_df = pd.DataFrame(correlations)
    
    # FDR correction
    from statsmodels.stats.multitest import multipletests
    _, corr_df['fdr'], _, _ = multipletests(corr_df['pvalue'], method='fdr_bh')
    
    return corr_df.sort_values('pvalue')

def compare_with_oocyte_signatures(gc_correlations, oocyte_correlations_path):
    """
    Compare granulosa cell aging genes with oocyte aging genes.
    Identifies shared vs cell-type-specific aging signatures.
    """
    # Load oocyte correlations
    oocyte_corr = pd.read_csv(oocyte_correlations_path)
    
    # Define significant genes
    gc_sig = set(gc_correlations[
        (gc_correlations['fdr'] < 0.1) & 
        (abs(gc_correlations['spearman_rho']) > 0.3)
    ]['gene'])
    
    oocyte_sig = set(oocyte_corr[
        (oocyte_corr['fdr'] < 0.1) & 
        (abs(oocyte_corr['spearman_rho']) > 0.3)
    ]['gene'])
    
    # Compute overlaps
    shared = gc_sig & oocyte_sig
    gc_specific = gc_sig - oocyte_sig
    oocyte_specific = oocyte_sig - gc_sig
    
    comparison = {
        'shared_aging_genes': len(shared),
        'granulosa_specific': len(gc_specific),
        'oocyte_specific': len(oocyte_specific),
        'jaccard_index': len(shared) / len(gc_sig | oocyte_sig) if (gc_sig | oocyte_sig) else 0,
        'shared_genes': list(shared)[:50],  # Top 50 for display
        'gc_specific_genes': list(gc_specific)[:50],
        'oocyte_specific_genes': list(oocyte_specific)[:50]
    }
    
    return comparison

def generate_clinical_biomarker_candidates(gc_correlations, min_rho=0.5, max_fdr=0.05):
    """
    Identify granulosa cell genes suitable as clinical biomarkers.
    
    Criteria:
    - Strong correlation with aging trajectory
    - Highly expressed (detectable in clinical samples)
    - Known biological relevance to ovarian function
    """
    
    # Known ovarian/reproductive genes for prioritization
    reproductive_genes = {
        'AMH', 'FSHR', 'LHCGR', 'CYP19A1', 'CYP11A1', 'STAR', 'INHA', 'INHBA', 
        'FOXL2', 'BMP15', 'GDF9', 'KIT', 'KITLG', 'WNT4', 'RSPO1', 'CTNNB1',
        'TP53', 'CDKN1A', 'CDKN2A', 'SIRT1', 'SIRT3', 'SOD1', 'SOD2', 'CAT',
        'GPX1', 'PRDX1', 'NRF2', 'KEAP1', 'TFAM', 'PGC1A', 'NRF1'
    }
    
    candidates = gc_correlations[
        (abs(gc_correlations['spearman_rho']) >= min_rho) &
        (gc_correlations['fdr'] <= max_fdr) &
        (gc_correlations['pct_expressed'] >= 20)  # Detectable in most cells
    ].copy()
    
    # Flag known reproductive genes
    candidates['known_reproductive'] = candidates['gene'].isin(reproductive_genes)
    
    # Sort by absolute correlation
    candidates['abs_rho'] = abs(candidates['spearman_rho'])
    candidates = candidates.sort_values(['known_reproductive', 'abs_rho'], ascending=[False, False])
    
    return candidates

if __name__ == "__main__":
    print("GSE202601 Granulosa Cell Analysis")
    print("=" * 50)
    print("\nThis script requires the processed h5ad file from GSE202601.")
    print("Please download from GEO and place in gse202601_data/ directory.")
    print("\nExpected file: gse202601_data/GSE202601_snRNA_processed.h5ad")
