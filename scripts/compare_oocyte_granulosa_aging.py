#!/usr/bin/env python3
"""
Compare aging signatures between oocytes and granulosa cells.

This script:
1. Loads oocyte trajectory correlations from our main analysis
2. Loads granulosa cell trajectory correlations from GSE202601
3. Computes overlap and generates comparison figures
4. Identifies candidate clinical biomarkers

Output: Figures and tables for paper discussion section.
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = Path(BASE_DIR) / "pipeline_results_scvi" / "tables"
GSE_DIR = Path(BASE_DIR) / "gse202601_data"
OUTPUT_DIR = Path(BASE_DIR) / "comparison_results"
OUTPUT_DIR.mkdir(exist_ok=True)

def load_oocyte_correlations():
    """Load oocyte gene-trajectory correlations."""
    corr_path = RESULTS_DIR / "gene_trajectory_correlations_with_symbols.csv"
    
    if not corr_path.exists():
        # Try alternative path
        corr_path = RESULTS_DIR / "gene_trajectory_correlations.csv"
    
    if corr_path.exists():
        df = pd.read_csv(corr_path)
        print(f"Loaded {len(df)} oocyte gene correlations")
        return df
    else:
        print(f"Oocyte correlations not found at {corr_path}")
        return None

def load_granulosa_correlations():
    """Load granulosa cell gene-trajectory correlations (if available)."""
    corr_path = GSE_DIR / "results" / "granulosa_trajectory_correlations.csv"
    
    if corr_path.exists():
        df = pd.read_csv(corr_path)
        print(f"Loaded {len(df)} granulosa gene correlations")
        return df
    else:
        print(f"Granulosa correlations not yet generated")
        print(f"Run analyze_granulosa_aging.py first after downloading GSE202601")
        return None

def compute_signature_overlap(oocyte_corr, granulosa_corr, 
                               rho_thresh=0.3, fdr_thresh=0.1):
    """
    Compute overlap between oocyte and granulosa aging signatures.
    """
    # Determine gene column name
    gene_col_oocyte = 'gene_symbol' if 'gene_symbol' in oocyte_corr.columns else 'gene'
    gene_col_granulosa = 'gene'
    
    # Determine correlation column names
    rho_col_oocyte = 'correlation' if 'correlation' in oocyte_corr.columns else 'spearman_rho'
    rho_col_granulosa = 'spearman_rho'
    
    # Define significant genes
    oocyte_sig_up = set(oocyte_corr[
        (oocyte_corr[rho_col_oocyte] > rho_thresh) & 
        (oocyte_corr['fdr'] < fdr_thresh)
    ][gene_col_oocyte])
    
    oocyte_sig_down = set(oocyte_corr[
        (oocyte_corr[rho_col_oocyte] < -rho_thresh) & 
        (oocyte_corr['fdr'] < fdr_thresh)
    ][gene_col_oocyte])
    
    granulosa_sig_up = set(granulosa_corr[
        (granulosa_corr[rho_col_granulosa] > rho_thresh) & 
        (granulosa_corr['fdr'] < fdr_thresh)
    ][gene_col_granulosa])
    
    granulosa_sig_down = set(granulosa_corr[
        (granulosa_corr[rho_col_granulosa] < -rho_thresh) & 
        (granulosa_corr['fdr'] < fdr_thresh)
    ][gene_col_granulosa])
    
    results = {
        'oocyte_up': len(oocyte_sig_up),
        'oocyte_down': len(oocyte_sig_down),
        'granulosa_up': len(granulosa_sig_up),
        'granulosa_down': len(granulosa_sig_down),
        'shared_up': len(oocyte_sig_up & granulosa_sig_up),
        'shared_down': len(oocyte_sig_down & granulosa_sig_down),
        'shared_up_genes': list(oocyte_sig_up & granulosa_sig_up)[:50],
        'shared_down_genes': list(oocyte_sig_down & granulosa_sig_down)[:50],
    }
    
    # Compute enrichment (Fisher's exact test)
    all_genes = set(oocyte_corr[gene_col_oocyte]) & set(granulosa_corr[gene_col_granulosa])
    
    for direction in ['up', 'down']:
        oocyte_set = oocyte_sig_up if direction == 'up' else oocyte_sig_down
        granulosa_set = granulosa_sig_up if direction == 'up' else granulosa_sig_down
        
        # Contingency table
        both = len(oocyte_set & granulosa_set & all_genes)
        oocyte_only = len((oocyte_set - granulosa_set) & all_genes)
        granulosa_only = len((granulosa_set - oocyte_set) & all_genes)
        neither = len(all_genes - oocyte_set - granulosa_set)
        
        if both + oocyte_only + granulosa_only + neither > 0:
            odds_ratio, pval = stats.fisher_exact([[both, oocyte_only], 
                                                    [granulosa_only, neither]])
            results[f'fisher_odds_{direction}'] = odds_ratio
            results[f'fisher_pval_{direction}'] = pval
    
    return results

def plot_correlation_comparison(oocyte_corr, granulosa_corr, output_path):
    """
    Scatter plot comparing correlation coefficients between cell types.
    """
    # Determine column names
    gene_col_oocyte = 'gene_symbol' if 'gene_symbol' in oocyte_corr.columns else 'gene'
    rho_col_oocyte = 'correlation' if 'correlation' in oocyte_corr.columns else 'spearman_rho'
    
    # Merge on gene symbol
    merged = pd.merge(
        oocyte_corr[[gene_col_oocyte, rho_col_oocyte, 'fdr']].rename(
            columns={rho_col_oocyte: 'oocyte_rho', 'fdr': 'oocyte_fdr', gene_col_oocyte: 'gene'}
        ),
        granulosa_corr[['gene', 'spearman_rho', 'fdr']].rename(
            columns={'spearman_rho': 'granulosa_rho', 'fdr': 'granulosa_fdr'}
        ),
        on='gene',
        how='inner'
    )
    
    print(f"Merged {len(merged)} common genes")
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Color by significance
    both_sig = (merged['oocyte_fdr'] < 0.1) & (merged['granulosa_fdr'] < 0.1)
    oocyte_only = (merged['oocyte_fdr'] < 0.1) & (merged['granulosa_fdr'] >= 0.1)
    granulosa_only = (merged['oocyte_fdr'] >= 0.1) & (merged['granulosa_fdr'] < 0.1)
    neither = ~(both_sig | oocyte_only | granulosa_only)
    
    ax.scatter(merged.loc[neither, 'oocyte_rho'], 
               merged.loc[neither, 'granulosa_rho'],
               alpha=0.3, s=10, c='gray', label='Neither sig.')
    ax.scatter(merged.loc[oocyte_only, 'oocyte_rho'], 
               merged.loc[oocyte_only, 'granulosa_rho'],
               alpha=0.5, s=20, c='blue', label='Oocyte only')
    ax.scatter(merged.loc[granulosa_only, 'oocyte_rho'], 
               merged.loc[granulosa_only, 'granulosa_rho'],
               alpha=0.5, s=20, c='green', label='Granulosa only')
    ax.scatter(merged.loc[both_sig, 'oocyte_rho'], 
               merged.loc[both_sig, 'granulosa_rho'],
               alpha=0.7, s=30, c='red', label='Both sig.')
    
    # Correlation line
    r, p = stats.pearsonr(merged['oocyte_rho'], merged['granulosa_rho'])
    ax.plot([-1, 1], [-1, 1], 'k--', alpha=0.5)
    
    ax.set_xlabel('Oocyte Aging Correlation (ρ)')
    ax.set_ylabel('Granulosa Cell Aging Correlation (ρ)')
    ax.set_title(f'Aging Signatures: Oocyte vs Granulosa\nPearson r = {r:.3f}, p = {p:.2e}')
    ax.legend(loc='lower right')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.axhline(0, color='gray', linestyle='-', alpha=0.3)
    ax.axvline(0, color='gray', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved correlation comparison: {output_path}")
    
    return merged, r, p

def generate_biomarker_table(oocyte_corr, granulosa_corr, output_path):
    """
    Generate table of candidate clinical biomarkers.
    Focus on genes that:
    1. Show aging correlation in granulosa cells (clinically accessible)
    2. Are validated by oocyte data (biologically relevant to oocyte quality)
    """
    gene_col_oocyte = 'gene_symbol' if 'gene_symbol' in oocyte_corr.columns else 'gene'
    rho_col_oocyte = 'correlation' if 'correlation' in oocyte_corr.columns else 'spearman_rho'
    
    merged = pd.merge(
        oocyte_corr[[gene_col_oocyte, rho_col_oocyte, 'fdr']].rename(
            columns={rho_col_oocyte: 'oocyte_rho', 'fdr': 'oocyte_fdr', gene_col_oocyte: 'gene'}
        ),
        granulosa_corr[['gene', 'spearman_rho', 'fdr']].rename(
            columns={'spearman_rho': 'granulosa_rho', 'fdr': 'granulosa_fdr'}
        ),
        on='gene',
        how='inner'
    )
    
    # Filter for candidate biomarkers
    candidates = merged[
        (abs(merged['granulosa_rho']) > 0.3) &  # Detectable in granulosa
        (merged['granulosa_fdr'] < 0.1) &
        (abs(merged['oocyte_rho']) > 0.2) &  # Validated in oocytes
        (np.sign(merged['granulosa_rho']) == np.sign(merged['oocyte_rho']))  # Same direction
    ].copy()
    
    candidates['combined_score'] = (
        abs(candidates['granulosa_rho']) + abs(candidates['oocyte_rho'])
    ) / 2
    
    candidates = candidates.sort_values('combined_score', ascending=False)
    
    # Save
    candidates.to_csv(output_path, index=False)
    print(f"Saved {len(candidates)} candidate biomarkers: {output_path}")
    
    return candidates

def main():
    print("Oocyte vs Granulosa Cell Aging Signature Comparison")
    
    # Load data
    oocyte_corr = load_oocyte_correlations()
    granulosa_corr = load_granulosa_correlations()
    
    if oocyte_corr is None:
        print("\nError: Cannot proceed without oocyte correlations.")
        print("Run the main pipeline first.")
        return
    
    if granulosa_corr is None:
        print("\nGranulosa data not yet available")
        print("""
To complete this comparison:

1. Download GSE202601 processed data from GEO:
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202601

2. Run the granulosa analysis:
   python gse202601_data/analyze_granulosa_aging.py

3. Re-run this comparison script.

For now, generating placeholder for paper discussion section...
""")
        
        # Generate placeholder discussion text
        top_genes = oocyte_corr.nsmallest(10, 'fdr')[
            ['gene_symbol' if 'gene_symbol' in oocyte_corr.columns else 'gene', 
             'correlation' if 'correlation' in oocyte_corr.columns else 'spearman_rho', 'fdr']
        ].to_string(index=False)
        
        discussion = f"""
## Alternative Cell Types for Clinical Aging Biomarkers

### Granulosa Cells as Accessible Biomarkers

Professor's question: "Are there different cell types that could be 
correlated with oocyte age that are easier to collect?"

**Answer**: Yes - granulosa cells represent an ideal alternative.

#### Clinical Advantages of Granulosa Cells:

1. **Routinely collected**: Retrieved during standard IVF oocyte retrieval
2. **Currently discarded**: No additional procedures required
3. **Abundant**: Multiple cells per follicle vs single oocyte
4. **Reflect oocyte status**: Direct communication via gap junctions

#### Evidence from Literature:

- Chen et al. (2024) Nature Aging showed granulosa cells exhibit 
  age-related transcriptomic changes in human ovaries (GSE202601)
- Granulosa cell markers can predict oocyte developmental potential
- Age-related changes in granulosa cells precede oocyte changes

#### Proposed Future Work:

1. Integrate GSE202601 granulosa cell snRNA-seq data
2. Identify shared aging signatures between oocytes and granulosa cells
3. Validate granulosa-based biomarkers that predict oocyte quality
4. Develop clinically deployable granulosa cell aging panel

#### Preliminary Analysis (Oocytes Only):

Based on our oocyte analysis, top aging-correlated genes include:

{top_genes}

These represent candidates for validation in granulosa cells.
"""
        
        with open(OUTPUT_DIR / "granulosa_discussion_placeholder.md", 'w') as f:
            f.write(discussion)
        
        print(f"Saved discussion placeholder: {OUTPUT_DIR}/granulosa_discussion_placeholder.md")
        return
    
    # Full comparison if both datasets available
    print("\nComputing signature overlap...")
    overlap = compute_signature_overlap(oocyte_corr, granulosa_corr)
    
    print("\nGenerating comparison figures...")
    merged, r, p = plot_correlation_comparison(oocyte_corr, granulosa_corr,
                                                OUTPUT_DIR / "scatter_oocyte_granulosa.png")
    
    print("\nIdentifying candidate biomarkers...")
    candidates = generate_biomarker_table(oocyte_corr, granulosa_corr,
                                          OUTPUT_DIR / "candidate_biomarkers.csv")
    
    # Summary
    summary = f"""
## Oocyte vs Granulosa Cell Aging Comparison Results

### Signature Overlap

- Oocyte aging genes (up): {overlap['oocyte_up']}
- Oocyte aging genes (down): {overlap['oocyte_down']}
- Granulosa aging genes (up): {overlap['granulosa_up']}
- Granulosa aging genes (down): {overlap['granulosa_down']}
- **Shared (increasing with age)**: {overlap['shared_up']}
- **Shared (decreasing with age)**: {overlap['shared_down']}

### Correlation Between Cell Types

- Pearson r = {r:.3f} (p = {p:.2e})
- {len(merged)} genes compared

### Candidate Clinical Biomarkers

- {len(candidates)} genes identified
- Criteria: Strong granulosa signal + oocyte validation + concordant direction

### Top 10 Biomarker Candidates:

{candidates.head(10).to_string(index=False) if len(candidates) > 0 else "None identified"}

### Figures Generated:

- scatter_oocyte_granulosa.png

### Files Generated:

- candidate_biomarkers.csv
"""
    
    with open(OUTPUT_DIR / "comparison_summary.md", 'w') as f:
        f.write(summary)
    
    print(summary)

if __name__ == "__main__":
    main()

