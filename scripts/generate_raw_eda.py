"""
Generate Raw EDA Visualizations

This script creates comprehensive exploratory data analysis plots from raw,
untransformed expression data including:
- Gene expression heatmaps (before transformation)
- Library size distributions
- Mitochondrial gene percentages
- Top expressed genes by stage
- Raw gene-gene correlation matrices
- Basic boxplots of raw counts

Author: Agna Chan, Aniqa Nayim, Rimjhim Singh
Date: November 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

def load_raw_expression_data():
    """
    Load raw expression data from abundance.tsv files.
    
    Returns
    -------
    expr_df : DataFrame
        Expression matrix (genes x samples)
    metadata : DataFrame
        Sample metadata
    """
    # Try to load from AnnData if it exists
    adata_paths = [
        'adata_zenodo.h5ad',
        '../adata_zenodo.h5ad',
        'data/adata_zenodo.h5ad'
    ]
    
    for path in adata_paths:
        if os.path.exists(path):
            try:
                import scanpy as sc
                adata = sc.read_h5ad(path)
                print(f"Loaded AnnData from {path}")
                
                # Get raw counts if available
                if 'counts' in adata.layers:
                    expr_df = pd.DataFrame(
                        adata.layers['counts'].T,
                        index=adata.var_names,
                        columns=adata.obs_names
                    )
                else:
                    # Use X as raw data
                    expr_df = pd.DataFrame(
                        adata.X.T,
                        index=adata.var_names,
                        columns=adata.obs_names
                    )
                
                metadata = adata.obs.copy()
                return expr_df, metadata
            except Exception as e:
                print(f"Could not load from {path}: {e}")
                break
    
    # Fallback: Load from abundance.tsv files
    print("Loading from abundance.tsv files...")
    kallisto_dir = 'zenodo_data/final_code/kallisto'
    if not os.path.exists(kallisto_dir):
        kallisto_dir = '../zenodo_data/final_code/kallisto'
    
    if not os.path.exists(kallisto_dir):
        raise FileNotFoundError(f"Could not find kallisto directory: {kallisto_dir}")
    
    # Find all abundance.tsv files
    abundance_files = glob.glob(os.path.join(kallisto_dir, '*/abundance.tsv'))
    
    if len(abundance_files) == 0:
        raise FileNotFoundError("No abundance.tsv files found")
    
    print(f"Found {len(abundance_files)} abundance files")
    
    # Load each file and combine
    expr_dict = {}
    for file in abundance_files:
        sample_name = os.path.basename(os.path.dirname(file))
        df = pd.read_csv(file, sep='\t')
        # Use TPM as raw expression (or est_counts if preferred)
        expr_dict[sample_name] = df.set_index('target_id')['tpm'].to_dict()
    
    # Convert to DataFrame
    all_genes = set()
    for expr in expr_dict.values():
        all_genes.update(expr.keys())
    
    expr_df = pd.DataFrame(expr_dict)
    expr_df = expr_df.fillna(0)
    
    # Load metadata
    metadata_paths = [
        'sample_metadata_with_age.csv',
        '../data/sample_metadata_with_age.csv',
        'data/sample_metadata_with_age.csv'
    ]
    
    metadata = None
    for path in metadata_paths:
        if os.path.exists(path):
            metadata = pd.read_csv(path)
            if 'sample' in metadata.columns:
                metadata = metadata.set_index('sample')
            break
    
    if metadata is None:
        # Create minimal metadata from sample names
        metadata = pd.DataFrame(index=expr_df.columns)
        metadata['stage'] = 'Unknown'
        for idx in metadata.index:
            if 'VG' in idx or 'GV' in idx:
                metadata.loc[idx, 'stage'] = 'GV'
            elif 'MI' in idx and 'MII' not in idx:
                metadata.loc[idx, 'stage'] = 'MI'
            elif 'MII' in idx:
                metadata.loc[idx, 'stage'] = 'MII'
    
    return expr_df, metadata

def identify_mitochondrial_genes(gene_names):
    """
    Identify mitochondrial genes.
    
    Parameters
    ----------
    gene_names : array-like
        Gene identifiers (may be transcript IDs or gene symbols)
    
    Returns
    -------
    mt_mask : array
        Boolean mask for mitochondrial genes
    """
    gene_names_str = [str(g).upper() for g in gene_names]
    gene_series = pd.Series(gene_names_str)
    
    mt_mask = np.zeros(len(gene_names), dtype=bool)
    
    # Pattern 1: Direct MT- prefix (human mitochondrial genes)
    mt_mask |= gene_series.str.contains('^MT-', regex=True, na=False).values
    mt_mask |= gene_series.str.contains('^MT_', regex=True, na=False).values
    mt_mask |= gene_series.str.contains('^MT\.', regex=True, na=False).values
    
    # Pattern 2: Contains MT in transcript ID (e.g., ENST...MT...)
    mt_mask |= gene_series.str.contains('MT', regex=False, na=False).values & \
               gene_series.str.contains('ENST', regex=False, na=False).values
    
    # Pattern 3: Common mitochondrial gene symbols/names
    mt_gene_patterns = [
        'COX1', 'COX2', 'COX3', 'COX4', 'COX5', 'COX6', 'COX7', 'COX8',
        'CYTB', 'CYTC',
        'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
        'ATP6', 'ATP8', 'ATP5',
        '12S', '16S', 'RNR1', 'RNR2',
        'TRNA', 'TRN'
    ]
    for pattern in mt_gene_patterns:
        mt_mask |= gene_series.str.contains(pattern, regex=False, na=False).values
    
    # Pattern 4: Chromosome M (mitochondrial chromosome)
    mt_mask |= gene_series.str.contains('CHRM', regex=False, na=False).values
    mt_mask |= gene_series.str.contains('CHRMT', regex=False, na=False).values
    
    return mt_mask

def plot_library_size_distribution(expr_df, metadata, output_dir='visualizations'):
    """
    Plot library size (total counts) distribution per cell.
    
    Parameters
    ----------
    expr_df : DataFrame
        Expression matrix
    metadata : DataFrame
        Sample metadata
    output_dir : str
        Output directory
    """
    # Calculate library size (total expression per sample)
    library_sizes = expr_df.sum(axis=0)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Overall distribution
    ax1 = axes[0]
    ax1.hist(library_sizes.values, bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(library_sizes.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {library_sizes.mean():.0f}')
    ax1.axvline(library_sizes.median(), color='green', linestyle='--', linewidth=2, label=f'Median: {library_sizes.median():.0f}')
    ax1.set_xlabel('Library Size (Total Expression)', fontweight='bold')
    ax1.set_ylabel('Number of Samples', fontweight='bold')
    ax1.set_title('Library Size Distribution (All Samples)', fontweight='bold', fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # By stage
    ax2 = axes[1]
    if 'stage' in metadata.columns:
        stages = metadata['stage'].dropna().unique()
        stage_data = []
        stage_labels = []
        for stage in stages:
            stage_samples = metadata[metadata['stage'] == stage].index
            stage_lib_sizes = library_sizes[stage_samples]
            stage_data.append(stage_lib_sizes.values)
            stage_labels.append(f'{stage}\n(n={len(stage_samples)})')
        
        bp = ax2.boxplot(stage_data, labels=stage_labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
            patch.set_alpha(0.7)
        ax2.set_ylabel('Library Size (Total Expression)', fontweight='bold')
        ax2.set_xlabel('Developmental Stage', fontweight='bold')
        ax2.set_title('Library Size by Stage', fontweight='bold', fontsize=12)
        ax2.grid(True, alpha=0.3, axis='y')
    else:
        ax2.text(0.5, 0.5, 'Stage information\nnot available', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Library Size by Stage', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_eda_library_size.png')
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_file}")

def plot_mitochondrial_percentage(expr_df, metadata, output_dir='visualizations'):
    """
    Plot mitochondrial gene percentage per cell.
    
    Parameters
    ----------
    expr_df : DataFrame
        Expression matrix
    metadata : DataFrame
        Sample metadata
    output_dir : str
        Output directory
    """
    # Identify mitochondrial genes
    mt_mask = identify_mitochondrial_genes(expr_df.index)
    n_mt_genes = mt_mask.sum()
    
    if n_mt_genes == 0:
        print("Warning: No mitochondrial genes identified")
        return
    
    # Calculate mitochondrial percentage per sample
    mt_expr = expr_df.loc[mt_mask].sum(axis=0)
    total_expr = expr_df.sum(axis=0)
    mt_pct = (mt_expr / total_expr * 100).fillna(0)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Overall distribution
    ax1 = axes[0]
    ax1.hist(mt_pct.values, bins=20, edgecolor='black', alpha=0.7, color='coral')
    ax1.axvline(mt_pct.mean(), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mt_pct.mean():.2f}%')
    ax1.axvline(mt_pct.median(), color='green', linestyle='--', linewidth=2, 
                label=f'Median: {mt_pct.median():.2f}%')
    ax1.set_xlabel('Mitochondrial Gene Percentage (%)', fontweight='bold')
    ax1.set_ylabel('Number of Samples', fontweight='bold')
    ax1.set_title(f'Mitochondrial Gene Percentage\n({n_mt_genes} MT genes identified)', 
                  fontweight='bold', fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # By stage
    ax2 = axes[1]
    if 'stage' in metadata.columns:
        stages = metadata['stage'].dropna().unique()
        stage_data = []
        stage_labels = []
        for stage in stages:
            stage_samples = metadata[metadata['stage'] == stage].index
            stage_mt_pct = mt_pct[stage_samples]
            stage_data.append(stage_mt_pct.values)
            stage_labels.append(f'{stage}\n(n={len(stage_samples)})')
        
        bp = ax2.boxplot(stage_data, labels=stage_labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightcoral')
            patch.set_alpha(0.7)
        ax2.set_ylabel('Mitochondrial Gene Percentage (%)', fontweight='bold')
        ax2.set_xlabel('Developmental Stage', fontweight='bold')
        ax2.set_title('Mitochondrial Percentage by Stage', fontweight='bold', fontsize=12)
        ax2.grid(True, alpha=0.3, axis='y')
    else:
        ax2.text(0.5, 0.5, 'Stage information\nnot available', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Mitochondrial Percentage by Stage', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_eda_mitochondrial.png')
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_file}")

def plot_top_genes_by_stage(expr_df, metadata, n_top=20, output_dir='visualizations'):
    """
    Plot top expressed genes by stage.
    
    Parameters
    ----------
    expr_df : DataFrame
        Expression matrix
    metadata : DataFrame
        Sample metadata
    n_top : int
        Number of top genes to show
    output_dir : str
        Output directory
    """
    if 'stage' not in metadata.columns:
        print("Warning: Stage information not available")
        return
    
    stages = metadata['stage'].dropna().unique()
    n_stages = len(stages)
    
    fig, axes = plt.subplots(1, n_stages, figsize=(6*n_stages, 8))
    if n_stages == 1:
        axes = [axes]
    
    for idx, stage in enumerate(stages):
        stage_samples = metadata[metadata['stage'] == stage].index
        stage_expr = expr_df[stage_samples]
        
        # Calculate mean expression per gene
        mean_expr = stage_expr.mean(axis=1).sort_values(ascending=False)
        top_genes = mean_expr.head(n_top)
        
        ax = axes[idx]
        y_pos = np.arange(len(top_genes))
        bars = ax.barh(y_pos, top_genes.values, color='steelblue', alpha=0.7, edgecolor='black')
        ax.set_yticks(y_pos)
        ax.set_yticklabels([str(g)[:30] for g in top_genes.index], fontsize=8)
        ax.set_xlabel('Mean Expression (TPM)', fontweight='bold')
        ax.set_title(f'Top {n_top} Genes: {stage}\n(n={len(stage_samples)} samples)', 
                    fontweight='bold', fontsize=11)
        ax.invert_yaxis()
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add value labels
        for i, (gene, val) in enumerate(zip(top_genes.index, top_genes.values)):
            ax.text(val, i, f' {val:.0f}', va='center', fontsize=7)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_eda_top_genes_by_stage.png')
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_file}")

def plot_gene_expression_heatmap(expr_df, metadata, n_genes=50, output_dir='visualizations'):
    """
    Plot gene expression heatmap of top variable genes.
    
    Parameters
    ----------
    expr_df : DataFrame
        Expression matrix
    metadata : DataFrame
        Sample metadata
    n_genes : int
        Number of top variable genes to show
    output_dir : str
        Output directory
    """
    # Calculate coefficient of variation for each gene
    gene_cv = expr_df.std(axis=1) / (expr_df.mean(axis=1) + 1e-8)
    top_var_genes = gene_cv.nlargest(n_genes).index
    
    # Subset expression matrix
    heatmap_data = expr_df.loc[top_var_genes]
    
    # Order samples by stage if available
    if 'stage' in metadata.columns:
        stage_order = {'GV': 0, 'MI': 1, 'MII': 2}
        metadata_sorted = metadata.sort_values('stage', key=lambda x: x.map(stage_order))
        heatmap_data = heatmap_data[metadata_sorted.index]
        col_colors = metadata_sorted['stage'].map({'GV': '#F8766D', 'MI': '#39B600', 'MII': '#00B0F6'})
    else:
        col_colors = None
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Log transform for visualization (but keep raw data)
    heatmap_data_log = np.log1p(heatmap_data)
    
    sns.heatmap(heatmap_data_log, 
                cmap='viridis',
                cbar_kws={'label': 'Log(Expression + 1)'},
                xticklabels=[s[:20] for s in heatmap_data.columns],
                yticklabels=[str(g)[:30] for g in top_var_genes],
                ax=ax,
                linewidths=0.5,
                linecolor='gray')
    
    ax.set_xlabel('Samples', fontweight='bold', fontsize=11)
    ax.set_ylabel('Top Variable Genes', fontweight='bold', fontsize=11)
    ax.set_title(f'Gene Expression Heatmap\n(Top {n_genes} Most Variable Genes, Raw Data)', 
                fontweight='bold', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=7)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_eda_heatmap.png')
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_file}")

def plot_gene_correlation_matrix(expr_df, n_genes=30, output_dir='visualizations'):
    """
    Plot correlation matrix for top expressed genes.
    
    Parameters
    ----------
    expr_df : DataFrame
        Expression matrix
    n_genes : int
        Number of top genes to include
    output_dir : str
        Output directory
    """
    # Select top expressed genes
    mean_expr = expr_df.mean(axis=1).sort_values(ascending=False)
    top_genes = mean_expr.head(n_genes).index
    
    # Calculate correlation matrix
    corr_matrix = expr_df.loc[top_genes].T.corr()
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sns.heatmap(corr_matrix,
                cmap='coolwarm',
                center=0,
                vmin=-1, vmax=1,
                square=True,
                cbar_kws={'label': 'Pearson Correlation'},
                xticklabels=[str(g)[:20] for g in top_genes],
                yticklabels=[str(g)[:20] for g in top_genes],
                ax=ax,
                linewidths=0.5,
                linecolor='gray')
    
    ax.set_xlabel('Genes', fontweight='bold', fontsize=11)
    ax.set_ylabel('Genes', fontweight='bold', fontsize=11)
    ax.set_title(f'Gene-Gene Correlation Matrix\n(Top {n_genes} Most Expressed Genes, Raw Data)', 
                fontweight='bold', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=7)
    plt.yticks(rotation=0, fontsize=7)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_eda_correlation_matrix.png')
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_file}")

def plot_raw_counts_boxplots(expr_df, metadata, n_genes=20, output_dir='visualizations'):
    """
    Plot boxplots of raw counts for top expressed genes.
    
    Parameters
    ----------
    expr_df : DataFrame
        Expression matrix
    metadata : DataFrame
        Sample metadata
    n_genes : int
        Number of top genes to plot
    output_dir : str
        Output directory
    """
    # Select top expressed genes
    mean_expr = expr_df.mean(axis=1).sort_values(ascending=False)
    top_genes = mean_expr.head(n_genes).index
    
    # Prepare data for plotting
    plot_data = []
    for gene in top_genes:
        gene_expr = expr_df.loc[gene]
        if 'stage' in metadata.columns:
            for stage in metadata['stage'].dropna().unique():
                stage_samples = metadata[metadata['stage'] == stage].index
                stage_expr = gene_expr[stage_samples]
                for val in stage_expr.values:
                    plot_data.append({
                        'Gene': str(gene)[:30],
                        'Stage': stage,
                        'Expression': val
                    })
        else:
            for val in gene_expr.values:
                plot_data.append({
                    'Gene': str(gene)[:30],
                    'Stage': 'All',
                    'Expression': val
                })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create figure
    n_cols = 4
    n_rows = (n_genes + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    for idx, gene in enumerate(top_genes):
        ax = axes[idx]
        gene_data = plot_df[plot_df['Gene'] == str(gene)[:30]]
        
        if 'stage' in metadata.columns and len(gene_data['Stage'].unique()) > 1:
            stages = gene_data['Stage'].unique()
            stage_data = [gene_data[gene_data['Stage'] == s]['Expression'].values for s in stages]
            bp = ax.boxplot(stage_data, labels=stages, patch_artist=True)
            for patch in bp['boxes']:
                patch.set_facecolor('lightblue')
                patch.set_alpha(0.7)
        else:
            ax.boxplot([gene_data['Expression'].values], labels=['All'])
        
        ax.set_title(str(gene)[:25], fontsize=9, fontweight='bold')
        ax.set_ylabel('Expression (TPM)', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Hide extra subplots
    for idx in range(n_genes, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(f'Raw Expression Boxplots\n(Top {n_genes} Most Expressed Genes)', 
                fontweight='bold', fontsize=14, y=0.995)
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'raw_eda_boxplots.png')
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_file}")

def main():
    """Main function to generate all raw EDA plots."""
    
    print("="*70)
    print("Raw EDA Visualization Generation")
    print("="*70)
    
    # Create output directory
    output_dir = 'visualizations'
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("\n1. Loading raw expression data...")
    try:
        expr_df, metadata = load_raw_expression_data()
        print(f"   Loaded: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")
        print(f"   Expression range: {expr_df.min().min():.2f} - {expr_df.max().max():.2f}")
    except Exception as e:
        print(f"   ERROR: Could not load data: {e}")
        return
    
    # Generate plots
    print("\n2. Generating library size distribution...")
    try:
        plot_library_size_distribution(expr_df, metadata, output_dir)
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n3. Generating mitochondrial gene percentage...")
    try:
        plot_mitochondrial_percentage(expr_df, metadata, output_dir)
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n4. Generating top genes by stage...")
    try:
        plot_top_genes_by_stage(expr_df, metadata, n_top=20, output_dir=output_dir)
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n5. Generating gene expression heatmap...")
    try:
        plot_gene_expression_heatmap(expr_df, metadata, n_genes=50, output_dir=output_dir)
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n6. Generating gene correlation matrix...")
    try:
        plot_gene_correlation_matrix(expr_df, n_genes=30, output_dir=output_dir)
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n7. Generating raw counts boxplots...")
    try:
        plot_raw_counts_boxplots(expr_df, metadata, n_genes=20, output_dir=output_dir)
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n" + "="*70)
    print("Raw EDA visualization generation complete!")
    print(f"All plots saved to: {output_dir}/")
    print("="*70)

if __name__ == '__main__':
    main()

