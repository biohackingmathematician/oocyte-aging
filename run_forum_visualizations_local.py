#!/usr/bin/env python3
"""
Standalone script to generate project forum visualizations locally.
Works without Jupyter - just run: python3 run_forum_visualizations_local.py
"""

import sys
import os

# Check Python version
print("="*70)
print("PROJECT FORUM VISUALIZATIONS - Local Execution")
print("="*70)
print(f"Python version: {sys.version}")

# Try to import required packages
try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    print("✓ Basic packages loaded")
except ImportError as e:
    print(f"❌ ERROR: Missing basic package: {e}")
    sys.exit(1)

# Try scanpy (may fail on Python 3.14)
try:
    import scanpy as sc
    HAS_SCANPY = True
    print("✓ scanpy available")
except ImportError as e:
    HAS_SCANPY = False
    print(f"⚠ scanpy not available: {e}")
    print("  Will try alternative methods...")

print("\n" + "="*70)

# ============================================================================
# Find and load data
# ============================================================================

def find_h5ad_file():
    """Find available h5ad file"""
    h5ad_files = [
        './adata_trajectory_complete.h5ad',
        './adata_with_pathway_scores.h5ad',
        './adata_final_with_intervention.h5ad',
        './adata_final_with_gplvm_and_risk.h5ad'
    ]
    
    for f in h5ad_files:
        if os.path.exists(f):
            return f
    return None

def load_data_scanpy(filename):
    """Load data using scanpy"""
    if not HAS_SCANPY:
        return None
    try:
        adata = sc.read_h5ad(filename)
        return adata
    except Exception as e:
        print(f"  ⚠ Error loading with scanpy: {e}")
        return None

def load_data_pandas():
    """Load data from CSV files as fallback"""
    data = {}
    
    # Load sample metadata
    if os.path.exists('./sample_metadata_with_age.csv'):
        sample_df = pd.read_csv('./sample_metadata_with_age.csv')
        data['samples'] = sample_df
        print(f"✓ Loaded sample metadata: {len(sample_df)} samples")
    
    # Load clinical data
    if os.path.exists('./clinical_decision_framework_final.csv'):
        clinical_df = pd.read_csv('./clinical_decision_framework_final.csv', index_col=0)
        data['clinical'] = clinical_df
        print(f"✓ Loaded clinical data: {len(clinical_df)} samples")
    
    return data if data else None

# Try to load h5ad file
h5ad_file = find_h5ad_file()
adata = None

if h5ad_file:
    print(f"\nFound h5ad file: {h5ad_file}")
    if HAS_SCANPY:
        print("Loading with scanpy...")
        adata = load_data_scanpy(h5ad_file)
        if adata:
            print(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    else:
        print("⚠ scanpy not available - cannot read h5ad file")
        print("  Falling back to CSV data...")
        csv_data = load_data_pandas()
else:
    print("\n⚠ No h5ad file found. Trying CSV files...")
    csv_data = load_data_pandas()

# ============================================================================
# VISUALIZATION 1: UMAP Plot Colored by Stage
# ============================================================================

if adata is not None:
    print("\n[1/2] Generating UMAP plot colored by stage...")
    
    # Check if UMAP exists
    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP embeddings...")
        if 'X_pca' not in adata.obsm:
            print("    Computing PCA first...")
            sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.umap(adata)
        print("  ✓ UMAP computed")
    
    if 'stage' in adata.obs.columns and 'X_umap' in adata.obsm:
        fig, ax = plt.subplots(1, 1, figsize=(8, 7))
        
        umap_coords = adata.obsm['X_umap']
        stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
        stages = adata.obs['stage'].unique()
        
        for stage in stages:
            if stage not in stage_colors:
                stage_colors[stage] = '#95a5a6'
            
            mask = adata.obs['stage'] == stage
            if mask.sum() > 0:
                ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                          c=stage_colors[stage], label=stage, s=150, alpha=0.7,
                          edgecolors='black', linewidth=1.5)
        
        ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
        ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
        ax.set_title('UMAP Visualization: Oocyte Developmental Stages', 
                     fontsize=14, fontweight='bold', pad=15)
        ax.legend(title='Stage', title_fontsize=11, fontsize=10, 
                  loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig('forum_umap_by_stage.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: forum_umap_by_stage.png")
        print(f"\nStage distribution:")
        print(adata.obs['stage'].value_counts().to_string())
    else:
        missing = []
        if 'X_umap' not in adata.obsm: missing.append('X_umap')
        if 'stage' not in adata.obs.columns: missing.append('stage')
        print(f"⚠ Skipping UMAP plot - missing: {', '.join(missing)}")
else:
    print("\n[1/2] ⚠ Skipping UMAP plot - requires h5ad file with scanpy")
    print("  To generate: Need to run notebook cells first to create h5ad files")

# ============================================================================
# VISUALIZATION 2: Health Score Boxplot
# ============================================================================

print("\n[2/2] Generating health score boxplot...")

health_scores = None
stages = None

# Try to get health scores from adata
if adata is not None:
    for col in ['oocyte_health_score', 'health_score', 'composite_health_score']:
        if col in adata.obs.columns:
            health_scores = adata.obs[col].values
            if 'stage' in adata.obs.columns:
                stages = adata.obs['stage'].values
            print(f"✓ Found health scores in column: {col}")
            break

# If not in adata, try CSV data
if health_scores is None:
    if 'csv_data' in locals() and csv_data:
        print("  Trying to extract from CSV files...")
        if 'samples' in csv_data and 'clinical' in csv_data:
            # Merge data
            merged = csv_data['clinical'].merge(
                csv_data['samples'], 
                left_index=True, 
                right_on='sample', 
                how='inner'
            )
            if 'stage' in merged.columns:
                stages = merged['stage'].values
                print(f"  ✓ Found stages for {len(stages)} samples")
            # Check if we can reconstruct health scores somehow
            print("  ⚠ Health scores not directly available in CSV")

# Try to compute health scores from pathway scores if available
if health_scores is None and adata is not None:
    pathway_scores = [col for col in adata.obs.columns 
                      if col.startswith('score_') or 'pathway' in col.lower()]
    if len(pathway_scores) > 0:
        print(f"  Computing health score from {len(pathway_scores)} pathway scores...")
        weights = {
            'Cell_Cycle': 0.20, 
            'Mitochondrial_OXPHOS': 0.30, 
            'DNA_Damage': 0.15, 
            'Spindle_Assembly': 0.20, 
            'Oocyte_Quality': 0.15
        }
        composite_score = np.zeros(adata.n_obs)
        for pathway, weight in weights.items():
            score_col = f'score_{pathway}'
            if score_col in adata.obs.columns:
                composite_score += weight * adata.obs[score_col].values
        
        adata.obs['oocyte_health_score'] = (
            (composite_score - composite_score.min()) / 
            (composite_score.max() - composite_score.min() + 1e-8)
        ) * 100
        health_scores = adata.obs['oocyte_health_score'].values
        if 'stage' in adata.obs.columns:
            stages = adata.obs['stage'].values
        print(f"  ✓ Computed health score (range: {health_scores.min():.1f}-{health_scores.max():.1f})")

if health_scores is None or stages is None:
    print("\n❌ ERROR: Cannot create boxplot - missing required data")
    print("  Required: health scores and stage information")
    print("\n  Solutions:")
    print("    1. Run notebook cells first to generate h5ad files")
    print("    2. Ensure h5ad file has 'oocyte_health_score' and 'stage' columns")
    sys.exit(1)

# Filter to GV and MI only
mask = np.isin(stages, ['GV', 'MI'])
stages_filtered = stages[mask]
health_scores_filtered = health_scores[mask]

# Prepare data
gv_scores = health_scores_filtered[stages_filtered == 'GV']
mi_scores = health_scores_filtered[stages_filtered == 'MI']

health_scores_by_stage = []
labels = []

if len(gv_scores) > 0:
    health_scores_by_stage.append(gv_scores)
    labels.append(f'GV\n(n={len(gv_scores)})')

if len(mi_scores) > 0:
    health_scores_by_stage.append(mi_scores)
    labels.append(f'MI\n(n={len(mi_scores)})')

# Statistical test
if len(health_scores_by_stage) == 2:
    stat, pval = stats.mannwhitneyu(health_scores_by_stage[0], 
                                     health_scores_by_stage[1],
                                     alternative='two-sided')
    print(f"\nStatistical test (Mann-Whitney U):")
    print(f"  U-statistic: {stat:.2f}")
    print(f"  p-value: {pval:.4f}")
    significance = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"  Significance: {significance}")
else:
    pval = 1.0
    significance = "ns"

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(6, 7))

stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
box_colors = [stage_colors.get(labels[i].split('\n')[0], '#95a5a6') for i in range(len(labels))]

bp = ax.boxplot(health_scores_by_stage, 
                labels=[l.split('\n')[0] for l in labels],
                patch_artist=True, widths=0.6, showmeans=True, meanline=False)

# Color boxes
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.5)

# Style elements
for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
    if element in bp:
        if element == 'medians':
            for line in bp[element]:
                line.set_color('black')
                line.set_linewidth(2)
        elif element == 'means':
            for marker in bp[element]:
                marker.set_markerfacecolor('black')
                marker.set_markeredgecolor('black')
                marker.set_markersize(8)
        else:
            for line in bp[element]:
                line.set_color('black')
                line.set_linewidth(1.5)

# Add significance bracket
if len(health_scores_by_stage) == 2 and pval < 0.05:
    y_max = max([np.max(scores) for scores in health_scores_by_stage])
    y_min = min([np.min(scores) for scores in health_scores_by_stage])
    y_range = y_max - y_min
    bracket_height = y_max + y_range * 0.1
    ax.plot([1, 1, 2, 2], [y_max + y_range*0.05, bracket_height, bracket_height, y_max + y_range*0.05],
            'k-', linewidth=1.5)
    ax.text(1.5, bracket_height + y_range*0.02, significance,
            ha='center', va='bottom', fontsize=12, fontweight='bold')

# Add annotations
for i, (scores, label) in enumerate(zip(health_scores_by_stage, labels)):
    n = len(scores)
    mean_val = np.mean(scores)
    ax.text(i+1, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05,
            f'n={n}\nμ={mean_val:.1f}', ha='center', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

ax.set_ylabel('Oocyte Health Score', fontsize=12, fontweight='bold')
ax.set_title('Health Score Comparison: GV vs MI Stages', 
             fontsize=14, fontweight='bold', pad=15)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('forum_health_score_boxplot.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: forum_health_score_boxplot.png")

print(f"\nSummary statistics:")
for i, (scores, label) in enumerate(zip(health_scores_by_stage, labels)):
    stage_name = label.split('\n')[0]
    print(f"\n{stage_name}:")
    print(f"  n = {len(scores)}")
    print(f"  Mean = {np.mean(scores):.2f}")
    print(f"  Median = {np.median(scores):.2f}")
    print(f"  SD = {np.std(scores):.2f}")
    print(f"  Range = [{np.min(scores):.2f}, {np.max(scores):.2f}]")

print("\n" + "="*70)
print("✓ Visualization generation complete!")
print("="*70)
print("\nGenerated files:")
if os.path.exists('forum_umap_by_stage.png'):
    print("  ✓ forum_umap_by_stage.png")
if os.path.exists('forum_health_score_boxplot.png'):
    print("  ✓ forum_health_score_boxplot.png")

