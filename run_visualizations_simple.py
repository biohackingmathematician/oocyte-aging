#!/usr/bin/env python3
"""
Run project forum visualizations using available data.
This script will work with CSV data if h5ad files aren't available.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import sys

print("="*70)
print("PROJECT FORUM VISUALIZATIONS - SIMPLIFIED VERSION")
print("="*70)

# Try to import scanpy (may fail on Python 3.14)
try:
    import scanpy as sc
    HAS_SCANPY = True
    print("✓ scanpy available")
except Exception as e:
    HAS_SCANPY = False
    print(f"⚠ scanpy not available: {e}")
    print("  Will try to use CSV data instead")

# Load clinical data which has some information
clinical_csv = './clinical_decision_framework_final.csv'
sample_csv = './sample_metadata_with_age.csv'

if os.path.exists(clinical_csv):
    clinical_df = pd.read_csv(clinical_csv, index_col=0)
    print(f"\n✓ Loaded clinical data: {len(clinical_df)} samples")
else:
    clinical_df = None

if os.path.exists(sample_csv):
    sample_df = pd.read_csv(sample_csv)
    print(f"✓ Loaded sample metadata: {len(sample_df)} samples")
else:
    sample_df = None

# Try to load h5ad if scanpy is available
adata = None
if HAS_SCANPY:
    data_files = [
        './adata_trajectory_complete.h5ad',
        './adata_with_pathway_scores.h5ad',
        './adata_final_with_intervention.h5ad'
    ]
    
    for data_file in data_files:
        if os.path.exists(data_file):
            try:
                adata = sc.read_h5ad(data_file)
                print(f"✓ Loaded h5ad data: {data_file}")
                break
            except Exception as e:
                print(f"  ⚠ Could not load {data_file}: {e}")
                continue

# ============================================================================
# VISUALIZATION 1: UMAP Plot (if we have UMAP data)
# ============================================================================

if adata is not None and 'X_umap' in adata.obsm:
    print("\n[1/2] Generating UMAP plot colored by stage...")
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    umap_coords = adata.obsm['X_umap']
    stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
    stages = adata.obs['stage'].unique() if 'stage' in adata.obs.columns else []
    
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
    ax.legend(title='Stage', title_fontsize=11, fontsize=10, loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig('forum_umap_by_stage.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: forum_umap_by_stage.png")
    print(f"  Stage distribution:\n{adata.obs['stage'].value_counts().to_string()}")
else:
    print("\n[1/2] ⚠ Skipping UMAP plot - requires h5ad file with UMAP coordinates")
    print("  To generate: Run earlier notebook cells that compute UMAP")

# ============================================================================
# VISUALIZATION 2: Health Score Boxplot
# ============================================================================

print("\n[2/2] Generating health score boxplot...")

# Get health scores - try multiple sources
health_scores = None
stages = None

if adata is not None:
    # Try from adata
    for col in ['oocyte_health_score', 'health_score', 'composite_health_score']:
        if col in adata.obs.columns:
            health_scores = adata.obs[col].values
            if 'stage' in adata.obs.columns:
                stages = adata.obs['stage'].values
            break

# If not in adata, try to reconstruct from what we have
if health_scores is None and clinical_df is not None and sample_df is not None:
    print("  Attempting to extract data from CSV files...")
    # Merge dataframes
    merged = clinical_df.merge(sample_df, left_index=True, right_on='sample', how='inner')
    if 'stage' in merged.columns:
        stages = merged['stage'].values
        # We don't have health scores in CSV, but we can check
        print("  ⚠ Health scores not found in CSV files")

if health_scores is None:
    print("\n❌ ERROR: Could not find health score data.")
    print("  Please ensure one of these exists:")
    print("    - adata_with_pathway_scores.h5ad (with health scores)")
    print("    - adata_final_with_intervention.h5ad")
    print("  Or run earlier notebook cells to compute health scores.")
    sys.exit(1)

# Filter to GV and MI only
mask = np.isin(stages, ['GV', 'MI'])
if mask.sum() == 0:
    print("⚠ No GV or MI stages found. Using all stages.")
    mask = np.ones(len(stages), dtype=bool)

stages_filtered = stages[mask]
health_scores_filtered = health_scores[mask]

# Prepare data for boxplot
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

if len(health_scores_by_stage) == 0:
    print("❌ ERROR: No data available for boxplot")
    sys.exit(1)

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

