#!/usr/bin/env python3
"""
PROJECT FORUM VISUALIZATIONS
Generate UMAP plot and health score boxplot for project forum presentation.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import sys

print("="*70)
print("PROJECT FORUM VISUALIZATIONS")
print("="*70)

# ============================================================================
# VISUALIZATION 1: UMAP Plot Colored by Stage (GV vs MI)
# ============================================================================

print("\n[1/2] Generating UMAP plot colored by stage...")

# Try to load data
adata = None
data_files = [
    './adata_trajectory_complete.h5ad',
    './adata_with_pathway_scores.h5ad',
    './adata_final_with_intervention.h5ad'
]

for data_file in data_files:
    if os.path.exists(data_file):
        try:
            adata = sc.read_h5ad(data_file)
            print(f"✓ Loaded data: {data_file}")
            print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")
            break
        except Exception as e:
            print(f"  ⚠ Could not load {data_file}: {e}")
            continue

if adata is None:
    print("\n❌ ERROR: Could not load any data file.")
    print("Available files in directory:")
    for f in os.listdir('.'):
        if f.endswith('.h5ad') or f.endswith('.csv'):
            print(f"  - {f}")
    sys.exit(1)

# Check if UMAP exists
if 'X_umap' not in adata.obsm:
    print("\nComputing UMAP embeddings...")
    if 'X_pca' not in adata.obsm:
        print("  Computing PCA...")
        sc.pp.pca(adata, n_comps=10)
    print("  Computing neighbors and UMAP...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
    sc.tl.umap(adata)
    print("  ✓ UMAP computed")

# Check if stage column exists
if 'stage' not in adata.obs.columns:
    print(f"\n⚠ WARNING: 'stage' column not found in data.")
    print(f"Available obs columns: {list(adata.obs.columns)[:10]}...")
    print("Cannot create stage-colored plot without stage information.")
    sys.exit(1)

# Create UMAP figure
fig, ax = plt.subplots(1, 1, figsize=(8, 7))

umap_coords = adata.obsm['X_umap']
stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
stages = adata.obs['stage'].unique()

for stage in stages:
    if stage not in stage_colors:
        stage_colors[stage] = '#95a5a6'
    
    mask = adata.obs['stage'] == stage
    if mask.sum() > 0:
        ax.scatter(umap_coords[mask, 0], 
                  umap_coords[mask, 1],
                  c=stage_colors[stage],
                  label=stage,
                  s=150,
                  alpha=0.7,
                  edgecolors='black',
                  linewidth=1.5)

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
output_file1 = 'forum_umap_by_stage.png'
plt.savefig(output_file1, dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_file1}")

print(f"\nStage distribution:")
print(adata.obs['stage'].value_counts().to_string())

# ============================================================================
# VISUALIZATION 2: Side-by-Side Boxplot - Health Score GV vs MI
# ============================================================================

print("\n[2/2] Generating health score boxplot...")

# Check if health score exists
health_score_col = None
for col in ['oocyte_health_score', 'health_score', 'composite_health_score']:
    if col in adata.obs.columns:
        health_score_col = col
        break

if health_score_col is None:
    print("\nComputing health score from pathway scores...")
    pathway_scores = [col for col in adata.obs.columns if col.startswith('score_')]
    if len(pathway_scores) > 0:
        weights = {'Cell_Cycle': 0.20, 'Mitochondrial_OXPHOS': 0.30, 
                  'DNA_Damage': 0.15, 'Spindle_Assembly': 0.20, 
                  'Oocyte_Quality': 0.15}
        composite_score = np.zeros(adata.n_obs)
        for pathway, weight in weights.items():
            score_col = f'score_{pathway}'
            if score_col in adata.obs.columns:
                composite_score += weight * adata.obs[score_col].values
        adata.obs['oocyte_health_score'] = (
            (composite_score - composite_score.min()) / 
            (composite_score.max() - composite_score.min() + 1e-8)
        ) * 100
        health_score_col = 'oocyte_health_score'
        print(f"✓ Computed health score (range: {adata.obs[health_score_col].min():.1f}-{adata.obs[health_score_col].max():.1f})")
    else:
        print("❌ ERROR: Cannot compute health score - no pathway scores found")
        sys.exit(1)

# Filter to GV and MI only
mask = adata.obs['stage'].isin(['GV', 'MI'])
stages_to_plot = ['GV', 'MI']

# Prepare data for boxplot
health_scores_by_stage = []
labels = []

for stage in stages_to_plot:
    stage_mask = adata.obs['stage'] == stage
    scores = adata.obs.loc[stage_mask, health_score_col].values
    if len(scores) > 0:
        health_scores_by_stage.append(scores)
        labels.append(f'{stage}\n(n={len(scores)})')

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

# Create boxplot
bp = ax.boxplot(health_scores_by_stage, 
                labels=[l.split('\n')[0] for l in labels],
                patch_artist=True,
                widths=0.6,
                showmeans=True,
                meanline=False)

# Color the boxes
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.5)

# Style the plot elements
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

# Add significance bracket if applicable
if len(health_scores_by_stage) == 2 and pval < 0.05:
    y_max = max([np.max(scores) for scores in health_scores_by_stage])
    y_min = min([np.min(scores) for scores in health_scores_by_stage])
    y_range = y_max - y_min
    bracket_height = y_max + y_range * 0.1
    
    ax.plot([1, 1, 2, 2], [y_max + y_range*0.05, bracket_height, bracket_height, y_max + y_range*0.05],
            'k-', linewidth=1.5)
    ax.text(1.5, bracket_height + y_range*0.02, significance,
            ha='center', va='bottom', fontsize=12, fontweight='bold')

# Add sample size annotations
for i, (scores, label) in enumerate(zip(health_scores_by_stage, labels)):
    n = len(scores)
    mean_val = np.mean(scores)
    ax.text(i+1, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05,
            f'n={n}\nμ={mean_val:.1f}',
            ha='center', va='top', fontsize=9, 
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

ax.set_ylabel('Oocyte Health Score', fontsize=12, fontweight='bold')
ax.set_title('Health Score Comparison: GV vs MI Stages', 
             fontsize=14, fontweight='bold', pad=15)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
output_file2 = 'forum_health_score_boxplot.png'
plt.savefig(output_file2, dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_file2}")

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
print("✓ All visualizations completed successfully!")
print("="*70)

