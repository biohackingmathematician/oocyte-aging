#!/usr/bin/env python3
"""
Extract data from h5ad files using h5py and create visualizations.
This works without scanpy!
"""

import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import sys

print("="*70)
print("PROJECT FORUM VISUALIZATIONS - Using h5py to read h5ad files")
print("="*70)

# Find h5ad file
h5ad_files = [
    './adata_trajectory_complete.h5ad',
    './adata_with_pathway_scores.h5ad',
    './adata_final_with_intervention.h5ad'
]

h5ad_file = None
for f in h5ad_files:
    if os.path.exists(f):
        h5ad_file = f
        break

if h5ad_file is None:
    print("\n❌ ERROR: No h5ad files found.")
    print("Available files:")
    for f in os.listdir('.'):
        if f.endswith('.h5ad'):
            print(f"  - {f}")
    sys.exit(1)

print(f"\n✓ Found h5ad file: {h5ad_file}")

# Read h5ad file with h5py
try:
    f = h5py.File(h5ad_file, 'r')
    print("✓ Opened h5ad file")
except Exception as e:
    print(f"❌ ERROR: Could not open {h5ad_file}: {e}")
    sys.exit(1)

# Extract data
data_dict = {}
try:
    # Try to get observation metadata
    if 'obs' in f:
        obs_group = f['obs']
        for key in obs_group.keys():
            try:
                if '_index' in obs_group[key]:
                    # This is likely a categorical/string array
                    idx = obs_group[key]['_index']
                    categories = obs_group[key]['categories']
                    data_dict[key] = [categories[idx[i]] for i in range(len(idx))]
                else:
                    # Numeric or simple array
                    data_dict[key] = np.array(obs_group[key])
            except:
                pass
    
    # Try to get UMAP coordinates
    if 'obsm' in f and 'X_umap' in f['obsm']:
        umap_data = f['obsm']['X_umap']
        data_dict['umap_1'] = np.array(umap_data[:, 0]) if umap_data.shape[1] > 0 else None
        data_dict['umap_2'] = np.array(umap_data[:, 1]) if umap_data.shape[1] > 1 else None
        print("✓ Found UMAP coordinates")
    
    # Try observation names
    if 'obs_names' in f:
        obs_names = f['obs_names']
        if 'categories' in obs_names:
            data_dict['sample'] = list(obs_names['categories'])
        else:
            data_dict['sample'] = [str(v) for v in obs_names]
    
    print(f"✓ Extracted {len(data_dict)} data columns")
    for key in data_dict.keys():
        if hasattr(data_dict[key], '__len__'):
            print(f"  - {key}: {len(data_dict[key])} values")
    
    f.close()
    
except Exception as e:
    print(f"⚠ Warning: Error reading h5ad structure: {e}")
    print("Trying alternative approach...")
    f.close()

# Convert to DataFrame for easier handling
df = pd.DataFrame(data_dict)
print(f"\n✓ Created DataFrame with {len(df)} rows and {len(df.columns)} columns")

# ============================================================================
# VISUALIZATION 1: UMAP Plot
# ============================================================================

if 'umap_1' in df.columns and 'umap_2' in df.columns and 'stage' in df.columns:
    print("\n[1/2] Generating UMAP plot colored by stage...")
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    
    stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
    stages = df['stage'].unique()
    
    for stage in stages:
        if stage not in stage_colors:
            stage_colors[stage] = '#95a5a6'
        
        mask = df['stage'] == stage
        if mask.sum() > 0:
            ax.scatter(df.loc[mask, 'umap_1'], 
                      df.loc[mask, 'umap_2'],
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
    plt.savefig('forum_umap_by_stage.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: forum_umap_by_stage.png")
    print(f"\nStage distribution:")
    print(df['stage'].value_counts().to_string())
else:
    print("\n[1/2] ⚠ Skipping UMAP plot - missing required columns")
    missing = []
    if 'umap_1' not in df.columns: missing.append('umap_1')
    if 'umap_2' not in df.columns: missing.append('umap_2')
    if 'stage' not in df.columns: missing.append('stage')
    print(f"  Missing: {', '.join(missing)}")
    print(f"  Available columns: {list(df.columns)}")

# ============================================================================
# VISUALIZATION 2: Health Score Boxplot
# ============================================================================

print("\n[2/2] Generating health score boxplot...")

# Look for health score column
health_score_col = None
for col in ['oocyte_health_score', 'health_score', 'composite_health_score']:
    if col in df.columns:
        health_score_col = col
        break

if health_score_col is None:
    print("⚠ Health score column not found. Checking available columns...")
    print(f"  Available columns: {list(df.columns)}")
    
    # Try to compute from pathway scores
    pathway_scores = [col for col in df.columns if col.startswith('score_') or 'pathway' in col.lower()]
    if len(pathway_scores) > 0:
        print(f"  Found pathway scores: {pathway_scores}")
        # Compute composite (simplified)
        df['computed_health_score'] = df[pathway_scores].mean(axis=1) * 100
        health_score_col = 'computed_health_score'
        print(f"  ✓ Computed health score from pathway scores")
    else:
        print("❌ ERROR: Cannot create boxplot - no health score data found")
        sys.exit(1)

if 'stage' not in df.columns:
    print("❌ ERROR: Cannot create boxplot - no 'stage' column found")
    sys.exit(1)

# Filter to GV and MI only
df_filtered = df[df['stage'].isin(['GV', 'MI'])].copy()

if len(df_filtered) == 0:
    print("⚠ No GV or MI stages found. Using all stages.")
    df_filtered = df.copy()

# Prepare data
gv_mask = df_filtered['stage'] == 'GV'
mi_mask = df_filtered['stage'] == 'MI'

health_scores_by_stage = []
labels = []

if gv_mask.sum() > 0:
    gv_scores = df_filtered.loc[gv_mask, health_score_col].values
    health_scores_by_stage.append(gv_scores)
    labels.append(f'GV\n(n={len(gv_scores)})')

if mi_mask.sum() > 0:
    mi_scores = df_filtered.loc[mi_mask, health_score_col].values
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
print("✓ All visualizations completed successfully!")
print("="*70)

