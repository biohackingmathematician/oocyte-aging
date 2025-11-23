#!/usr/bin/env python3
"""
Generate visualizations from available CSV data.
Works locally without requiring h5ad files or scanpy.
This creates simplified versions based on available data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

print("="*70)
print("PROJECT FORUM VISUALIZATIONS - From CSV Data")
print("="*70)

# Load available CSV files
sample_csv = '../data/sample_metadata_with_age.csv'
clinical_csv = '../data/clinical_decision_framework_final.csv'

if not os.path.exists(sample_csv):
    print(f"❌ ERROR: {sample_csv} not found")
    exit(1)

if not os.path.exists(clinical_csv):
    print(f"❌ ERROR: {clinical_csv} not found")
    exit(1)

# Load data
sample_df = pd.read_csv(sample_csv)
clinical_df = pd.read_csv(clinical_csv, index_col=0)

print(f"✓ Loaded sample metadata: {len(sample_df)} samples")
print(f"✓ Loaded clinical data: {len(clinical_df)} samples")

# Merge data
merged = clinical_df.merge(sample_df, left_index=True, right_on='sample', how='inner')
print(f"✓ Merged data: {len(merged)} samples")

# Check what we have
print(f"\nAvailable columns: {list(merged.columns)}")
print(f"\nStage distribution:")
print(merged['stage'].value_counts())

# ============================================================================
# VISUALIZATION 1: Simplified Stage Distribution
# ============================================================================

print("\n[1/2] Creating stage distribution visualization...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Bar chart of stages
stage_counts = merged['stage'].value_counts()
colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
bar_colors = [colors.get(stage, '#95a5a6') for stage in stage_counts.index]

bars = ax1.bar(stage_counts.index, stage_counts.values, color=bar_colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax1.set_xlabel('Developmental Stage', fontsize=12, fontweight='bold')
ax1.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
ax1.set_title('Sample Distribution by Stage', fontsize=14, fontweight='bold', pad=15)
ax1.grid(True, alpha=0.3, axis='y')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add count labels on bars
for bar in bars:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}', ha='center', va='bottom', fontweight='bold')

# Plot 2: Risk group by stage
if 'risk_group' in merged.columns:
    risk_by_stage = pd.crosstab(merged['stage'], merged['risk_group'])
    # Only plot stages that exist
    risk_by_stage = risk_by_stage.loc[risk_by_stage.sum(axis=1) > 0]
    if len(risk_by_stage) > 0:
        risk_by_stage.plot(kind='bar', ax=ax2, 
                          color=['#2ecc71', '#f39c12', '#e74c3c'], 
                          alpha=0.7, edgecolor='black', width=0.8)
        ax2.set_xlabel('Developmental Stage', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
        ax2.set_title('Risk Groups by Stage', fontsize=14, fontweight='bold', pad=15)
        ax2.legend(title='Risk Group', fontsize=9, framealpha=0.9)
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=0)
        # Add value labels on bars
        for container in ax2.containers:
            ax2.bar_label(container, fmt='%d', label_type='edge', fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'No risk group data available', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Risk Groups by Stage', fontsize=14, fontweight='bold', pad=15)
        ax2.axis('off')
else:
    # Age distribution by stage
    age_col = 'age_x' if 'age_x' in merged.columns else ('age_y' if 'age_y' in merged.columns else 'age')
    if age_col in merged.columns:
        gv_ages = merged[merged['stage'] == 'GV'][age_col].dropna()
        mi_ages = merged[merged['stage'] == 'MI'][age_col].dropna()
        
        if len(gv_ages) > 0 and len(mi_ages) > 0:
            bp2 = ax2.boxplot([gv_ages, mi_ages], labels=['GV', 'MI'], patch_artist=True,
                             boxprops=dict(facecolor='lightblue', alpha=0.7, edgecolor='black', linewidth=1.5),
                             medianprops=dict(color='black', linewidth=2))
            ax2.set_ylabel('Age (years)', fontsize=12, fontweight='bold')
            ax2.set_title('Age Distribution by Stage', fontsize=14, fontweight='bold', pad=15)
            ax2.grid(True, alpha=0.3, axis='y')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
        else:
            ax2.text(0.5, 0.5, 'Insufficient age data', 
                    ha='center', va='center', transform=ax2.transAxes, fontsize=12)
            ax2.set_title('Age Distribution by Stage', fontsize=14, fontweight='bold', pad=15)
            ax2.axis('off')
    else:
        ax2.text(0.5, 0.5, 'No age data available', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Additional Information', fontsize=14, fontweight='bold', pad=15)
        ax2.axis('off')

# Adjust spacing to avoid blank spaces
plt.tight_layout(pad=2.0, w_pad=3.0, h_pad=2.0)
plt.savefig('../visualizations/forum_stage_overview.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("✓ Saved: forum_stage_overview.png")

# ============================================================================
# VISUALIZATION 2: Boxplot Comparison (using available metrics)
# ============================================================================

print("\n[2/2] Creating comparison boxplot...")

# Try different metrics we might have
metrics_to_try = [
    ('oocyte_health_score', 'Health Score'),
    ('cellular_age_z', 'Cellular Age (Z)'),
    ('risk_score', 'Risk Score'),
    ('cellular_age_uncertainty', 'Uncertainty')
]

available_metric = None
metric_name = None

for col, name in metrics_to_try:
    if col in merged.columns:
        available_metric = col
        metric_name = name
        break

if available_metric is None:
    print("⚠ No health score or comparable metric found")
    print(f"  Available numeric columns: {merged.select_dtypes(include=[np.number]).columns.tolist()}")
    print("  Creating age-based comparison instead...")
    
    # Use age as fallback
    available_metric = 'age'
    metric_name = 'Chronological Age (years)'

# Filter to GV and MI only
df_plot = merged[merged['stage'].isin(['GV', 'MI'])].copy()

if len(df_plot) == 0:
    print("⚠ No GV or MI samples found")
    df_plot = merged.copy()

# Prepare data
gv_values = df_plot[df_plot['stage'] == 'GV'][available_metric].dropna().values
mi_values = df_plot[df_plot['stage'] == 'MI'][available_metric].dropna().values

values_by_stage = []
labels = []

if len(gv_values) > 0:
    values_by_stage.append(gv_values)
    labels.append(f'GV\n(n={len(gv_values)})')

if len(mi_values) > 0:
    values_by_stage.append(mi_values)
    labels.append(f'MI\n(n={len(mi_values)})')

if len(values_by_stage) == 0:
    print("❌ ERROR: No data available for boxplot")
    exit(1)

# Statistical test
if len(values_by_stage) == 2:
    stat, pval = stats.mannwhitneyu(values_by_stage[0], 
                                     values_by_stage[1],
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

bp = ax.boxplot(values_by_stage, 
                tick_labels=[l.split('\n')[0] for l in labels],
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
if len(values_by_stage) == 2 and pval < 0.05:
    y_max = max([np.max(vals) for vals in values_by_stage])
    y_min = min([np.min(vals) for vals in values_by_stage])
    y_range = y_max - y_min
    bracket_height = y_max + y_range * 0.1
    ax.plot([1, 1, 2, 2], [y_max + y_range*0.05, bracket_height, bracket_height, y_max + y_range*0.05],
            'k-', linewidth=1.5)
    ax.text(1.5, bracket_height + y_range*0.02, significance,
            ha='center', va='bottom', fontsize=12, fontweight='bold')

# Add annotations
for i, (vals, label) in enumerate(zip(values_by_stage, labels)):
    n = len(vals)
    mean_val = np.mean(vals)
    ax.text(i+1, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05,
            f'n={n}\nμ={mean_val:.2f}', ha='center', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

ax.set_ylabel(metric_name, fontsize=12, fontweight='bold')
title = f'{metric_name} Comparison: GV vs MI Stages'
ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Determine output filename
output_file = 'forum_health_score_boxplot.png' if 'health' in metric_name.lower() else f'forum_{available_metric}_boxplot.png'

# Adjust spacing and ensure no blank areas
plt.tight_layout(pad=2.0)
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.1)
plt.close()
print(f"✓ Saved: {output_file}")

print(f"\nSummary statistics for {metric_name}:")
for i, (vals, label) in enumerate(zip(values_by_stage, labels)):
    stage_name = label.split('\n')[0]
    print(f"\n{stage_name}:")
    print(f"  n = {len(vals)}")
    print(f"  Mean = {np.mean(vals):.2f}")
    print(f"  Median = {np.median(vals):.2f}")
    print(f"  SD = {np.std(vals):.2f}")
    print(f"  Range = [{np.min(vals):.2f}, {np.max(vals):.2f}]")

print("\n" + "="*70)
print("✓ Visualizations generated from CSV data!")
print("="*70)
print("\nGenerated files:")
if os.path.exists('forum_stage_overview.png'):
    print("  ✓ forum_stage_overview.png")
if os.path.exists(output_file):
    print(f"  ✓ {output_file}")
print("\nNote: UMAP plot requires h5ad files. Run notebook cells first to generate them.")

