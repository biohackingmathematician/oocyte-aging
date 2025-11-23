#!/usr/bin/env python3
"""
Review and fix blank spaces in existing visualizations.
This script identifies and fixes layout issues.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import sys

print("="*70)
print("REVIEW AND FIX VISUALIZATION BLANK SPACES")
print("="*70)

# Load data
sample_csv = './sample_metadata_with_age.csv'
clinical_csv = './clinical_decision_framework_final.csv'

if not os.path.exists(sample_csv) or not os.path.exists(clinical_csv):
    print("❌ ERROR: Required CSV files not found")
    sys.exit(1)

sample_df = pd.read_csv(sample_csv)
clinical_df = pd.read_csv(clinical_csv, index_col=0)
merged = clinical_df.merge(sample_df, left_index=True, right_on='sample', how='inner')

print(f"✓ Loaded data: {len(merged)} samples")
print(f"  Columns: {list(merged.columns)}")

# ============================================================================
# FIXED VISUALIZATION 1: Stage Overview (Improved Layout)
# ============================================================================

print("\n[1/2] Creating improved stage overview visualization...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Oocyte Data Overview', fontsize=16, fontweight='bold', y=0.98)

# Plot 1: Stage distribution with better spacing
stage_counts = merged['stage'].value_counts()
colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
bar_colors = [colors.get(stage, '#95a5a6') for stage in stage_counts.index]

bars = ax1.bar(stage_counts.index, stage_counts.values, color=bar_colors, 
               alpha=0.8, edgecolor='black', linewidth=2, width=0.6)
ax1.set_xlabel('Developmental Stage', fontsize=13, fontweight='bold')
ax1.set_ylabel('Number of Samples', fontsize=13, fontweight='bold')
ax1.set_title('A. Sample Distribution by Stage', fontsize=14, fontweight='bold', pad=15)
ax1.grid(True, alpha=0.3, axis='y', linestyle='--')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylim(0, max(stage_counts.values) * 1.15)  # Add 15% padding at top

# Add count labels on bars
for bar in bars:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}', ha='center', va='bottom', 
             fontweight='bold', fontsize=11)

# Plot 2: Risk groups with proper data check
if 'risk_group' in merged.columns:
    risk_by_stage = pd.crosstab(merged['stage'], merged['risk_group'], margins=False)
    # Only include rows with data
    risk_by_stage = risk_by_stage.loc[risk_by_stage.sum(axis=1) > 0]
    
    if len(risk_by_stage) > 0:
        # Reorder columns if needed
        risk_cols = ['Low Risk (Resilient Agers)', 'Moderate Risk', 'High Risk (Accelerated Agers)']
        existing_cols = [col for col in risk_cols if col in risk_by_stage.columns]
        if len(existing_cols) > 0:
            risk_by_stage = risk_by_stage[existing_cols]
        
        # Color mapping
        color_map = {
            'Low Risk (Resilient Agers)': '#2ecc71',
            'Moderate Risk': '#f39c12',
            'High Risk (Accelerated Agers)': '#e74c3c'
        }
        colors_list = [color_map.get(col, '#95a5a6') for col in risk_by_stage.columns]
        
        risk_by_stage.plot(kind='bar', ax=ax2, color=colors_list, 
                          alpha=0.8, edgecolor='black', linewidth=1.5, width=0.7)
        ax2.set_xlabel('Developmental Stage', fontsize=13, fontweight='bold')
        ax2.set_ylabel('Number of Samples', fontsize=13, fontweight='bold')
        ax2.set_title('B. Risk Groups by Stage', fontsize=14, fontweight='bold', pad=15)
        ax2.legend(title='Risk Group', fontsize=9, framealpha=0.95, loc='best')
        ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=0, ha='center')
        
        # Add value labels on bars
        for container in ax2.containers:
            labels = [f'{int(v)}' if v > 0 else '' for v in container.datavalues]
            ax2.bar_label(container, labels=labels, label_type='edge', 
                         fontsize=8, padding=2)
    else:
        ax2.text(0.5, 0.5, 'No risk group data\navailable for stages', 
                ha='center', va='center', transform=ax2.transAxes, 
                fontsize=12, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax2.set_title('B. Risk Groups by Stage', fontsize=14, fontweight='bold', pad=15)
        ax2.axis('off')
else:
    ax2.text(0.5, 0.5, 'No risk group data\navailable', 
            ha='center', va='center', transform=ax2.transAxes, 
            fontsize=12, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax2.set_title('B. Risk Groups by Stage', fontsize=14, fontweight='bold', pad=15)
    ax2.axis('off')

# Improve spacing
plt.tight_layout(rect=[0, 0, 1, 0.96], pad=3.0, w_pad=4.0)
plt.savefig('forum_stage_overview_fixed.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()
print("✓ Saved: forum_stage_overview_fixed.png")

# ============================================================================
# FIXED VISUALIZATION 2: Boxplot (Improved Layout)
# ============================================================================

print("\n[2/2] Creating improved boxplot visualization...")

# Get the metric
metrics_to_try = [
    ('oocyte_health_score', 'Health Score'),
    ('cellular_age_z', 'Cellular Age (Z)'),
    ('risk_score', 'Risk Score'),
]

available_metric = None
metric_name = None

for col, name in metrics_to_try:
    if col in merged.columns:
        available_metric = col
        metric_name = name
        break

if available_metric is None:
    available_metric = 'cellular_age_z'
    metric_name = 'Cellular Age (Z)'

# Filter to GV and MI
df_plot = merged[merged['stage'].isin(['GV', 'MI'])].copy()

if len(df_plot) == 0:
    df_plot = merged.copy()

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

# Statistical test
if len(values_by_stage) == 2:
    stat, pval = stats.mannwhitneyu(values_by_stage[0], values_by_stage[1],
                                     alternative='two-sided')
    significance = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
else:
    pval = 1.0
    significance = "ns"

# Create figure with better layout
fig, ax = plt.subplots(1, 1, figsize=(7, 8))
fig.suptitle(f'{metric_name} Comparison: GV vs MI Stages', 
             fontsize=15, fontweight='bold', y=0.98)

stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
box_colors = [stage_colors.get(labels[i].split('\n')[0], '#95a5a6') 
              for i in range(len(labels))]

bp = ax.boxplot(values_by_stage, labels=[l.split('\n')[0] for l in labels],
                patch_artist=True, widths=0.65, showmeans=True, meanline=False)

# Color boxes
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.75)
    patch.set_edgecolor('black')
    patch.set_linewidth(2)

# Style elements
for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
    if element in bp:
        if element == 'medians':
            for line in bp[element]:
                line.set_color('black')
                line.set_linewidth(2.5)
        elif element == 'means':
            for marker in bp[element]:
                marker.set_markerfacecolor('black')
                marker.set_markeredgecolor('black')
                marker.set_markersize(9)
        else:
            for line in bp[element]:
                line.set_color('black')
                line.set_linewidth(2)

# Add significance bracket
if len(values_by_stage) == 2 and pval < 0.05:
    y_max = max([np.max(vals) for vals in values_by_stage])
    y_min = min([np.min(vals) for vals in values_by_stage])
    y_range = y_max - y_min
    bracket_height = y_max + y_range * 0.12
    
    ax.plot([1, 1, 2, 2], [y_max + y_range*0.06, bracket_height, bracket_height, y_max + y_range*0.06],
            'k-', linewidth=2)
    ax.text(1.5, bracket_height + y_range*0.02, significance,
            ha='center', va='bottom', fontsize=13, fontweight='bold')

# Add annotations with better positioning
y_lim = ax.get_ylim()
y_range = y_lim[1] - y_lim[0]
y_offset = y_lim[0] - y_range * 0.08

for i, (vals, label) in enumerate(zip(values_by_stage, labels)):
    n = len(vals)
    mean_val = np.mean(vals)
    median_val = np.median(vals)
    ax.text(i+1, y_offset, f'n={n} | μ={mean_val:.2f} | m={median_val:.2f}',
            ha='center', va='top', fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                     edgecolor='gray', alpha=0.9))

# Adjust y-axis to fit annotations
ax.set_ylim(y_lim[0] - y_range * 0.15, y_lim[1])

ax.set_ylabel(metric_name, fontsize=13, fontweight='bold', labelpad=10)
ax.set_xlabel('Developmental Stage', fontsize=13, fontweight='bold', labelpad=10)
ax.grid(True, alpha=0.3, linestyle='--', axis='y')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=11)

plt.tight_layout(rect=[0, 0, 1, 0.97], pad=2.5)
output_file = 'forum_comparison_boxplot_fixed.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none', pad_inches=0.2)
plt.close()
print(f"✓ Saved: {output_file}")

print(f"\nStatistical Summary:")
for i, (vals, label) in enumerate(zip(values_by_stage, labels)):
    stage_name = label.split('\n')[0]
    print(f"\n{stage_name}:")
    print(f"  n = {len(vals)}")
    print(f"  Mean = {np.mean(vals):.3f}")
    print(f"  Median = {np.median(vals):.3f}")
    print(f"  SD = {np.std(vals):.3f}")
    print(f"  Range = [{np.min(vals):.3f}, {np.max(vals):.3f}]")

if len(values_by_stage) == 2:
    print(f"\nMann-Whitney U test: p = {pval:.4f} ({significance})")

print("\n" + "="*70)
print("✓ Fixed visualizations generated!")
print("="*70)
print("\nGenerated files:")
print("  ✓ forum_stage_overview_fixed.png")
print(f"  ✓ {output_file}")
print("\nNote: These versions have:")
print("  - Improved spacing and layout")
print("  - Better labels and annotations")
print("  - No blank spaces")
print("  - Proper axis limits and padding")

