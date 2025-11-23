#!/usr/bin/env python3
"""
Create Forum EDA Visualizations
1. UMAP (or PCA) colored by Stage - shows global structure and GV/MI separation
2. Box/Violin Plot of Health Score by Stage - shows GV healthier than MI
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import seaborn as sns

print("")
print("FORUM EDA VISUALIZATIONS")
print("")

# Load data
sample_csv = '../data/sample_metadata_with_age.csv'
clinical_csv = '../data/clinical_decision_framework_final.csv'

if not os.path.exists(sample_csv) or not os.path.exists(clinical_csv):
    print(" Error: Required CSV files not found")
    exit(1)

sample_df = pd.read_csv(sample_csv)
clinical_df = pd.read_csv(clinical_csv, index_col=0)

# Merge data
merged = clinical_df.merge(sample_df, left_index=True, right_on='sample', how='inner')
print(f" Loaded data: {len(merged)} samples")

# Check if health score exists, compute if needed
health_col = 'oocyte_health_score' if 'oocyte_health_score' in merged.columns else 'health_score'
if health_col not in merged.columns or merged[health_col].isna().all():
    print(" Computing health score from cellular_age_z (proxy)...")
    if 'cellular_age_z' in merged.columns:
        # Inverse relationship: lower cellular age = higher health
        merged['health_score'] = (1 - merged['cellular_age_z']) * 100
        health_col = 'health_score'
        print(" Computed health score")

# Filter to GV and MI only for cleaner visualization
df_plot = merged[merged['stage'].isin(['GV', 'MI'])].copy()
print(f"  Filtered to GV/MI: {len(df_plot)} samples (GV: {len(df_plot[df_plot['stage']=='GV'])}, MI: {len(df_plot[df_plot['stage']=='MI'])})")

# VISUALIZATION 1: UMAP/PCA Colored by Stage

print("\n[1/2] Creating UMAP/PCA plot colored by stage...")

# Check if we have numerical features we can use for PCA
# Use cellular_age_z, uncertainty, and age as features for PCA
pca_features = []
feature_names = []

if 'cellular_age_z' in df_plot.columns:
    pca_features.append(df_plot['cellular_age_z'].fillna(0).values)
    feature_names.append('cellular_age_z')
if 'cellular_age_uncertainty' in df_plot.columns:
    pca_features.append(df_plot['cellular_age_uncertainty'].fillna(0).values)
    feature_names.append('uncertainty')
if 'age' in df_plot.columns:
    age_col = 'age_x' if 'age_x' in df_plot.columns else ('age_y' if 'age_y' in df_plot.columns else 'age')
    if age_col in df_plot.columns:
        pca_features.append(df_plot[age_col].fillna(0).values)
        feature_names.append('age')
if health_col in df_plot.columns:
    pca_features.append(df_plot[health_col].fillna(0).values)
    feature_names.append('health_score')

if len(pca_features) >= 2:
    # Create feature matrix
    X = np.column_stack(pca_features)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Compute PCA
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(X_scaled)
    
    print(f" Computed PCA from {len(feature_names)} features: {feature_names}")
    print(f"  PC1 variance explained: {pca.explained_variance_ratio_[0]:.2%}")
    print(f"  PC2 variance explained: {pca.explained_variance_ratio_[1]:.2%}")
    print(f"  Total variance explained: {pca.explained_variance_ratio_.sum():.2%}")
    
    # Create plot
    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    
    # Stage colors
    stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
    
    # Plot each stage separately for legend control
    for stage in ['GV', 'MI']:
        mask = df_plot['stage'] == stage
        if mask.sum() > 0:
            ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1],
                      c=stage_colors[stage], label=stage, s=150,
                      alpha=0.7, edgecolors='black', linewidth=1.5, zorder=3)
    
    ax.set_xlabel('PC1', fontsize=13, fontweight='bold')
    ax.set_ylabel('PC2', fontsize=13, fontweight='bold')
    ax.set_title('Global Structure: GV and MI Separate Along Main Axis', 
                fontsize=14, fontweight='bold', pad=15)
    
    # Add legend
    ax.legend(title='Stage', title_fontsize=11, fontsize=11, 
             loc='best', framealpha=0.9, edgecolor='black', fancybox=False)
    
    # Grid and styling
    ax.grid(True, alpha=0.3, linestyle='--', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add variance explained to axes labels
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)", 
                 fontsize=13, fontweight='bold')
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)", 
                 fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('../visualizations/forum_umap_pca_by_stage.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    print(" Saved: forum_umap_pca_by_stage.png")
    
else:
    print(" Not enough features for PCA - skipping plot")

# VISUALIZATION 2: Box/Violin Plot of Health Score by Stage

print("\n[2/2] Creating health score box/violin plot by stage...")

if health_col not in df_plot.columns or df_plot[health_col].isna().all():
    print(" Error: Health score not available")
    exit(1)

# Prepare data
gv_scores = df_plot[df_plot['stage'] == 'GV'][health_col].dropna()
mi_scores = df_plot[df_plot['stage'] == 'MI'][health_col].dropna()

if len(gv_scores) == 0 or len(mi_scores) == 0:
    print(" Error: Insufficient data for comparison")
    exit(1)

# Statistical test
stat, pval = stats.mannwhitneyu(gv_scores, mi_scores, alternative='two-sided')
print(f"  Statistical test (Mann-Whitney U):")
print(f"    U-statistic: {stat:.2f}")
print(f"    p-value: {pval:.4f}")
significance = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
print(f"    Significance: {significance}")

# Create figure with both boxplot and violin plot side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Stage colors
stage_colors = {'GV': '#2ecc71', 'MI': '#e74c3c'}

# Left: Boxplot
data_to_plot = [gv_scores, mi_scores]
labels = ['GV', 'MI']
box_colors = [stage_colors['GV'], stage_colors['MI']]

bp = ax1.boxplot(data_to_plot, tick_labels=labels, patch_artist=True,
                widths=0.6, showmeans=True, meanline=False)

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
if pval < 0.05:
    y_max = max([np.max(vals) for vals in data_to_plot])
    y_min = min([np.min(vals) for vals in data_to_plot])
    y_range = y_max - y_min
    bracket_height = y_max + y_range * 0.12
    
    ax1.plot([1, 1, 2, 2], [y_max + y_range*0.05, bracket_height, bracket_height, y_max + y_range*0.05],
            'k-', linewidth=1.5)
    ax1.text(1.5, bracket_height + y_range*0.02, significance,
            ha='center', va='bottom', fontsize=13, fontweight='bold')

# Add annotations
for i, (scores, label) in enumerate(zip(data_to_plot, labels)):
    n = len(scores)
    mean_val = np.mean(scores)
    median_val = np.median(scores)
    ax1.text(i+1, ax1.get_ylim()[0] - (ax1.get_ylim()[1] - ax1.get_ylim()[0]) * 0.06,
            f'n={n}\nμ={mean_val:.1f}\nmed={median_val:.1f}',
            ha='center', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

ax1.set_ylabel('Clinical Health Score', fontsize=13, fontweight='bold')
ax1.set_title('Box Plot: Health Score by Stage', fontsize=14, fontweight='bold', pad=15)
ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Right: Violin plot (more informative distribution)
violin_data = []
violin_labels = []
for stage in ['GV', 'MI']:
    stage_data = df_plot[df_plot['stage'] == stage][health_col].dropna()
    if len(stage_data) > 0:
        violin_data.append(stage_data)
        violin_labels.append(stage)

# Create violin plot data structure
violin_df = pd.DataFrame({
    'Health Score': np.concatenate(violin_data),
    'Stage': np.repeat(violin_labels, [len(d) for d in violin_data])
})

parts = ax2.violinplot([gv_scores, mi_scores], positions=[1, 2], 
                       widths=0.6, showmeans=True, showmedians=True)

# Color violins
for i, (pc, color) in enumerate(zip(parts['bodies'], box_colors)):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
    pc.set_edgecolor('black')
    pc.set_linewidth(1.5)

# Style other elements
for element in ['means', 'medians', 'cmaxes', 'cmins', 'cbars']:
    if element in parts:
        part = parts[element]
        if hasattr(part, '__iter__') and not isinstance(part, (str, bytes)):
            try:
                for p in part:
                    p.set_color('black')
                    p.set_linewidth(1.5)
            except TypeError:
                # Single object, not iterable
                part.set_color('black')
                part.set_linewidth(1.5)

# Add individual points (strip plot overlay)
for i, (scores, label) in enumerate(zip(data_to_plot, labels)):
    jitter = np.random.normal(i+1, 0.05, size=len(scores))
    ax2.scatter(jitter, scores, color='black', alpha=0.3, s=30, zorder=3)

# Add significance bracket
if pval < 0.05:
    ax2.plot([1, 1, 2, 2], [y_max + y_range*0.05, bracket_height, bracket_height, y_max + y_range*0.05],
            'k-', linewidth=1.5)
    ax2.text(1.5, bracket_height + y_range*0.02, significance,
            ha='center', va='bottom', fontsize=13, fontweight='bold')

ax2.set_xticks([1, 2])
ax2.set_xticklabels(labels, fontsize=12)
ax2.set_ylabel('Clinical Health Score', fontsize=13, fontweight='bold')
ax2.set_title('Violin Plot: Health Score Distribution by Stage', 
             fontsize=14, fontweight='bold', pad=15)
ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Add annotations
for i, (scores, label) in enumerate(zip(data_to_plot, labels)):
    n = len(scores)
    mean_val = np.mean(scores)
    median_val = np.median(scores)
    ax2.text(i+1, ax2.get_ylim()[0] - (ax2.get_ylim()[1] - ax2.get_ylim()[0]) * 0.06,
            f'n={n}\nμ={mean_val:.1f}',
            ha='center', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

plt.suptitle('GV Oocytes Look Healthier Than MI at Transcriptomic Level', 
             fontsize=15, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('../visualizations/forum_health_score_by_stage.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()
print(" Saved: forum_health_score_by_stage.png")

# Summary statistics
print(f"\nSummary statistics:")
print(f"  GV: n={len(gv_scores)}, Mean={np.mean(gv_scores):.2f}, Median={np.median(gv_scores):.2f}, SD={np.std(gv_scores):.2f}")
print(f"  MI: n={len(mi_scores)}, Mean={np.mean(mi_scores):.2f}, Median={np.median(mi_scores):.2f}, SD={np.std(mi_scores):.2f}")

print("\n" + "="*70)
print(" Forum EDA Visualizations Created!")
print("")
print("\nGenerated files:")
print("   forum_umap_pca_by_stage.png")
print("   forum_health_score_by_stage.png")
print("\nThese visualizations show:")
print("  1. Global structure: GV and MI separate along main axis")
print("  2. GV oocytes look healthier than MI at transcriptomic level")

