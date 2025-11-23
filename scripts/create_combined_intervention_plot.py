#!/usr/bin/env python3
"""
Combined Intervention Plot: Pseudotime/Latent Age vs Health Score
- x-axis: pseudotime / latent age z
- y-axis: health score
- dot size/color: uncertainty
- horizontal lines: optimal/critical thresholds
- shaded region: optimal intervention window (high CHS + low σ)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import os

print("")
print("COMBINED INTERVENTION PLOT: Pseudotime/Latent Age vs Health Score")
print("")

# Load Data

clinical_csv = '../data/clinical_decision_framework_final.csv'
sample_csv = '../data/sample_metadata_with_age.csv'

if not os.path.exists(clinical_csv):
    print(" Error: clinical_decision_framework_final.csv not found")
    exit(1)

df = pd.read_csv(clinical_csv, index_col=0)
print(f" Loaded clinical data: {len(df)} samples")

# Load sample metadata for stage and potential pseudotime
if os.path.exists(sample_csv):
    sample_df = pd.read_csv(sample_csv)
    # Merge on sample name
    if 'sample' in sample_df.columns:
        df = df.merge(sample_df[['sample', 'stage']], 
                     left_index=True, right_on='sample', how='left')
        df.set_index('sample', inplace=True, drop=False)
        print(f" Merged stage information")
else:
    df['stage'] = 'Unknown'

# Compute/Estimate Health Score if Missing

if 'health_score' not in df.columns and 'oocyte_health_score' not in df.columns:
    print("\n Health score not found. Computing proxy from available data...")
    
    # Proxy health score based on:
    # 1. Inverse of risk_score (lower risk = higher health)
    # 2. Stage (GV = higher health than MI)
    # 3. Inverse of uncertainty (lower uncertainty = higher confidence = higher health)
    
    # Normalize risk_score (inverse)
    if 'risk_score' in df.columns:
        risk_norm = 1.0 - (df['risk_score'] - df['risk_score'].min()) / (df['risk_score'].max() - df['risk_score'].min() + 1e-8)
    else:
        risk_norm = np.ones(len(df)) * 0.5
    
    # Stage contribution (GV = 1.0, MI = 0.6, MII = 0.4)
    stage_map = {'GV': 1.0, 'MI': 0.6, 'MII': 0.4, 'Unknown': 0.5}
    stage_score = df['stage'].map(stage_map).fillna(0.5).values
    
    # Uncertainty contribution (inverse, normalized)
    if 'cellular_age_uncertainty' in df.columns:
        unc_norm = 1.0 - (df['cellular_age_uncertainty'] - df['cellular_age_uncertainty'].min()) / (df['cellular_age_uncertainty'].max() - df['cellular_age_uncertainty'].min() + 1e-8)
    else:
        unc_norm = np.ones(len(df)) * 0.5
    
    # Weighted combination
    health_score = (0.4 * risk_norm + 0.4 * stage_score + 0.2 * unc_norm) * 100
    
    df['health_score'] = health_score
    print(f" Computed proxy health score (range: {df['health_score'].min():.1f}-{df['health_score'].max():.1f})")
else:
    health_col = 'health_score' if 'health_score' in df.columns else 'oocyte_health_score'
    df['health_score'] = df[health_col]
    print(f" Using existing health score (range: {df['health_score'].min():.1f}-{df['health_score'].max():.1f})")

# Get X-axis: Pseudotime or Latent Age Z

# Try to find pseudotime first, fallback to cellular_age_z
if 'dpt_pseudotime' in df.columns:
    x_data = df['dpt_pseudotime'].values
    x_label = 'Diffusion Pseudotime (τ)'
    print(f" Using DPT pseudotime (range: {x_data.min():.3f}-{x_data.max():.3f})")
elif 'pseudotime' in df.columns:
    x_data = df['pseudotime'].values
    x_label = 'Pseudotime'
    print(f" Using pseudotime (range: {x_data.min():.3f}-{x_data.max():.3f})")
elif 'cellular_age_z' in df.columns:
    x_data = df['cellular_age_z'].values
    x_label = 'Cellular Age (Latent Z)'
    print(f" Using cellular age Z (range: {x_data.min():.3f}-{x_data.max():.3f})")
else:
    print(" Error: No pseudotime or cellular_age_z found")
    exit(1)

# Get Y-axis: Health Score

y_data = df['health_score'].values
y_label = 'Clinical Health Score (CHS)'

# Get Uncertainty for Size/Color

if 'cellular_age_uncertainty' in df.columns:
    uncertainty = df['cellular_age_uncertainty'].values
    print(f" Using uncertainty (range: {uncertainty.min():.2f}-{uncertainty.max():.2f})")
else:
    print(" Uncertainty not found. Using constant uncertainty.")
    uncertainty = np.ones(len(df)) * 0.5

# Define Thresholds

# Define thresholds based on percentiles
critical_threshold = np.percentile(y_data, 25)  # Bottom 25% = critical
warning_threshold = np.percentile(y_data, 50)   # Bottom 50% = warning
optimal_threshold = np.percentile(y_data, 75)   # Top 25% = optimal

print(f"\nThresholds:")
print(f"  Critical: {critical_threshold:.1f} (bottom 25%)")
print(f"  Warning: {warning_threshold:.1f} (median)")
print(f"  Optimal: {optimal_threshold:.1f} (top 25%)")

# Define Optimal Intervention Window

# Optimal window: High CHS (>= optimal_threshold) AND Low uncertainty (bottom 33%)
uncertainty_low_threshold = np.percentile(uncertainty, 33)  # Bottom 33% = low uncertainty

print(f"  Low uncertainty threshold: {uncertainty_low_threshold:.2f} (bottom 33%)")

# Create Combined Plot

print("\n[Creating combined intervention plot...]")

fig, ax = plt.subplots(figsize=(14, 10))

# Normalize uncertainty for size (min size = 50, max size = 500)
uncertainty_normalized = (uncertainty - uncertainty.min()) / (uncertainty.max() - uncertainty.min() + 1e-8)
dot_sizes = 50 + uncertainty_normalized * 450  # Range: 50-500

# Color by uncertainty (red = high uncertainty, blue = low uncertainty)
# Use a colormap that goes from blue (low) to red (high)
colors = plt.cm.RdYlBu_r(uncertainty_normalized)

# Shade optimal intervention window
# Define the region: high CHS (>= optimal_threshold) AND low uncertainty (<= uncertainty_low_threshold)
x_min, x_max = x_data.min(), x_data.max()
y_min, y_max = y_data.min(), y_data.max()

# Create a polygon for the optimal intervention window
optimal_x_range = [x_min, x_max, x_max, x_min]
optimal_y_range = [optimal_threshold, optimal_threshold, y_max, y_max]

# But we also want to constrain by uncertainty, so we'll shade the region
# where uncertainty is low AND health score is high
optimal_mask = (y_data >= optimal_threshold) & (uncertainty <= uncertainty_low_threshold)

# Shade the high CHS region (above optimal threshold)
ax.fill_between([x_min, x_max], optimal_threshold, y_max, 
                alpha=0.15, color='green', label='High CHS Region')

# Additionally shade points that are in optimal window (high CHS + low uncertainty)
if optimal_mask.sum() > 0:
    optimal_x = x_data[optimal_mask]
    optimal_y = y_data[optimal_mask]
    # Create a convex hull or simple bounding box
    if len(optimal_x) > 0:
        opt_x_min, opt_x_max = optimal_x.min(), optimal_x.max()
        opt_y_min, opt_y_max = optimal_y.min(), optimal_y.max()
        # Shade a tighter region
        ax.fill_between([opt_x_min, opt_x_max], opt_y_min, opt_y_max,
                        alpha=0.25, color='lime', 
                        label=f'Optimal Intervention Window\n(High CHS + Low σ, n={optimal_mask.sum()})')

# Plot scatter points
scatter = ax.scatter(x_data, y_data, s=dot_sizes, c=colors, 
                       alpha=0.7, edgecolors='black', linewidths=1.5,
                       zorder=5)

# Add colorbar for uncertainty
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.RdYlBu_r, 
                                         norm=plt.Normalize(vmin=uncertainty.min(), 
                                                           vmax=uncertainty.max())),
                   ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Uncertainty (σ)', fontsize=12, fontweight='bold')

# Add horizontal threshold lines
ax.axhline(y=critical_threshold, color='darkred', linestyle='--', 
          linewidth=2.5, label=f'Critical Threshold ({critical_threshold:.1f})', zorder=3)
ax.axhline(y=warning_threshold, color='orange', linestyle='--', 
          linewidth=2.5, label=f'Warning Threshold ({warning_threshold:.1f})', zorder=3)
ax.axhline(y=optimal_threshold, color='green', linestyle='--', 
          linewidth=2.5, label=f'Optimal Threshold ({optimal_threshold:.1f})', zorder=3)

# Add text annotations for thresholds
ax.text(x_max * 0.98, critical_threshold, 'Critical', 
       ha='right', va='bottom', fontsize=10, fontweight='bold',
       bbox=dict(boxstyle='round,pad=0.3', facecolor='red', alpha=0.3))
ax.text(x_max * 0.98, warning_threshold, 'Warning', 
       ha='right', va='bottom', fontsize=10, fontweight='bold',
       bbox=dict(boxstyle='round,pad=0.3', facecolor='orange', alpha=0.3))
ax.text(x_max * 0.98, optimal_threshold, 'Optimal', 
       ha='right', va='bottom', fontsize=10, fontweight='bold',
       bbox=dict(boxstyle='round,pad=0.3', facecolor='green', alpha=0.3))

# Labels and title
ax.set_xlabel(x_label, fontsize=14, fontweight='bold')
ax.set_ylabel(y_label, fontsize=14, fontweight='bold')
ax.set_title('Intervention Decision Framework:\nPseudotime/Latent Age vs Health Score with Uncertainty', 
            fontsize=16, fontweight='bold', pad=15)

# Add legend
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
              markersize=10, label='Low Uncertainty', markeredgecolor='black'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
              markersize=10, label='High Uncertainty', markeredgecolor='black'),
    mpatches.Patch(facecolor='green', alpha=0.15, label='High CHS Region'),
    mpatches.Patch(facecolor='lime', alpha=0.25, label='Optimal Window (High CHS + Low σ)')
]

# Combine with threshold lines
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, loc='upper left', fontsize=10, 
         framealpha=0.95, frameon=True, fancybox=True, shadow=True)

# Grid
ax.grid(True, alpha=0.3, linestyle='--', zorder=0)

# Add statistics text box
stats_text = f"n = {len(df)}\n"
stats_text += f"Optimal Window: {optimal_mask.sum()} cells ({optimal_mask.sum()/len(df)*100:.1f}%)\n"
stats_text += f"Mean CHS: {y_data.mean():.1f} ± {y_data.std():.1f}\n"
stats_text += f"Mean Uncertainty: {uncertainty.mean():.2f} ± {uncertainty.std():.2f}"

ax.text(0.02, 0.02, stats_text, transform=ax.transAxes,
       fontsize=10, verticalalignment='bottom',
       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Stage annotation if available
if 'stage' in df.columns:
    stage_counts = df['stage'].value_counts()
    stage_text = " | ".join([f"{stage}: {count}" for stage, count in stage_counts.items()])
    ax.text(0.5, 0.98, f"Stages: {stage_text}", transform=ax.transAxes,
           ha='center', va='top', fontsize=9, style='italic',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

plt.tight_layout()
plt.savefig('../visualizations/combined_intervention_plot.png', dpi=300, bbox_inches='tight',
           facecolor='white', edgecolor='none')
plt.close()

print(" Saved: combined_intervention_plot.png")

# Print summary
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"Total cells: {len(df)}")
print(f"Optimal intervention window (High CHS + Low σ): {optimal_mask.sum()} cells ({optimal_mask.sum()/len(df)*100:.1f}%)")
print(f"\nHealth Score Statistics:")
print(f"  Mean: {y_data.mean():.1f} ± {y_data.std():.1f}")
print(f"  Range: [{y_data.min():.1f}, {y_data.max():.1f}]")
print(f"\nUncertainty Statistics:")
print(f"  Mean: {uncertainty.mean():.2f} ± {uncertainty.std():.2f}")
print(f"  Range: [{uncertainty.min():.2f}, {uncertainty.max():.2f}]")
print(f"\nX-axis ({x_label}):")
print(f"  Range: [{x_data.min():.3f}, {x_data.max():.3f}]")

print(f"\n{'='*70}")
print(" Combined intervention plot complete!")
print(f"{'='*70}")

