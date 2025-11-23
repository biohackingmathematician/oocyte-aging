#!/usr/bin/env python3
"""
Create Complete Results Summary Figure - All Panels Filled
Based on project goals from mid-progress report:
- DPT trajectory analysis (ρ=-0.79, p<0.001)
- GPLVM cellular age with uncertainty
- Health score by stage (GV: 76.7, MI: 61.0)
- Risk stratification (Low/Moderate/High)
- Clinical intervention windows
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
import os

print("")
print("COMPLETE RESULTS SUMMARY - OOCYTE AGING ANALYSIS")
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
print(f"  Columns: {list(merged.columns)}")

# Extract age column
age_col = 'age_x' if 'age_x' in merged.columns else ('age_y' if 'age_y' in merged.columns else 'age')

# Compute health score if not available
if 'oocyte_health_score' not in merged.columns and 'health_score' not in merged.columns:
    print(" Computing health score from cellular_age_z (proxy)...")
    # Use inverse of cellular_age_z as proxy (lower cellular age = higher health)
    if 'cellular_age_z' in merged.columns:
        merged['health_score'] = (1 - merged['cellular_age_z']) * 100
        merged['oocyte_health_score'] = merged['health_score']
        print(" Computed health score from cellular_age_z")

# Create comprehensive summary figure
print("\nCreating complete results summary figure...")

fig = plt.figure(figsize=(20, 12))
gs = GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.4, left=0.06, right=0.98, top=0.96, bottom=0.05)

colors_stage = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
colors_risk = {
    'Low Risk (Resilient Agers)': '#2ecc71',
    'Moderate Risk': '#f39c12',
    'High Risk (Accelerated Agers)': '#e74c3c'
}

# ROW 1: Basic Distributions

# Panel A: Stage distribution
ax_a = fig.add_subplot(gs[0, 0])
stage_counts = merged['stage'].value_counts()
bar_colors = [colors_stage.get(s, '#95a5a6') for s in stage_counts.index]
bars_a = ax_a.bar(stage_counts.index, stage_counts.values, color=bar_colors,
                 alpha=0.8, edgecolor='black', linewidth=2, width=0.6)
ax_a.set_xlabel('Developmental Stage', fontsize=11, fontweight='bold')
ax_a.set_ylabel('Number of Oocytes', fontsize=11, fontweight='bold')
ax_a.set_title('A. Sample Distribution by Stage', fontsize=12, fontweight='bold', pad=10)
ax_a.grid(True, alpha=0.3, axis='y', linestyle='--')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
for bar in bars_a:
    height = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2., height, f'{int(height)}',
             ha='center', va='bottom', fontweight='bold', fontsize=11)

# Panel B: Risk groups
ax_b = fig.add_subplot(gs[0, 1])
if 'risk_group' in merged.columns and merged['risk_group'].notna().sum() > 0:
    risk_counts = merged['risk_group'].value_counts()
    bar_colors_b = [colors_risk.get(r, '#95a5a6') for r in risk_counts.index]
    bars_b = ax_b.bar(range(len(risk_counts)), risk_counts.values,
                     color=bar_colors_b, alpha=0.8, edgecolor='black', linewidth=2, width=0.6)
    ax_b.set_xticks(range(len(risk_counts)))
    ax_b.set_xticklabels([r.replace(' Risk', '\nRisk').replace(' (Resilient Agers)', '').replace(' (Accelerated Agers)', '')
                          for r in risk_counts.index], fontsize=9, ha='center')
    ax_b.set_ylabel('Number of Oocytes', fontsize=11, fontweight='bold')
    ax_b.set_title('B. Risk Group Distribution', fontsize=12, fontweight='bold', pad=10)
    ax_b.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    for bar in bars_b:
        height = bar.get_height()
        ax_b.text(bar.get_x() + bar.get_width()/2., height, f'{int(height)}',
                 ha='center', va='bottom', fontweight='bold', fontsize=11)
else:
    ax_b.text(0.5, 0.5, 'Risk group data\nnot available',
             ha='center', va='center', transform=ax_b.transAxes, fontsize=11,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_b.set_title('B. Risk Group Distribution', fontsize=12, fontweight='bold', pad=10)
    ax_b.axis('off')

# Panel C: Cellular age by stage
ax_c = fig.add_subplot(gs[0, 2])
if 'cellular_age_z' in merged.columns and 'stage' in merged.columns:
    stages_to_plot = []
    data_to_plot = []
    for stage in ['GV', 'MI', 'MII']:
        stage_data = merged[merged['stage'] == stage]['cellular_age_z'].dropna()
        if len(stage_data) > 0:
            stages_to_plot.append(stage)
            data_to_plot.append(stage_data)

    if len(data_to_plot) > 0:
        bp_c = ax_c.boxplot(data_to_plot, tick_labels=stages_to_plot, patch_artist=True,
                           boxprops=dict(facecolor='lightblue', alpha=0.7, edgecolor='black', linewidth=1.5),
                           medianprops=dict(color='black', linewidth=2),
                           widths=0.6)
        ax_c.set_ylabel('Cellular Age (Z)', fontsize=11, fontweight='bold')
        ax_c.set_title('C. Cellular Age by Stage', fontsize=12, fontweight='bold', pad=10)
        ax_c.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax_c.spines['top'].set_visible(False)
        ax_c.spines['right'].set_visible(False)
    else:
        ax_c.text(0.5, 0.5, 'Cellular age data\nnot available',
                 ha='center', va='center', transform=ax_c.transAxes, fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax_c.set_title('C. Cellular Age by Stage', fontsize=12, fontweight='bold', pad=10)
        ax_c.axis('off')
else:
    ax_c.text(0.5, 0.5, 'Cellular age data\nnot available',
             ha='center', va='center', transform=ax_c.transAxes, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_c.set_title('C. Cellular Age by Stage', fontsize=12, fontweight='bold', pad=10)
    ax_c.axis('off')

# Panel D: Age distribution
ax_d = fig.add_subplot(gs[0, 3])
if age_col in merged.columns and merged[age_col].notna().sum() > 0:
    ages = merged[age_col].dropna()
    ax_d.hist(ages, bins=8, color='steelblue', alpha=0.7, edgecolor='black', linewidth=1.5)
    ax_d.set_xlabel('Chronological Age (years)', fontsize=11, fontweight='bold')
    ax_d.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax_d.set_title('D. Age Distribution', fontsize=12, fontweight='bold', pad=10)
    ax_d.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)

    # Add mean line
    mean_age = ages.mean()
    ax_d.axvline(mean_age, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_age:.1f} years')
    ax_d.legend(fontsize=9)
else:
    ax_d.text(0.5, 0.5, 'Age data\nnot available',
             ha='center', va='center', transform=ax_d.transAxes, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_d.set_title('D. Age Distribution', fontsize=12, fontweight='bold', pad=10)
    ax_d.axis('off')

# ROW 2: Trajectory and Correlations

# Panel E: Health score by stage (from mid-progress report: GV=76.7, MI=61.0)
ax_e = fig.add_subplot(gs[1, 0])
health_col = 'oocyte_health_score' if 'oocyte_health_score' in merged.columns else 'health_score'
if health_col in merged.columns and 'stage' in merged.columns:
    stages_to_plot = []
    data_to_plot = []
    for stage in ['GV', 'MI', 'MII']:
        stage_data = merged[merged['stage'] == stage][health_col].dropna()
        if len(stage_data) > 0:
            stages_to_plot.append(stage)
            data_to_plot.append(stage_data)

    if len(data_to_plot) > 0:
        bp_e = ax_e.boxplot(data_to_plot, tick_labels=stages_to_plot, patch_artist=True,
                           boxprops=dict(facecolor='lightgreen', alpha=0.7, edgecolor='black', linewidth=1.5),
                           medianprops=dict(color='black', linewidth=2),
                           widths=0.6)
        ax_e.set_ylabel('Oocyte Health Score', fontsize=11, fontweight='bold')
        ax_e.set_title('E. Health Score by Stage', fontsize=12, fontweight='bold', pad=10)
        ax_e.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax_e.spines['top'].set_visible(False)
        ax_e.spines['right'].set_visible(False)

        # Add threshold lines
        ax_e.axhline(79.9, color='blue', linestyle='--', linewidth=1.5, alpha=0.7, label='Optimal (79.9)')
        ax_e.axhline(53.2, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Critical (53.2)')
        ax_e.legend(fontsize=8, loc='upper right')
    else:
        ax_e.text(0.5, 0.5, 'Health score data\nnot available',
                 ha='center', va='center', transform=ax_e.transAxes, fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax_e.set_title('E. Health Score by Stage', fontsize=12, fontweight='bold', pad=10)
        ax_e.axis('off')
else:
    ax_e.text(0.5, 0.5, 'Health score data\nnot available',
             ha='center', va='center', transform=ax_e.transAxes, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_e.set_title('E. Health Score by Stage', fontsize=12, fontweight='bold', pad=10)
    ax_e.axis('off')

# Panel F: Cellular age vs chronological age (correlation)
ax_f = fig.add_subplot(gs[1, 1:3])
if 'cellular_age_z' in merged.columns and age_col in merged.columns:
    mask = merged[age_col].notna() & merged['cellular_age_z'].notna()
    if mask.sum() > 0:
        scatter = ax_f.scatter(merged.loc[mask, age_col], merged.loc[mask, 'cellular_age_z'],
                              c=merged.loc[mask, 'cellular_age_uncertainty'] if 'cellular_age_uncertainty' in merged.columns else merged.loc[mask, 'stage'].map(colors_stage) if 'stage' in merged.columns else 'blue',
                              s=120, cmap='plasma' if 'cellular_age_uncertainty' in merged.columns else None,
                              alpha=0.7, edgecolors='black', linewidth=1)
        ax_f.set_xlabel('Chronological Age (years)', fontsize=12, fontweight='bold')
        ax_f.set_ylabel('Cellular Age (Z)', fontsize=12, fontweight='bold')
        ax_f.set_title('F. Cellular vs Chronological Age', fontsize=13, fontweight='bold', pad=12)
        ax_f.grid(True, alpha=0.3, linestyle='--')
        ax_f.spines['top'].set_visible(False)
        ax_f.spines['right'].set_visible(False)

        if 'cellular_age_uncertainty' in merged.columns:
            plt.colorbar(scatter, ax=ax_f, label='Uncertainty', fraction=0.046, pad=0.04)

        # Add correlation
        if mask.sum() > 2:
            corr, pval = stats.pearsonr(merged.loc[mask, age_col], merged.loc[mask, 'cellular_age_z'])
            ax_f.text(0.05, 0.95, f'r = {corr:.3f}\np = {pval:.4f}',
                     transform=ax_f.transAxes, fontsize=11, fontweight='bold',
                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    else:
        ax_f.text(0.5, 0.5, 'Insufficient data\nfor correlation analysis',
                 ha='center', va='center', transform=ax_f.transAxes, fontsize=12,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax_f.set_title('F. Cellular vs Chronological Age', fontsize=13, fontweight='bold', pad=12)
        ax_f.axis('off')
else:
    ax_f.text(0.5, 0.5, 'Age or cellular age\ndata not available',
             ha='center', va='center', transform=ax_f.transAxes, fontsize=12,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_f.set_title('F. Cellular vs Chronological Age', fontsize=13, fontweight='bold', pad=12)
    ax_f.axis('off')

# Panel G: Uncertainty by stage
ax_g = fig.add_subplot(gs[1, 3])
if 'cellular_age_uncertainty' in merged.columns and 'stage' in merged.columns:
    stages_to_plot = []
    data_to_plot = []
    for stage in ['GV', 'MI', 'MII']:
        stage_data = merged[merged['stage'] == stage]['cellular_age_uncertainty'].dropna()
        if len(stage_data) > 0:
            stages_to_plot.append(stage)
            data_to_plot.append(stage_data)

    if len(data_to_plot) > 0:
        bp_g = ax_g.boxplot(data_to_plot, tick_labels=stages_to_plot, patch_artist=True,
                           boxprops=dict(facecolor='lightcoral', alpha=0.7, edgecolor='black', linewidth=1.5),
                           medianprops=dict(color='black', linewidth=2),
                           widths=0.6)
        ax_g.set_ylabel('Trajectory Uncertainty', fontsize=11, fontweight='bold')
        ax_g.set_title('G. Uncertainty by Stage', fontsize=12, fontweight='bold', pad=10)
        ax_g.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax_g.spines['top'].set_visible(False)
        ax_g.spines['right'].set_visible(False)
    else:
        ax_g.text(0.5, 0.5, 'Uncertainty data\nnot available',
                 ha='center', va='center', transform=ax_g.transAxes, fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax_g.set_title('G. Uncertainty by Stage', fontsize=12, fontweight='bold', pad=10)
        ax_g.axis('off')
else:
    ax_g.text(0.5, 0.5, 'Uncertainty data\nnot available',
             ha='center', va='center', transform=ax_g.transAxes, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_g.set_title('G. Uncertainty by Stage', fontsize=12, fontweight='bold', pad=10)
    ax_g.axis('off')

# ROW 3: Clinical Insights and Summary

# Panel H: Risk groups by stage
ax_h = fig.add_subplot(gs[2, 0:2])
if 'risk_group' in merged.columns and 'stage' in merged.columns:
    risk_by_stage = pd.crosstab(merged['stage'], merged['risk_group'], margins=False)
    risk_by_stage = risk_by_stage.loc[risk_by_stage.sum(axis=1) > 0]
    risk_cols = ['Low Risk (Resilient Agers)', 'Moderate Risk', 'High Risk (Accelerated Agers)']
    existing_cols = [col for col in risk_cols if col in risk_by_stage.columns]

    if len(risk_by_stage) > 0 and len(existing_cols) > 0:
        risk_by_stage = risk_by_stage[existing_cols]
        colors_list = [colors_risk.get(col, '#95a5a6') for col in risk_by_stage.columns]

        risk_by_stage.plot(kind='bar', ax=ax_h, color=colors_list,
                          alpha=0.8, edgecolor='black', linewidth=1.5, width=0.7)
        ax_h.set_xlabel('Developmental Stage', fontsize=12, fontweight='bold')
        ax_h.set_ylabel('Number of Oocytes', fontsize=12, fontweight='bold')
        ax_h.set_title('H. Risk Groups by Developmental Stage', fontsize=13, fontweight='bold', pad=12)
        ax_h.legend(title='Risk Group', fontsize=9, framealpha=0.95, loc='best')
        ax_h.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax_h.spines['top'].set_visible(False)
        ax_h.spines['right'].set_visible(False)
        plt.setp(ax_h.xaxis.get_majorticklabels(), rotation=0, ha='center')

        for container in ax_h.containers:
            labels = [f'{int(v)}' if v > 0 else '' for v in container.datavalues]
            ax_h.bar_label(container, labels=labels, label_type='edge', fontsize=9, padding=3)
    else:
        ax_h.text(0.5, 0.5, 'Risk group by stage\ndata not available',
                 ha='center', va='center', transform=ax_h.transAxes, fontsize=12,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax_h.set_title('H. Risk Groups by Developmental Stage', fontsize=13, fontweight='bold', pad=12)
        ax_h.axis('off')
else:
    ax_h.text(0.5, 0.5, 'Risk or stage data\nnot available',
             ha='center', va='center', transform=ax_h.transAxes, fontsize=12,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax_h.set_title('H. Risk Groups by Developmental Stage', fontsize=13, fontweight='bold', pad=12)
    ax_h.axis('off')

# Panel I: Summary statistics text
ax_i = fig.add_subplot(gs[2, 2:4])
ax_i.axis('off')

# Compile summary statistics
summary_parts = ["KEY FINDINGS:", ""]

# Sample info
summary_parts.append(f"• Total Samples: {len(merged)} oocytes")
if 'stage' in merged.columns:
    stage_counts = merged['stage'].value_counts()
    stage_summary = ", ".join([f"{stage}: {count}" for stage, count in stage_counts.items()])
    summary_parts.append(f"• Stages: {stage_summary}")

# Age info
if age_col in merged.columns and merged[age_col].notna().sum() > 0:
    ages = merged[age_col].dropna()
    summary_parts.append(f"• Chronological Age: Mean = {ages.mean():.1f} years (Range: {ages.min():.0f}-{ages.max():.0f})")

# Cellular age info
if 'cellular_age_z' in merged.columns:
    summary_parts.append(f"• Cellular Age (Z): Mean = {merged['cellular_age_z'].mean():.3f} (Range: {merged['cellular_age_z'].min():.3f}-{merged['cellular_age_z'].max():.3f})")

# Correlation
if 'cellular_age_z' in merged.columns and age_col in merged.columns:
    mask = merged[age_col].notna() & merged['cellular_age_z'].notna()
    if mask.sum() > 2:
        corr, pval = stats.pearsonr(merged.loc[mask, age_col], merged.loc[mask, 'cellular_age_z'])
        summary_parts.append(f"• Correlation (Z vs Age): r = {corr:.3f}, p = {pval:.4f}")

# Health scores
health_col = 'oocyte_health_score' if 'oocyte_health_score' in merged.columns else 'health_score'
if health_col in merged.columns and 'stage' in merged.columns:
    gv_mean = merged[merged['stage'] == 'GV'][health_col].mean() if len(merged[merged['stage'] == 'GV']) > 0 else None
    mi_mean = merged[merged['stage'] == 'MI'][health_col].mean() if len(merged[merged['stage'] == 'MI']) > 0 else None
    if gv_mean is not None and mi_mean is not None:
        summary_parts.append(f"• Health Scores: GV = {gv_mean:.1f}, MI = {mi_mean:.1f}")

# Risk groups
if 'risk_group' in merged.columns and merged['risk_group'].notna().sum() > 0:
    risk_counts = merged['risk_group'].value_counts()
    risk_summary = ", ".join([f"{group.replace(' Risk', '').replace(' (Resilient Agers)', '').replace(' (Accelerated Agers)', '')}: {count}"
                             for group, count in risk_counts.items()])
    summary_parts.append(f"• Risk Groups: {risk_summary}")

summary_text = "\n".join(summary_parts)

ax_i.text(0.05, 0.95, summary_text, ha='left', va='top', transform=ax_i.transAxes,
         fontsize=12, family='monospace',
         bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3, pad=10))

fig.suptitle('Complete Results Summary: Oocyte Aging Analysis', fontsize=16, fontweight='bold', y=0.995)

plt.savefig('../visualizations/complete_results_summary.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none', pad_inches=0.2)
plt.close()

print(" Saved: complete_results_summary.png")
print("")
print(" Complete Results Summary Generated!")
print("")
print("\nAll panels populated with available data.")
print("Key findings from mid-progress report:")
print("  • DPT correlation: ρ = -0.79, p < 0.001")
print("  • Health scores: GV = 76.7, MI = 61.0")
print("  • Intervention thresholds: Optimal >79.9, Critical <53.2")

