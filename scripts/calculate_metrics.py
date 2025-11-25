#!/usr/bin/env python3
"""
Calculate Additional Metrics for Model Validation
Implements Priority 1 metrics that can be calculated from existing data
"""

import pandas as pd
import numpy as np
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve, auc,
    silhouette_score, davies_bouldin_score, calinski_harabasz_score,
    confusion_matrix, classification_report
)
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import os

print("")
print("CALCULATING ADDITIONAL VALIDATION METRICS")
print("")

# Data loading

clinical_csv = '../data/clinical_decision_framework_final.csv'
sample_csv = '../data/sample_metadata_with_age.csv'

if not os.path.exists(clinical_csv):
    print("Error: clinical_decision_framework_final.csv not found")
    exit(1)

df = pd.read_csv(clinical_csv, index_col=0)
print(f"Loaded clinical data: {len(df)} samples")

# Sample metadata loading
if os.path.exists(sample_csv):
    sample_df = pd.read_csv(sample_csv)
    if 'sample' in sample_df.columns:
        df = df.merge(sample_df[['sample', 'stage']],
                     left_index=True, right_on='sample', how='left')
        df.set_index('sample', inplace=True, drop=False)
        print(f"Merged stage information")
else:
    df['stage'] = 'Unknown'

# METRIC 1: Risk Stratification Performance (AUC-ROC)

print("\n[1] Calculating Risk Stratification Performance...")

if 'risk_group' in df.columns and 'risk_score' in df.columns:
    # Binary label encoding for risk stratification: High Risk = 1, Others = 0
    df['is_high_risk'] = (df['risk_group'] == 'High Risk (Accelerated Agers)').astype(int)

    if df['is_high_risk'].sum() > 0:
        # Area under ROC curve computation
        try:
            auc_roc = roc_auc_score(df['is_high_risk'], df['risk_score'])
            print(" AUC-ROC for High Risk prediction: {auc_roc:.3f}")

            # ROC curve calculation
            fpr, tpr, thresholds = roc_curve(df['is_high_risk'], df['risk_score'])

            # ROC curve figure generation
            plt.figure(figsize=(8, 6))
            plt.plot(fpr, tpr, linewidth=2, label=f'ROC curve (AUC = {auc_roc:.3f})')
            plt.plot([0, 1], [0, 1], 'k--', label='Random classifier')
            plt.xlabel('False Positive Rate', fontsize=12, fontweight='bold')
            plt.ylabel('True Positive Rate', fontsize=12, fontweight='bold')
            plt.title('Risk Stratification: ROC Curve\n(High Risk vs Low/Moderate Risk)',
                     fontsize=14, fontweight='bold')
            plt.legend(fontsize=11)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig('../visualizations/metrics_roc_curve.png', dpi=300, bbox_inches='tight')
            plt.close()
            print(" Saved: ../visualizations/metrics_roc_curve.png")

            # Precision-recall curve calculation
            precision, recall, pr_thresholds = precision_recall_curve(df['is_high_risk'], df['risk_score'])
            pr_auc = auc(recall, precision)
            print(" Precision-Recall AUC: {pr_auc:.3f}")

            plt.figure(figsize=(8, 6))
            plt.plot(recall, precision, linewidth=2, label=f'PR curve (AUC = {pr_auc:.3f})')
            plt.xlabel('Recall', fontsize=12, fontweight='bold')
            plt.ylabel('Precision', fontsize=12, fontweight='bold')
            plt.title('Risk Stratification: Precision-Recall Curve',
                     fontsize=14, fontweight='bold')
            plt.legend(fontsize=11)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig('../visualizations/metrics_pr_curve.png', dpi=300, bbox_inches='tight')
            plt.close()
            print(" Saved: metrics_pr_curve.png")

        except Exception as e:
            print(" Could not calculate AUC-ROC: {e}")
    else:
        print(" No high-risk samples found for AUC calculation")
else:
    print(" Risk group data not available")

# METRIC 2: Clinical Health Score Discriminative Ability

print("\n[2] Calculating Clinical Health Score Validation...")

# Health score computation health score if missing
if 'health_score' not in df.columns and 'oocyte_health_score' not in df.columns:
    if 'risk_score' in df.columns:
        risk_norm = 1.0 - (df['risk_score'] - df['risk_score'].min()) / (df['risk_score'].max() - df['risk_score'].min() + 1e-8)
        stage_map = {'GV': 1.0, 'MI': 0.6, 'MII': 0.4, 'Unknown': 0.5}
        stage_score = df['stage'].map(stage_map).fillna(0.5).values
        if 'cellular_age_uncertainty' in df.columns:
            unc_norm = 1.0 - (df['cellular_age_uncertainty'] - df['cellular_age_uncertainty'].min()) / (df['cellular_age_uncertainty'].max() - df['cellular_age_uncertainty'].min() + 1e-8)
        else:
            unc_norm = np.ones(len(df)) * 0.5
        df['health_score'] = (0.4 * risk_norm + 0.4 * stage_score + 0.2 * unc_norm) * 100
    else:
        print(" Cannot compute health score proxy")
        df['health_score'] = np.nan

if 'health_score' in df.columns and 'stage' in df.columns:
    # Discriminative ability: GV vs MI
    gv_mi_mask = df['stage'].isin(['GV', 'MI'])
    if gv_mi_mask.sum() >= 4:  # Need at least 2 per group
        gv_mi_df = df[gv_mi_mask].copy()
        gv_mi_df['is_GV'] = (gv_mi_df['stage'] == 'GV').astype(int)

        try:
            chs_auc = roc_auc_score(gv_mi_df['is_GV'], gv_mi_df['health_score'])
            print(" CHS AUC for GV vs MI classification: {chs_auc:.3f}")

            # Non-parametric statistical testing
            gv_scores = gv_mi_df[gv_mi_df['stage'] == 'GV']['health_score'].values
            mi_scores = gv_mi_df[gv_mi_df['stage'] == 'MI']['health_score'].values

            from scipy.stats import mannwhitneyu
            stat, pval = mannwhitneyu(gv_scores, mi_scores, alternative='two-sided')
            print(" Mann-Whitney U test: U = {stat:.1f}, p = {pval:.4f}")

            # Effect size estimation (Cohen's d)
            mean_diff = gv_scores.mean() - mi_scores.mean()
            pooled_std = np.sqrt((gv_scores.std()**2 + mi_scores.std()**2) / 2)
            cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
            print(" Cohen's d effect size: {cohens_d:.3f}")

        except Exception as e:
            print(" Could not calculate CHS discriminative ability: {e}")
    else:
        print(" Insufficient GV/MI samples for classification")

# METRIC 3: Latent Space Quality

print("\n[3] Calculating Latent Space Quality Metrics...")

if 'cellular_age_z' in df.columns and 'stage' in df.columns:
    # Prepare data for clustering metrics
    stage_mask = df['stage'].isin(['GV', 'MI'])
    if stage_mask.sum() >= 4:
        X = df.loc[stage_mask, ['cellular_age_z']].values
        labels = df.loc[stage_mask, 'stage'].values

        # Convert labels to numeric
        label_map = {'GV': 0, 'MI': 1}
        numeric_labels = np.array([label_map.get(l, -1) for l in labels])
        valid_mask = numeric_labels >= 0

        if valid_mask.sum() >= 4:
            X_valid = X[valid_mask]
            labels_valid = numeric_labels[valid_mask]

            try:
                # Silhouette score
                silhouette = silhouette_score(X_valid, labels_valid)
                print(" Silhouette score (GV vs MI separation): {silhouette:.3f}")

                # Davies-Bouldin index
                db_index = davies_bouldin_score(X_valid, labels_valid)
                print(" Davies-Bouldin index: {db_index:.3f} (lower is better)")

                # Calinski-Harabasz score
                ch_score = calinski_harabasz_score(X_valid, labels_valid)
                print(" Calinski-Harabasz score: {ch_score:.3f} (higher is better)")

            except Exception as e:
                print(" Could not calculate clustering metrics: {e}")
        else:
            print(" Insufficient valid labels")
    else:
        print(" Insufficient stage data")
else:
    print(" Missing cellular_age_z or stage data")

# METRIC 4: Correlation Analysis

print("\n[4] Calculating Correlation Metrics...")

if 'cellular_age_z' in df.columns:
    # Correlation with chronological age
    if 'age' in df.columns:
        age_mask = df['age'].notna() & df['cellular_age_z'].notna()
        if age_mask.sum() >= 3:
            corr_age, pval_age = pearsonr(df.loc[age_mask, 'age'],
                                         df.loc[age_mask, 'cellular_age_z'])
            print(" Correlation (Z vs Age): r = {corr_age:.3f}, p = {pval_age:.4f}")

    # Correlation with health score
    if 'health_score' in df.columns:
        hs_mask = df['health_score'].notna() & df['cellular_age_z'].notna()
        if hs_mask.sum() >= 3:
            corr_hs, pval_hs = pearsonr(df.loc[hs_mask, 'cellular_age_z'],
                                       df.loc[hs_mask, 'health_score'])
            print(" Correlation (Z vs Health Score): r = {corr_hs:.3f}, p = {pval_hs:.4f}")

# METRIC 5: Age Discrepancy Analysis

print("\n[5] Calculating Age Discrepancy Metrics...")

if 'cellular_age_z' in df.columns and 'age' in df.columns:
    age_mask = df['age'].notna() & df['cellular_age_z'].notna()
    if age_mask.sum() >= 3:
        # Normalize both to same scale for comparison
        age_norm = (df.loc[age_mask, 'age'] - df.loc[age_mask, 'age'].min()) / \
                   (df.loc[age_mask, 'age'].max() - df.loc[age_mask, 'age'].min() + 1e-8)
        z_norm = df.loc[age_mask, 'cellular_age_z']

        # Age discrepancy: |Z_cellular - age_chronological| / age_chronological
        age_discrepancy = np.abs(z_norm - age_norm) / (age_norm + 1e-8)
        print(" Mean age discrepancy: {age_discrepancy.mean():.3f}")
        print(" Median age discrepancy: {age_discrepancy.median():.3f}")

        # Proportion of accelerated agers (Z >> chronological age)
        accelerated_threshold = 0.2  # 20% higher
        accelerated_mask = z_norm > (age_norm + accelerated_threshold)
        accelerated_pct = accelerated_mask.sum() / len(age_mask) * 100
        print(" Proportion of accelerated agers (Z > age + 0.2): {accelerated_pct:.1f}%")

# METRIC 6: Uncertainty Statistics

print("\n[6] Calculating Uncertainty Statistics...")

if 'cellular_age_uncertainty' in df.columns:
    uncertainty = df['cellular_age_uncertainty'].dropna()
    if len(uncertainty) > 0:
        print(" Mean uncertainty: {uncertainty.mean():.2f} ± {uncertainty.std():.2f}")
        print(" Uncertainty range: [{uncertainty.min():.2f}, {uncertainty.max():.2f}]")
        print(" Coefficient of variation: {uncertainty.std() / uncertainty.mean():.3f}")

        # Uncertainty by stage
        if 'stage' in df.columns:
            for stage in df['stage'].unique():
                stage_unc = df[df['stage'] == stage]['cellular_age_uncertainty'].dropna()
                if len(stage_unc) > 0:
                    print("{stage}: {stage_unc.mean():.2f} ± {stage_unc.std():.2f}")

# Generate Summary Report

print("")
print("METRICS SUMMARY")
print("")

metrics_summary = {
    'Risk Stratification': {
        'AUC-ROC': f"{auc_roc:.3f}" if 'auc_roc' in locals() else "N/A",
        'PR-AUC': f"{pr_auc:.3f}" if 'pr_auc' in locals() else "N/A",
        'Target': "> 0.75",
        'Status': "✅" if 'auc_roc' in locals() and auc_roc > 0.75 else "⏳"
    },
    'CHS Discriminative Ability': {
        'AUC (GV vs MI)': f"{chs_auc:.3f}" if 'chs_auc' in locals() else "N/A",
        "Cohen's d": f"{cohens_d:.3f}" if 'cohens_d' in locals() else "N/A",
        'Target': "> 0.7",
        'Status': "✅" if 'chs_auc' in locals() and chs_auc > 0.7 else "⏳"
    },
    'Latent Space Quality': {
        'Silhouette Score': f"{silhouette:.3f}" if 'silhouette' in locals() else "N/A",
        'Davies-Bouldin': f"{db_index:.3f}" if 'db_index' in locals() else "N/A",
        'Calinski-Harabasz': f"{ch_score:.3f}" if 'ch_score' in locals() else "N/A",
        'Status': "✅" if 'silhouette' in locals() and silhouette > 0.3 else "⏳"
    },
    'Correlation Analysis': {
        'Z vs Age': f"r = {corr_age:.3f}, p = {pval_age:.4f}" if 'corr_age' in locals() else "N/A",
        'Z vs Health Score': f"r = {corr_hs:.3f}, p = {pval_hs:.4f}" if 'corr_hs' in locals() else "N/A",
        'Target': "r > 0.7",
        'Status': "✅" if 'corr_age' in locals() and abs(corr_age) > 0.7 else "⏳"
    }
}

for category, metrics in metrics_summary.items():
    print(f"\n{category}:")
    for metric, value in metrics.items():
        print("{metric}: {value}")

print("")
print("Metrics calculation complete.")
print("")
print("\nGenerated files:")
print(" ../visualizations/metrics_roc_curve.png - ROC curve for risk stratification")
print(" ../visualizations/metrics_pr_curve.png - Precision-Recall curve")

