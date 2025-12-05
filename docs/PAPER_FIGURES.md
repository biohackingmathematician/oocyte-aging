# Paper Figures Inventory

## Main Figures

### Figure 1: GPLVM Cellular Age Trajectory

- **File**: `visualizations/gplvm_trajectory_analysis.png`
- **Content**: UMAP visualization colored by cellular age with uncertainty
- **Status**: ✅ Ready

### Figure 2: Risk Stratification

- **File**: `visualizations/risk_stratification.png`
- **Content**: 
  - A) Risk group distribution (pie/bar chart)
  - B) Feature heatmap by risk group
- **Status**: ✅ Ready

### Figure 3: Gene-Trajectory Correlations (TO CREATE)

- **File**: `visualizations/gene_trajectory_volcano.png`
- **Content**: Volcano plot of gene correlations vs. pseudotime
- **Status**: 🔨 Create with script below

### Figure 4: Complete Results Summary

- **File**: `visualizations/complete_results_summary.png`
- **Content**: Multi-panel summary figure
- **Status**: ✅ Ready

## Supplementary Figures

### Fig S1: Metrics ROC Curve

- **File**: `visualizations/metrics_roc_curve.png`
- **Status**: ✅ Ready

### Fig S2: Precision-Recall Curve

- **File**: `visualizations/metrics_pr_curve.png`
- **Status**: ✅ Ready

### Fig S3: Hyperparameter Sensitivity (TO CREATE)

- **File**: `visualizations/hyperparameter_sensitivity.png`
- **Content**: Heatmap/line plot of metrics across configurations
- **Status**: 🔨 Create if sensitivity results exist

---

## Figure Generation Scripts

### Create Gene Correlation Volcano Plot

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load correlations
df = pd.read_csv('pipeline_results_scvi/tables/gene_trajectory_correlations_with_symbols.csv')

# Create volcano plot
fig, ax = plt.subplots(figsize=(10, 8))

# Color by significance
colors = np.where(
    (df['fdr'] < 0.1) & (abs(df['correlation']) > 0.5),
    np.where(df['correlation'] > 0, 'red', 'blue'),
    'gray'
)

ax.scatter(df['correlation'], -np.log10(df['fdr']), 
           c=colors, alpha=0.5, s=10)

# Add labels for top genes
top_genes = df.nsmallest(10, 'fdr')
for _, row in top_genes.iterrows():
    ax.annotate(row['gene_symbol'], 
                (row['correlation'], -np.log10(row['fdr'])),
                fontsize=8)

ax.axhline(-np.log10(0.1), color='gray', linestyle='--', alpha=0.5)
ax.axvline(0.5, color='gray', linestyle='--', alpha=0.5)
ax.axvline(-0.5, color='gray', linestyle='--', alpha=0.5)

ax.set_xlabel('Spearman ρ (correlation with cellular age)')
ax.set_ylabel('-log10(FDR)')
ax.set_title('Gene-Trajectory Correlations')

plt.tight_layout()
plt.savefig('visualizations/gene_trajectory_volcano.png', dpi=300)
plt.close()
```

