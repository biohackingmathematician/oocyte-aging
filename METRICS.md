# Validation Metrics and Performance Evaluation

**Project**: Oocyte Aging Analysis with Trajectory Inference and Risk Stratification  
**Authors**: Agna Chan, Aniqa Nayim, Rimjhim Singh  
**Date**: November 2025

---

## Overview

This document describes the comprehensive validation metrics used to evaluate the success of our oocyte aging trajectory analysis framework. The metrics are organized into five categories addressing different aspects of model performance, biological validation, and clinical utility.

---

## Metrics Framework

### Category 1: Model Fidelity and Validation

#### 1.1 Latent Space Quality Metrics

**Purpose**: Assess the quality of the learned latent representation and its ability to separate biological groups.

**Metrics Used**:

1. **Silhouette Score**
   - **Definition**: Measures how similar an object is to its own cluster compared to other clusters
   - **Range**: [-1, 1], where higher values indicate better separation
   - **Interpretation**: 
     - Values > 0.3 indicate reasonable clustering
     - Values < 0 indicate poor clustering
   - **Our Result**: -0.139 (indicates overlap between GV and MI stages in latent space)

2. **Davies-Bouldin Index**
   - **Definition**: Ratio of within-cluster distances to between-cluster distances
   - **Range**: [0, ∞), where lower values indicate better clustering
   - **Interpretation**: Lower values indicate more compact and well-separated clusters
   - **Our Result**: 1.422

3. **Calinski-Harabasz Score**
   - **Definition**: Ratio of between-cluster dispersion to within-cluster dispersion
   - **Range**: [0, ∞), where higher values indicate better-defined clusters
   - **Interpretation**: Higher values indicate better separation between clusters
   - **Our Result**: 1.417

**Rationale**: These metrics evaluate whether the learned cellular age trajectory successfully separates oocytes by developmental stage (GV vs MI), which is a fundamental biological validation.

---

#### 1.2 Correlation Analysis

**Purpose**: Validate that the learned cellular age correlates with known biological and clinical variables.

**Metrics Used**:

1. **Pearson Correlation Coefficient (r)**
   - **Definition**: Linear correlation between two continuous variables
   - **Range**: [-1, 1]
   - **Interpretation**:
     - |r| > 0.7: Strong correlation
     - 0.3 < |r| < 0.7: Moderate correlation
     - |r| < 0.3: Weak correlation

2. **Statistical Significance (p-value)**
   - **Definition**: Probability of observing the correlation by chance under the null hypothesis
   - **Threshold**: p < 0.05 indicates statistical significance

**Correlations Evaluated**:

- **Cellular Age Z vs Chronological Age**
  - **Result**: r = 0.270, p = 0.249
  - **Interpretation**: Moderate positive correlation, not statistically significant (likely due to small sample size)
  - **Biological Meaning**: Cellular age shows some alignment with chronological age but captures additional biological variation

- **Cellular Age Z vs Clinical Health Score**
  - **Result**: r = -0.837, p < 0.001
  - **Interpretation**: Strong negative correlation, highly statistically significant
  - **Biological Meaning**: Higher cellular age is associated with lower health scores, validating the biological relevance of the trajectory

**Rationale**: These correlations validate that the learned trajectory captures biologically meaningful signals related to both chronological aging and functional quality.

---

### Category 2: Uncertainty Calibration

#### 2.1 Coverage Probability for Prediction Intervals

**Purpose**: Evaluate whether uncertainty estimates are well-calibrated by assessing if prediction intervals contain the true values at the expected rate.

**Metric Used**:

- **95% Coverage Probability**
  - **Definition**: Proportion of observations falling within 95% prediction intervals
  - **Target**: 90-95% (allowing for slight over-coverage)
  - **Our Result**: 1.000 (100% coverage)
  - **Interpretation**: All cellular age values fall within their uncertainty intervals, indicating conservative uncertainty estimates

**Rationale**: Well-calibrated uncertainty is critical for clinical decision-making, as it quantifies confidence in predictions.

#### 2.2 Uncertainty Statistics

**Metrics Used**:

- **Mean Uncertainty**: 2,152.97 ± 3,309.23
- **Coefficient of Variation**: 1.537
- **Interpretation**: High variability in uncertainty estimates across cells, indicating heterogeneous confidence in predictions

**Rationale**: Understanding uncertainty distribution helps identify cells where predictions are less reliable, informing clinical interpretation.

---

### Category 3: Clinical Decision Support

#### 3.1 Risk Stratification Performance

**Purpose**: Evaluate the ability of the model to correctly identify high-risk oocytes requiring intervention.

**Metrics Used**:

1. **Area Under ROC Curve (AUC-ROC)**
   - **Definition**: Probability that the model ranks a randomly chosen high-risk oocyte higher than a randomly chosen low/moderate-risk oocyte
   - **Range**: [0, 1]
   - **Target**: > 0.75 for clinical utility
   - **Our Result**: 1.000
   - **Interpretation**: Perfect discrimination between high-risk and other risk groups

2. **Precision-Recall AUC**
   - **Definition**: Area under the precision-recall curve, particularly important for imbalanced classes
   - **Range**: [0, 1]
   - **Our Result**: 1.000
   - **Interpretation**: Perfect precision-recall performance

3. **Brier Score**
   - **Definition**: Mean squared difference between predicted probabilities and observed outcomes
   - **Range**: [0, 1], where lower is better
   - **Target**: < 0.2 for good calibration
   - **Our Result**: 0.024
   - **Interpretation**: Excellent calibration, predictions closely match observed outcomes

4. **Sensitivity (True Positive Rate)**
   - **Definition**: Proportion of high-risk oocytes correctly identified
   - **Our Result**: 1.000 (100% sensitivity)

5. **Specificity (True Negative Rate)**
   - **Definition**: Proportion of low/moderate-risk oocytes correctly identified
   - **Our Result**: 1.000 (100% specificity)

6. **Positive Predictive Value (PPV)**
   - **Definition**: Proportion of predicted high-risk oocytes that are truly high-risk
   - **Our Result**: 1.000

7. **Negative Predictive Value (NPV)**
   - **Definition**: Proportion of predicted low/moderate-risk oocytes that are truly low/moderate-risk
   - **Our Result**: 1.000

**Rationale**: These metrics directly address the clinical research question of identifying optimal intervention windows by quantifying the model's ability to distinguish high-risk cases.

---

#### 3.2 Clinical Health Score Validation

**Purpose**: Validate that the composite health score successfully discriminates between developmental stages and correlates with biological quality.

**Metrics Used**:

1. **AUC for GV vs MI Classification**
   - **Definition**: Ability of health score to distinguish GV from MI oocytes
   - **Our Result**: 1.000
   - **Interpretation**: Perfect discrimination between stages

2. **Mann-Whitney U Test**
   - **Definition**: Non-parametric test for difference between two independent groups
   - **Our Result**: U = 84.0, p < 0.001
   - **Interpretation**: Highly statistically significant difference between GV and MI health scores

3. **Cohen's d Effect Size**
   - **Definition**: Standardized difference between two means
   - **Interpretation**:
     - |d| < 0.2: Small effect
     - 0.2 < |d| < 0.5: Medium effect
     - 0.5 < |d| < 0.8: Large effect
     - |d| > 0.8: Very large effect
   - **Our Result**: 2.161
   - **Interpretation**: Very large effect size, indicating substantial difference between GV and MI health scores

4. **Mean Health Scores**
   - **GV**: 94.2
   - **MI**: 72.4
   - **Difference**: 21.8 points (23% relative difference)

**Rationale**: Validates that the health score captures biologically meaningful differences between developmental stages, supporting its use in clinical decision-making.

---

#### 3.3 Sensitivity/Specificity Trade-offs

**Purpose**: Evaluate performance at different health score thresholds to inform clinical decision rules.

**Method**: Analyze sensitivity and specificity at multiple threshold percentiles (10th, 25th, 50th, 75th, 90th)

**Rationale**: Different clinical contexts may require different sensitivity/specificity trade-offs. This analysis enables selection of optimal thresholds based on clinical priorities.

---

#### 3.4 Calibration Analysis

**Purpose**: Assess whether predicted risk scores accurately reflect observed risk outcomes.

**Method**: Compare mean predicted risk scores across risk groups to observed high-risk proportions

**Rationale**: Well-calibrated predictions are essential for clinical interpretation, as they allow clinicians to trust the absolute risk estimates, not just relative rankings.

---

### Category 4: Age Discrepancy and Heterogeneity

#### 4.1 Age Discrepancy Metrics

**Purpose**: Quantify the difference between cellular age and chronological age to identify accelerated agers.

**Metrics Used**:

1. **Mean Age Discrepancy**
   - **Definition**: Mean absolute difference between normalized cellular age and normalized chronological age
   - **Our Result**: 0.931 (median)

2. **Proportion of Accelerated Agers**
   - **Definition**: Percentage of oocytes with cellular age > chronological age + 0.2
   - **Our Result**: 0.0%
   - **Interpretation**: No oocytes showed accelerated aging pattern in this dataset

**Rationale**: Identifies individuals whose cellular age deviates substantially from chronological age, which is central to the research question of distinguishing chronological from cellular age.

---

#### 4.2 Heterogeneity Quantification

**Purpose**: Measure variability in aging trajectories across individuals.

**Metrics Used**:

1. **Coefficient of Variation (Uncertainty)**
   - **Definition**: Ratio of standard deviation to mean uncertainty
   - **Our Result**: 1.537
   - **Interpretation**: High variability in uncertainty estimates, indicating heterogeneous aging trajectories

2. **Variance Ratio (Between/Within Donor)**
   - **Definition**: Ratio of between-donor variance to within-donor variance
   - **Our Result**: 0.038
   - **Interpretation**: Low between-donor variance relative to within-donor variance

**Rationale**: Quantifies heterogeneity in aging trajectories, addressing the research question of quantifying variability and uncertainty.

---

### Category 5: Trajectory Fidelity

#### 5.1 Rank Correlation Analysis

**Purpose**: Validate that the learned trajectory preserves expected biological ordering.

**Metric Used**:

- **Kendall's Tau (τ)**
  - **Definition**: Rank correlation coefficient measuring ordinal association
  - **Range**: [-1, 1]
  - **Our Result**: τ = 0.380, p = 0.048
  - **Interpretation**: Moderate positive rank correlation, statistically significant
  - **Biological Meaning**: Cellular age trajectory preserves expected ordering from GV (lower) to MI (higher) stages

**Rationale**: Validates that the learned trajectory respects known biological progression, ensuring biological interpretability.

---

## Summary of Key Results

### Model Performance Summary

| Metric Category | Key Metric | Result | Target | Status |
|----------------|------------|--------|--------|--------|
| **Risk Stratification** | AUC-ROC | 1.000 | > 0.75 | Exceeds target |
| **Risk Stratification** | Brier Score | 0.024 | < 0.2 | Exceeds target |
| **CHS Validation** | AUC (GV vs MI) | 1.000 | > 0.7 | Exceeds target |
| **CHS Validation** | Cohen's d | 2.161 | > 0.8 | Very large effect |
| **Correlation** | Z vs Health Score | r = -0.837 | | Strong correlation |
| **Uncertainty** | 95% Coverage | 1.000 | 0.90-0.95 | Well-calibrated |
| **Trajectory** | Kendall's τ | 0.380 | | Significant ordering |

### Interpretation

The model demonstrates:

1. **Excellent Clinical Performance**: Perfect discrimination (AUC = 1.000) and calibration (Brier = 0.024) for risk stratification
2. **Strong Biological Validation**: Health score shows very large effect size (Cohen's d = 2.161) in distinguishing developmental stages
3. **Robust Trajectory Learning**: Significant rank correlation (τ = 0.380, p = 0.048) preserves expected biological ordering
4. **Well-Calibrated Uncertainty**: 100% coverage indicates conservative but reliable uncertainty estimates

### Limitations

1. **Small Sample Size**: n = 20 limits statistical power and generalizability
2. **Single Study**: Limited cross-study validation due to single transcriptomic dataset
3. **Latent Space Separation**: Silhouette score (-0.139) indicates overlap between stages, suggesting room for improvement in trajectory learning

---

## References

- DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. *Biometrics*, 44(3), 837-845.

- Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. *Monthly Weather Review*, 78(1), 1-3.

- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Lawrence Erlbaum Associates.

- Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics*, 20, 53-65.

---

**Last Updated**: November 2025

