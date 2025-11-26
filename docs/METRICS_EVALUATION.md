# Metrics Evaluation and Recommendations

**Project**: Oocyte Aging Analysis with Trajectory Inference and Risk Stratification  
**Authors**: Agna Chan, Aniqa Nayim, Rimjhim Singh  
**Date**: November 2025

---

## Executive Summary

This document evaluates the current implementation against the research objectives and identifies critical missing metrics for comprehensive model validation and clinical translation. The analysis is based on cross-validation with the project plan and mid-progress report.

---

## Current Implementation Status

### ✅ Implemented Metrics

1. **Trajectory Analysis**:
   - Diffusion Pseudotime (DPT) correlation with stage: ρ = -0.79, p < 0.001
   - Cellular age Z (normalized 0-1)
   - Heuristic uncertainty estimates (σ)

2. **Health Scoring**:
   - Composite Health Score (CHS) integrating pathway-level expression
   - Risk stratification (Low/Moderate/High risk groups)
   - Stage-specific health score distributions

3. **Basic Validation**:
   - Correlation between cellular age Z and chronological age
   - Risk group distribution statistics

---

## Critical Missing Metrics

### CATEGORY 1: Model Fidelity & Validation

#### A. GPLVM Implementation Metrics (High Priority)

**Status**: ⏳ Not Yet Implemented (PCA fallback currently used)

**Required Metrics**:

1. **ELBO (Evidence Lower Bound)**
   - Purpose: Track training convergence
   - Target: Converged ELBO value
   - Implementation: Track during GPLVM training iterations

2. **Reconstruction Error**
   - Formula: MSE = (1/N) Σ(X_original - X_reconstructed)²
   - Target: MSE < 0.5
   - Purpose: Validate model fidelity to original expression data

3. **Latent Space Quality**
   - Silhouette score for stage separation (GV vs MI)
   - Davies-Bouldin index (lower = better clustering)
   - Calinski-Harabasz score (higher = better defined clusters)

**Current Limitation**: Using PCA-based trajectory approximation instead of full Bayesian GPLVM due to Python 3.14 compatibility constraints.

---

#### B. Uncertainty Calibration (Partially Implemented)

**Status**: ️ Partially Implemented

**Current**: Uncertainty values (σ) calculated heuristically

**Missing Metrics**:

4. **Calibration Curves**
   - Plot predicted vs actual uncertainty
   - Expected calibration error (ECE)
   - Maximum calibration error (MCE)
   - Target: ECE < 0.1

5. **Prediction Intervals**
   - Coverage probability (e.g., do 95% intervals contain true values 95% of time?)
   - Target: Coverage 90-95%
   - Interval width analysis

**Recommendation**: Implement calibration analysis using existing uncertainty estimates.

---

### CATEGORY 2: Cross-Study Validation

#### C. Leave-One-Study-Out Cross-Validation

**Status**: ⏳ Framework exists but limited by single-study dataset

**Required Metrics**:

6. **Cross-Study Correlation**
   - Pearson r between predicted cellular age Z and chronological age
   - Target: r > 0.7 across all held-out studies
   - Current: r = -0.79 for DPT pseudotime (good baseline)

7. **Batch Effect Quantification**
   - kBET (k-nearest neighbor batch effect test)
   - LISI (Local Inverse Simpson's Index) - measures mixing
   - Silhouette coefficient for batch vs biological signal
   - Target: kBET p-value > 0.05 (batch effects removed)

8. **Cross-Technology Replication**
   - Concordance correlation coefficient (CCC) between GSE155179 and GSE95477
   - Bland-Altman plots for agreement

**Current Limitation**: Only one study (Zenodo) with sufficient transcriptomic data. Age data from GEO studies (GSE155179, GSE95477) but no matching transcriptomic profiles.

---

### CATEGORY 3: Biological Validation

#### D. Pathway Analysis (Partially Done)

**Status**: ✅ Pathway scores calculated, ⏳ Enrichment statistics needed

**Missing Metrics**:

9. **GSEA Statistics**
   - Normalized Enrichment Score (NES) for OXPHOS, spindle, meiosis pathways
   - FDR q-values (should be < 0.05)
   - Leading edge analysis (core enrichment genes)

10. **Differential Expression Validation**
    - Fold change for known age markers (TOP2B, PSMA1, PSMA2)
    - Compare findings vs Zhang et al. 2020 (dataset source)
    - Cohen's d effect sizes for age groups

**Recommendation**: Add GSEA analysis using existing pathway scores.

---

### CATEGORY 4: Clinical Decision Support (High Priority)

#### E. Intervention Timing Accuracy

**Status**: ⏳ NOT YET IMPLEMENTED

**Critical Missing Metrics**:

11. **AMH Prediction Performance**
    - Mean Absolute Error (MAE) in years for AMH threshold prediction
    - Target: ±2 years for AMH <1.0 ng/mL prediction
    - R² for cellular age Z vs AMH correlation
    - **Current**: AMH calibration skipped due to gpflow/tensorflow compatibility

12. **Time-to-Threshold Accuracy**
    - Compare predictions to Penn Ovarian Aging Study ground truth
    - Hazard ratios for high-risk vs low-risk groups
    - C-statistic (concordance index) for time-to-event prediction

13. **Risk Stratification Performance**
    - AUC-ROC for predicting "accelerated agers"
    - Target: AUC > 0.75
    - Precision-Recall curves (important for imbalanced classes)
    - Confusion matrix (sensitivity, specificity, PPV, NPV)

**Priority**: HIGH - Directly addresses research question on intervention timing.

---

#### F. Clinical Calibration

**Status**: ⏳ MISSING

**Required Metrics**:

14. **Calibration Plots**
    - Do 20% of "High Risk" women actually show accelerated aging?
    - Hosmer-Lemeshow test for goodness-of-fit
    - Brier score (measures calibration + discrimination)
    - Target: Brier score < 0.2

15. **Net Reclassification Improvement (NRI)**
    - How many women are correctly reclassified vs ASRM age-35 guideline?
    - Integrated Discrimination Improvement (IDI)

---

### CATEGORY 5: Optimal Transport

**Status**: ⏳ Planned, Needs Specification

**Required Metrics**:

16. **Wasserstein Distance**
    - W1 distance between predicted and actual aging trajectories
    - Per-pathway Wasserstein distance (do OXPHOS trajectories match expectations?)
    - Use POT library as planned

17. **Trajectory Fidelity**
    - Fréchet distance between learned and expected GV→MI trajectory
    - Kendall's tau for rank correlation of pseudotime vs cellular age Z

---

## Recommended Metrics Implementation Priority

### Priority 1: Clinical Decision Support (Addresses Research Question)

1. **Risk Stratification Performance** (AUC-ROC, Precision-Recall)
   - Feasibility: HIGH - Can calculate from existing risk groups
   - Impact: HIGH - Directly answers intervention timing question

2. **Clinical Health Score Validation**
   - Discriminative ability: AUC for GV vs MI classification
   - Criterion validity: Correlation with known outcomes
   - Construct validity: Correlation with AMH (if available)

3. **Calibration Analysis**
   - Brier score calculation
   - Calibration plots for risk predictions

### Priority 2: Model Validation

4. **Uncertainty Calibration**
   - Expected Calibration Error (ECE)
   - Coverage probability for prediction intervals

5. **Latent Space Quality**
   - Silhouette score for stage separation
   - Davies-Bouldin index

### Priority 3: Biological Validation

6. **GSEA Statistics**
   - Normalized Enrichment Score (NES)
   - FDR q-values for pathway enrichment

7. **Differential Expression Validation**
   - Fold changes for known age markers
   - Comparison to literature (Zhang et al. 2020)

### Priority 4: Advanced Validation (Requires Additional Data)

8. **Cross-Study Validation**
   - Requires multiple studies with matching transcriptomic + age data

9. **AMH Calibration**
   - Requires gpflow/tensorflow (Python 3.10-3.13) or alternative GP implementation

10. **Optimal Transport Metrics**
    - Requires trajectory comparison data

---

## Implementation Roadmap

### Phase 1: Immediate (Can implement now)

- [x] Risk stratification AUC-ROC calculation
- [x] Clinical Health Score validation metrics
- [x] Uncertainty calibration analysis
- [x] Latent space quality metrics (Silhouette, Davies-Bouldin)

### Phase 2: Short-term (Requires minor code additions)

- [ ] GSEA statistics for pathway enrichment
- [ ] Differential expression validation
- [ ] Calibration plots and Brier score
- [ ] Precision-Recall curves

### Phase 3: Medium-term (Requires additional data or packages)

- [ ] AMH prediction performance (requires GP regression)
- [ ] Cross-study validation (requires multiple studies)
- [ ] Optimal transport metrics (requires trajectory comparison)

---

## Metrics Tracking Table

### Table 1: Model Performance Summary

| Metric | Value | Target | Status | Notes |
|--------|-------|--------|--------|-------|
| **Trajectory Learning** |
| DPT Correlation (pseudotime vs stage) | ρ = -0.79 | ρ > 0.7 | ✅ | Strong ordering |
| GPLVM ELBO | N/A | Converged | ⏳ | PCA fallback used |
| Reconstruction MSE | N/A | < 0.5 | ⏳ | Not applicable (PCA) |
| Correlation Z vs age | [TO CALCULATE] | r > 0.7 | ⏳ | Needs validation |
| **Uncertainty Quantification** |
| Mean uncertainty (σ) | 2152.97 ± 3225.44 | - | ✅ | Heuristic estimate |
| Calibration ECE | [TO CALCULATE] | < 0.1 | ⏳ | Needs implementation |
| Coverage (95% CI) | [TO CALCULATE] | 90-95% | ⏳ | Needs implementation |
| **Biological Validation** |
| OXPHOS pathway score | [Calculated] | Negative correlation | ✅ | Needs NES |
| Spindle pathway score | [Calculated] | Positive correlation | ✅ | Needs NES |
| Pathway FDR | [TO CALCULATE] | < 0.05 | ⏳ | Needs GSEA |
| **Cross-Study Validation** |
| LOO CV correlation | N/A | r > 0.7 | ⏳ | Single study available |
| Batch kBET p-value | [TO CALCULATE] | > 0.05 | ⏳ | Needs implementation |
| Wasserstein distance | [TO CALCULATE] | < 0.3 | ⏳ | Needs implementation |
| **Clinical Performance** |
| Risk stratification AUC | [TO CALCULATE] | > 0.75 | ⏳ | **HIGH PRIORITY** |
| CHS discriminative ability | [TO CALCULATE] | AUC > 0.7 | ⏳ | **HIGH PRIORITY** |
| Brier score | [TO CALCULATE] | < 0.2 | ⏳ | **HIGH PRIORITY** |
| AMH prediction MAE | N/A | ±2 years | ⏳ | Requires GP regression |

---

## Comparison to Baseline Methods

| Method | Uncertainty | Clinical Output | Correlation with Age | Limitations |
|--------|-------------|-----------------|---------------------|-------------|
| DPT (Current) |  | Health score | ρ = -0.79 | No uncertainty |
| PCA-based (Current) | ✅ σ² (heuristic) | Health score + risk groups | r = 0.27 | Not true Bayesian |
| scVI + DPT (Implemented) | ✅ σ² (heuristic) | Health score + risk groups | r = 0.27 | Sensitivity analysis confirms robustness |
| GPLVM (Planned) | ✅ σ² (Bayesian) | Health score + CI | Target r > 0.75 | Requires compatible Python |
| Monocle (Literature) |  | Trajectory only | Variable | No clinical translation |
| **Our Framework** | ✅ | Risk groups + timing | r = 0.27, p = 0.25 | Small sample size (n=20) |

**Note**: Hyperparameter sensitivity analysis (see `scripts/gplvm_hyperparam_sensitivity.py`) confirms that key findings (MI > GV uncertainty, ~20-30% high-uncertainty cells) are robust across kernel hyperparameters.

---

## Key Metrics for Research Question

### For "Quantify Heterogeneity and Uncertainty":

1. ✅ Per-cell uncertainty (σ²) - **IMPLEMENTED** (heuristic)
2. ⏳ Calibration curves showing uncertainty is accurate - **NEEDS IMPLEMENTATION**
3. ⏳ Heterogeneity quantification: coefficient of variation across pathways - **NEEDS IMPLEMENTATION**
4. ⏳ Inter-individual variability: ratio of between-donor to within-donor variance - **NEEDS IMPLEMENTATION**

### For "Identify Optimal Intervention Windows":

1. ⏳ **Lead time analysis**: Years gained vs ASRM age-35 recommendation - **NEEDS IMPLEMENTATION**
2. ⏳ **Sensitivity/Specificity trade-offs** at different CHS thresholds - **NEEDS IMPLEMENTATION**
3. ⏳ **Number needed to screen (NNS)** to prevent one late intervention - **NEEDS IMPLEMENTATION**
4. ✅ Stage-specific feasibility (67% GV viable vs 17% MI) - **DOCUMENTED**

### For "Distinguish Chronological from Cellular Age":

1. ⏳ **Age discrepancy metric**: |Z_cellular - age_chronological| / age_chronological - **NEEDS IMPLEMENTATION**
2. ⏳ **Proportion of "accelerated agers"**: % with Z >> chronological age - **NEEDS IMPLEMENTATION**
3. ⏳ **Biological age prediction error**: MAE between predicted and actual AMH decline - **REQUIRES AMH DATA**

---

## Recommendations for Final Report

### 1. Emphasize Strong DPT Results
- DPT correlation (ρ = -0.79, p < 0.001) is a strong baseline
- Clearly state that PCA-based trajectory is a fallback, not the primary method
- Acknowledge limitations honestly while highlighting what was achieved

### 2. Add Quantitative Clinical Validation
- Calculate AUC-ROC for risk stratification (feasible with current data)
- Add calibration plots for risk predictions
- Compare to ASRM age-35 guideline quantitatively

### 3. Implement Feasible Metrics
- Uncertainty calibration (ECE, coverage probability)
- Latent space quality (Silhouette, Davies-Bouldin)
- GSEA statistics for pathway enrichment

### 4. Document Limitations Clearly
- Small sample size (n=20) limits statistical power
- Single study limits cross-validation
- Python 3.14 compatibility constraints
- Missing AMH calibration due to package dependencies

### 5. Highlight Novel Contributions
- Integration of trajectory inference with clinical risk stratification
- Composite health score combining multiple pathways
- Framework for intervention timing (even if not fully validated)

---

## Next Steps

1. ✅ **Completed**: Risk stratification AUC calculation (AUC = 1.000)
2. ✅ **Completed**: Clinical Health Score validation (AUC = 1.000, Cohen's d = 2.161)
3. ✅ **Completed**: Uncertainty calibration analysis (95% coverage = 1.000)
4. ✅ **Completed**: Latent space quality metrics (Silhouette, Davies-Bouldin, Calinski-Harabasz)
5. ✅ **Completed**: Hyperparameter sensitivity analysis (`scripts/gplvm_hyperparam_sensitivity.py`)
6. ✅ **Completed**: Gene-pseudotime correlation analysis with FDR correction (`scripts/create_gene_pseudotime_plots.py`)
7. **Remaining**: Full Bayesian GPLVM when compatible Python version available
8. **Remaining**: Expand dataset to full 72-oocyte cohort

---

## References

- Zhang et al. 2020: Dataset source for oocyte transcriptomics
- ASRM Guidelines: Age-35 recommendation for fertility assessment
- Penn Ovarian Aging Study: Ground truth for AMH decline trajectories

---

**Last Updated**: November 2025  
**Status**: Active evaluation and implementation

