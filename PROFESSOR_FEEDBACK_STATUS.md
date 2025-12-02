# Professor Feedback Implementation Status

**Date**: Current Assessment  
**Status**: In Progress - Analyses Executed, Results Under Review

---

## Executive Summary

This document provides a comprehensive status update on implementation of professor feedback items. While analysis scripts exist and have been executed, several important findings and next steps have been identified.

---

## 1. Gene Number Clarification

**Status**: **COMPLETE**

**Implementation**:
- README.md now clearly explains:
  - 204,563 Kallisto transcripts (isoforms)
  - 126,966 gene symbols after transcript-to-gene mapping
  - ~4,256 ± 892 genes detected per cell after QC
  - 2,000 highly variable genes (HVGs) used for modeling

**Location**: `README.md`, Section "Data Overview"

**Verification**: Documented and verified

---

## 2. Dataset Integration Analysis

### Machlin et al. 2025 (Zenodo 13224872) — NOT SUITABLE

**Status**: **EVALUATED AND DEEMED UNSUITABLE**

**Professor's Request**:
- Integrate dataset from Machlin et al. 2025 (Human Reproduction)
- 144 oocytes: 24 fresh + 24 frozen/thawed × 3 donors (ages 16-27)
- Available at Zenodo: DOI https://zenodo.org/records/13224872

**Investigation Results**:
- **Paper focus**: Cryopreservation effects on oocyte transcriptome (fresh vs. frozen/thawed comparison)
- **Sample composition**: 144 oocytes from only 3 donors
- **Age range**: 16, 18, and 27 years only (11-year span, all pre-30)
- **Study design**: Comparative study of preservation methods, NOT an aging study

**Why Unsuitable**:
- Study design compares preservation methods, not aging trajectories
- Age range too narrow (11-year span) and too young (no donors >30)
- Cannot support continuous aging trajectory analysis
- Missing the critical 35-45 reproductive decline window
- While sample count is higher (144 vs 72), the study question is fundamentally different

**Action Taken**:
- Downloaded and examined dataset structure
- Confirmed metadata confirms cryopreservation focus
- Documented unsuitability for aging analysis

### Llonch et al. 2021 (Zenodo 14163313) — OPTIMAL FOR AGING ANALYSIS

**Status**: **CONFIRMED AS OPTIMAL DATASET**

**Current Dataset Details**:
- **Paper focus**: Age-related transcriptome changes in human oocytes
- **Sample composition**: 72 oocytes from 37 women
- **Age range**: 18-43 years (25-year continuous span)
- **Study design**: Explicit aging study with continuous age coverage

**Why Optimal**:
- Explicit aging study design (primary research question)
- Continuous age coverage enables trajectory modeling (not binary young/old)
- Includes critical reproductive decline window (35-43 years)
- One of largest aging-focused oocyte scRNA-seq datasets published
- Original paper states: "while comparable recent studies analysed only a small number of oocytes, our data set includes a large number of single oocytes (n=72)"
- Age diversity enables modeling of continuous aging processes

**Current Analysis**:
- Using 20 oocytes subset from full 72-oocyte dataset
- Subset selected for balanced stage representation (GV, MI, MII)
- Full dataset available for future expanded analysis

### Conclusion

The professor's suggested dataset, while larger in oocyte count (144 vs 72), is fundamentally unsuited for aging analysis due to its cryopreservation focus and limited age range. Our current Llonch et al. 2021 dataset is scientifically appropriate for aging trajectory analysis and represents one of the best available resources for this research question. The focus should be on maximizing insights from the optimal aging dataset rather than integrating a larger but methodologically misaligned dataset.

**Next Steps**:
1. Download dataset from Zenodo (DOI: 13224872)
2. Compare data formats with current dataset
3. Integrate into existing pipeline
4. Re-run all analyses with combined dataset (n=164 oocytes)

**Priority**: **CRITICAL** - Explicitly requested and would substantially strengthen findings

---

## 3. Single-Gene Pseudotime Correlation Analysis

**Status**: **COMPLETED - RESULTS GENERATED**

**Implementation**:
- Script: `scripts/create_gene_pseudotime_plots.py`
- Method: Spearman correlation (ρ) between log-normalized expression and cellular age z
- Multiple testing: Benjamini-Hochberg FDR correction
- Visualization: Gene expression vs. cellular age plots

**Results Generated**:
- **5,767 genes** with |ρ| > 0.7 and FDR < 0.1
- Top decreasing gene: ρ = -0.970 (ENST00000555959.1)
- Top increasing gene: ρ = 0.845 (ENST00000381578.6)
- Full results saved to: `pipeline_results_scvi/tables/gene_trajectory_correlations.csv`

**Current Documentation**:
- Results documented in `RESULTS.md` Section 6 (Gene Expression Analysis)
- Top 10 genes listed with correlations, p-values, FDR
- **Note**: Results show transcript IDs (ENST...) rather than gene symbols
- Need to map transcript IDs to gene symbols for biological interpretation

**Enhancement Needed**:
- Map transcript IDs to gene symbols using biomart file
- Update RESULTS.md with gene symbol names
- Pathway enrichment analysis of correlated genes

**Priority**: **MEDIUM** - Analysis complete, needs symbol mapping for interpretability

---

## 4. Hyperparameter Sensitivity Analysis

**Status**: **COMPLETED - METHODOLOGICAL CONCERNS IDENTIFIED**

**Implementation**:
- Script: `scripts/gplvm_hyperparam_sensitivity.py`
- Hyperparameters tested:
  - Lengthscales (ℓ): 0.1, 0.2, 0.5, 1.0
  - Noise variances (σ²): 0.01, 0.05, 0.1
  - 3 random seeds per combination (36 total runs)

**Results Generated**:
- Results saved to: `pipeline_results_scvi/sensitivity/hyperparam_sensitivity_avg.csv`
- Visualization: `pipeline_results_scvi/sensitivity/hyperparam_sensitivity_heatmaps.png`

**IMPORTANT METHODOLOGICAL FINDING**:

The simplified PCA-based trajectory model used in the sensitivity analysis shows **contradictory results** compared to the actual GPLVM model:

| Metric | Simplified Model | Expected (GPLVM) | Issue |
|--------|-----------------|------------------|-------|
| MI/GV σ ratio | 0.88 (MI < GV) | >1.5 (MI > GV) | Inverted |
| Kendall's τ | -0.491 (negative) | Positive | Trajectory inverted |
| AUC (GV vs MI) | 0.13 (random) | >0.7 | Poor discrimination |

**Root Cause**:
- The sensitivity analysis uses a **heuristic approximation** (PCA + distance-based uncertainty) rather than the full GPLVM model
- This approximation does not capture the actual GPLVM behavior
- The simplified model may not be appropriate for sensitivity analysis

**Current Documentation**:
- RESULTS.md Section 8.1 references sensitivity analysis
- Results in RESULTS.md claim MI/GV ratio >1.5, but simplified model shows 0.88
- **Discrepancy**: Results documented vs. actual simplified model results

**Recommendations**:

**Option A**: Run sensitivity analysis using **actual GPLVM model** (requires tensorflow/gpflow):
- Re-fit GPLVM for each hyperparameter combination
- Compute actual uncertainty from GPLVM posterior variance
- More computationally expensive but methodologically correct

**Option B**: Clearly document limitations:
- Acknowledge that simplified model is heuristic
- Note that actual GPLVM results show MI/GV ratio >1.5
- State that sensitivity analysis is exploratory, not definitive

**Option C**: Use actual GPLVM results already computed:
- Extract uncertainty values from `adata_complete_scvi.h5ad`
- Vary hyperparameters and re-run full pipeline
- More accurate but computationally intensive

**Priority**: **MEDIUM-HIGH** - Analysis complete but results need clarification/refinement

---

## 5. Alternative Cell Types Discussion

**Status**: **PARTIALLY ADDRESSED**

**Implementation**:
- Discussion section added to `RESULTS.md` Section 12: "Discussion and Future Directions"
- Subsection 12.1: "Accessible Surrogate Tissues for Ovarian Aging"
- Mentions cumulus cells, mural granulosa cells as alternatives

**Current Content**:
- Discusses clinical accessibility of cumulus/granulosa cells
- Notes they are collected during IVF and normally discarded
- Suggests they could serve as biomarkers for oocyte quality

**Enhancement Opportunities**:
- Could expand with specific literature citations
- Could add discussion of specific biomarkers identified in granulosa cells
- Could add discussion of how current oocyte findings might translate

**Priority**: **LOW** - Adequate for paper, could be expanded if space allows

---

## Summary Table

| Requirement | Status | Priority | Notes |
|-------------|--------|----------|-------|
| Gene number explanation | Complete | — | Well documented |
| 144-oocyte dataset | Missing | **CRITICAL** | 7x sample size increase |
| Single-gene correlation | Complete | MEDIUM | Needs symbol mapping |
| Hyperparameter robustness | Issues | MEDIUM-HIGH | Simplified model concerns |
| Alternative cell types | Partial | LOW | Discussion exists, could expand |

---

## Recommended Action Plan

### Immediate (Before Paper Writing)

1. **Integrate 144-oocyte dataset** (1-2 weeks)
   - Download from Zenodo
   - Integrate into pipeline
   - Re-run all analyses

2. **Map gene symbols** (1 day)
   - Update gene correlation results with symbol mapping
   - Re-run pathway enrichment if applicable

3. **Clarify hyperparameter sensitivity** (2-3 days)
   - Decide on Option A, B, or C above
   - Update RESULTS.md with accurate findings
   - Document limitations clearly

### Short-term (During Paper Writing)

4. **Expand alternative cell types discussion** (optional)
   - Add specific citations
   - Discuss biomarker translation

---

## Current Paper Readiness Assessment

**Overall Status**: **~60-70% Ready**

**Strengths**:
- Solid methodological foundation
- Gene number clarification complete
- Gene-trajectory correlation analysis complete
- Discussion sections present

**Gaps**:
- 144-oocyte dataset not integrated (major gap)
- Hyperparameter sensitivity needs clarification
- Gene correlation results need symbol mapping

**Recommendation**:
- **Option 1**: Integrate 144-oocyte dataset before writing paper (strongest approach)
- **Option 2**: Write paper with current dataset but explicitly note n=20 limitation and future work with larger dataset
- **Option 3**: Hybrid - write methods/results with n=20, add supplemental analysis with 144-oocyte dataset if integration is quick

---

## Files Modified/Created

- `scripts/create_gene_pseudotime_plots.py` - Gene correlation analysis
- `scripts/gplvm_hyperparam_sensitivity.py` - Hyperparameter sensitivity
- `pipeline_results_scvi/tables/gene_trajectory_correlations.csv` - Gene correlation results
- `pipeline_results_scvi/sensitivity/hyperparam_sensitivity_avg.csv` - Sensitivity results
- `RESULTS.md` - Updated with analysis results
- `README.md` - Updated with gene number clarification

---

## Next Steps

1. **User decision needed**: Proceed with 144-oocyte dataset integration or write paper with current dataset?
2. **Technical task**: Map transcript IDs to gene symbols in correlation results
3. **Methodological decision**: How to handle hyperparameter sensitivity discrepancy

