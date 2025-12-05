# Paper Draft Outline

## Working Title

**Evaluating Oocyte Aging Uncertainty: A Multidimensional Bayesian Approach for Personalized Fertility Preservation Timing**

## Target

- Course: STAT GR5243 Applied Data Science
- Format: Research paper (~4000-6000 words)
- Deadline: [INSERT]

---

## Abstract (250 words)

- [ ] Problem statement (1-2 sentences)
- [ ] Methods overview (2-3 sentences)
- [ ] Key results (3-4 sentences)
- [ ] Conclusion (1-2 sentences)

## 1. Introduction (~800 words)

- [ ] 1.1 Oocyte aging and fertility decline
- [ ] 1.2 Limitations of chronological age as predictor
- [ ] 1.3 Single-cell transcriptomics opportunity
- [ ] 1.4 Gap: uncertainty quantification lacking
- [ ] 1.5 Our contribution (3 bullet points)

## 2. Methods (~1200 words)

- [ ] 2.1 Data Sources
  - Primary dataset: Llonch et al. 2021
  - Age integration: GSE155179, GSE95477
  - Dataset selection rationale (why not Machlin)
- [ ] 2.2 Preprocessing
  - Kallisto quantification
  - Gene symbol mapping
  - QC and filtering
  - HVG selection
- [ ] 2.3 Trajectory Analysis
  - GPLVM formulation
  - Uncertainty quantification
  - Hyperparameter choices
- [ ] 2.4 Clinical Risk Stratification
  - Health score calculation
  - Risk group definitions
- [ ] 2.5 Gene Correlation Analysis
  - Spearman correlation
  - FDR correction
- [ ] 2.6 Validation Approach
  - Metrics used
  - Cross-validation strategy

## 3. Results (~1200 words)

- [ ] 3.1 Cellular Age Trajectory (Figure 1)
  - Trajectory visualization
  - Correlation with chronological age
  - Uncertainty patterns
- [ ] 3.2 Risk Stratification (Figure 2)
  - Three risk groups
  - Group characteristics
  - Clinical health score validation
- [ ] 3.3 Aging-Associated Genes (Figure 3, Table 1)
  - 3,683 significant genes
  - Top genes and biological interpretation
  - Pathway enrichment
- [ ] 3.4 Model Validation (Table 2)
  - Discriminative metrics
  - Sensitivity analysis results

## 4. Discussion (~1000 words)

- [ ] 4.1 Summary of findings
- [ ] 4.2 Comparison with prior work
- [ ] 4.3 Clinical implications
- [ ] 4.4 Limitations
  - Sample size
  - Single dataset
  - Hyperparameter sensitivity scope
- [ ] 4.5 Future Directions
  - Granulosa cell biomarkers
  - Multi-study integration
  - Prospective validation

## 5. Conclusion (~200 words)

- [ ] Key takeaways
- [ ] Clinical relevance
- [ ] Call to action

## References

- [ ] Compile from docs/references/

## Supplementary Materials

- [ ] Full gene correlation table
- [ ] Additional figures
- [ ] Code availability statement

---

## Key Statistics to Include

| Metric | Value |
|--------|-------|
| Total oocytes | 72 |
| Analysis subset | 20 |
| Age range | 18-43 years |
| Genes tested | 126,966 |
| Significant genes | 3,683 |
| Risk groups | 3 |
| AUC-ROC | 1.000 |
| Cohen's d | 2.161 |
| Health score correlation | r = -0.79 |

---

## Figures Checklist

- [x] Figure 1: Trajectory visualization
- [x] Figure 2: Risk stratification
- [ ] Figure 3: Gene volcano plot
- [x] Figure 4: Summary panel
- [x] Fig S1: ROC curve
- [x] Fig S2: PR curve

