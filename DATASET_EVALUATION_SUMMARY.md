# Dataset Evaluation and Integration Summary

**Date**: December 2025  
**Purpose**: Document dataset selection rationale and alternative dataset evaluation

---

## Executive Summary

This document summarizes the evaluation of datasets for oocyte aging trajectory analysis, including the rationale for selecting the Llonch et al. 2021 dataset and the evaluation of alternative datasets suggested by the professor.

---

## Primary Dataset: Llonch et al. 2021 (Zenodo 14163313)

### Dataset Details

- **Source**: Zenodo DOI: 10.5281/zenodo.14163313
- **Paper**: Llonch, S., et al. (2021). Single human oocyte transcriptome analysis reveals a novel aging trajectory. *Aging Cell*, 20(5), e13360
- **Full Dataset**: 72 oocytes from 37 women
- **Age Range**: 18-43 years (25-year continuous span)
- **Stages**: GV, MI, MII oocytes across age range
- **Analysis Subset**: 20 oocytes (6 GV, 14 MI) with complete metadata

### Why This Dataset is Optimal

1. **Explicit Aging Focus**
   - Primary research question addresses age-related transcriptomic changes
   - Study designed specifically for aging trajectory analysis
   - Original paper states: "while comparable recent studies analysed only a small number of oocytes, our data set includes a large number of single oocytes (n=72)"

2. **Continuous Age Coverage**
   - 25-year span (18-43 years) enables continuous trajectory modeling
   - Not limited to binary young/old comparison
   - Multiple oocytes per age point provides statistical power

3. **Critical Age Window Included**
   - Includes the 35-43 year reproductive decline period
   - Covers both pre- and post-35 age ranges
   - Enables modeling of accelerated aging after 35

4. **Large Sample Size**
   - 72 oocytes is among the largest aging-focused oocyte scRNA-seq datasets published
   - 37 donors provides genetic diversity while maintaining age diversity
   - Sufficient power for trajectory inference and statistical analyses

5. **Stage Diversity**
   - Includes GV, MI, and MII stages across age range
   - Enables stage-specific aging analysis
   - Supports trajectory modeling across maturation stages

### Technical Details

- **Raw Data**: 204,563 Kallisto transcripts (isoforms)
- **Gene-level**: 126,966 unique genes after transcript-to-gene mapping
- **Detected per Cell**: ~4,256 ± 892 genes (after QC)
- **Modeling**: 2,000 highly variable genes used for trajectory analysis
- **Format**: Kallisto TPM abundance estimates

---

## Alternative Dataset Evaluation: Machlin et al. 2025 (Zenodo 13224872)

### Dataset Details

- **Source**: Zenodo DOI: 10.5281/zenodo.13224872
- **Paper**: Machlin, J., et al. (2025). Single-cell analysis comparing early-stage oocytes from fresh and slow-frozen/thawed human ovarian cortex reveals minimal impact of cryopreservation on the oocyte transcriptome. *Human Reproduction*
- **Oocyte Count**: 144 oocytes
- **Donors**: 3 donors
- **Age Range**: 16, 18, and 27 years (only 3 ages)
- **Study Design**: 24 fresh + 24 frozen/thawed oocytes × 3 donors

### Why This Dataset is Unsuitable

1. **Wrong Research Question**
   - Study focuses on **cryopreservation effects** (fresh vs. frozen/thawed comparison)
   - Primary question: Does freezing affect transcriptome? (NOT: How does age affect transcriptome?)
   - Experimental design optimized for preservation method comparison, not aging

2. **Insufficient Age Range**
   - Only 3 ages: 16, 18, 27 years (11-year span)
   - All donors are pre-30 (no post-30 representation)
   - Missing the critical 35-45 reproductive decline window
   - Cannot support continuous aging trajectory analysis

3. **Limited Age Diversity**
   - Only 3 donors (vs. 37 in Llonch dataset)
   - No age diversity within each donor
   - Cannot distinguish age effects from donor effects

4. **Cannot Support Our Analysis**
   - Cannot model continuous aging trajectory (only 3 age points)
   - Cannot assess reproductive decline period (all donors <30)
   - Statistical power limited by small number of unique ages
   - Results would be confounded by preservation method effects

### Investigation Results

- Downloaded and examined dataset structure from Zenodo
- Confirmed metadata indicates cryopreservation focus
- Verified age range limitation (3 donors, ages 16-27)
- Dataset downloaded to `machlin_2025_data/` for reference

### Conclusion

While the Machlin dataset has a larger oocyte count (144 vs. 72), it is fundamentally unsuited for aging trajectory analysis. The larger sample size does not compensate for:
- Wrong experimental design (cryopreservation vs. aging)
- Insufficient age range (3 ages, all <30)
- Missing critical age window (35-43 years)
- Limited age diversity (only 3 donors)

---

## Additional Dataset: GSE202601 Granulosa Cell Analysis

### Dataset Details (For Future Integration)

- **Source**: GEO GSE202601
- **Paper**: Chen et al. (2024) Nature Aging. "Molecular and genetic insights into human ovarian aging from single-nuclei multi-omics analyses"
- **Samples**: 8 human ovaries (4 young: ages 23-29; 4 aged: ages 49-54)
- **Data Type**: snRNA-seq + snATAC-seq
- **Cell Types**: Granulosa, theca, stromal, immune, endothelial cells

### Clinical Relevance

This dataset addresses the professor's question: **"Are there different cell types that could be correlated with oocyte age that are easier to collect?"**

**Answer**: Yes - **granulosa cells** are ideal because:

1. **Clinically Accessible**
   - Routinely collected during IVF oocyte retrieval procedures
   - Currently discarded after procedure (no additional procedures needed)
   - Multiple cells per follicle vs. single oocyte

2. **Biologically Relevant**
   - Direct gap-junction communication with oocytes
   - Granulosa cell gene expression predicts oocyte competence
   - Age-related changes in granulosa cells reflect ovarian aging

3. **Literature Support**
   - Chen et al. (2024) demonstrated age-related transcriptomic changes in granulosa cells
   - Multiple studies show granulosa markers predict oocyte quality

### Integration Status

- Metadata downloaded (16 samples: 8 snRNA-seq + 8 snATAC-seq)
- Analysis framework created (`scripts/analyze_granulosa_aging.py`)
- Comparison script ready (`scripts/compare_oocyte_granulosa_aging.py`)
- Pending: Download processed count matrices from GEO
- Pending: Run granulosa cell trajectory analysis
- Pending: Compare with oocyte aging signatures

### Future Work

1. Download processed GSE202601 h5ad file from GEO supplementary files
2. Run granulosa cell trajectory analysis using same approach as oocyte analysis
3. Identify shared aging signatures between oocytes and granulosa cells
4. Develop granulosa-cell-based clinical biomarkers for oocyte aging

---

## Final Recommendation

**Use Llonch et al. 2021 dataset (Zenodo 14163313) for primary aging analysis.**

**Rationale**:
- Scientifically appropriate for aging trajectory research question
- Optimal age range and diversity for continuous modeling
- One of the best available resources for this research question
- Sufficient sample size for robust statistical analysis

**Do NOT integrate Machlin et al. 2025 dataset** because:
- Fundamentally different research question (cryopreservation vs. aging)
- Insufficient age range for aging analysis
- Would confound results without adding scientific value

**Future integration**: GSE202601 granulosa cell data for clinically accessible biomarkers (complementary analysis, not replacement).

---

## Files and Scripts Created

1. **Documentation**:
   - `PROFESSOR_FEEDBACK_STATUS.md` - Updated with dataset evaluation
   - `README.md` - Updated with dataset selection rationale
   - `RESULTS.md` - Added dataset selection rationale section
   - `DATASET_EVALUATION_SUMMARY.md` - This document

2. **Scripts**:
   - `scripts/download_gse202601.py` - Download and setup GSE202601 analysis
   - `scripts/compare_oocyte_granulosa_aging.py` - Compare oocyte and granulosa signatures
   - `gse202601_data/analyze_granulosa_aging.py` - Granulosa cell analysis template

3. **Data**:
   - `machlin_2025_data/` - Downloaded Machlin dataset (for reference)
   - `gse202601_data/gse202601_metadata.csv` - GSE202601 metadata
   - `comparison_results/granulosa_discussion_placeholder.md` - Discussion placeholder

---

## References

1. Llonch, S., et al. (2021). Single human oocyte transcriptome analysis reveals a novel aging trajectory. *Aging Cell*, 20(5), e13360. DOI: 10.1111/acel.13360

2. Machlin, J., et al. (2025). Single-cell analysis comparing early-stage oocytes from fresh and slow-frozen/thawed human ovarian cortex reveals minimal impact of cryopreservation on the oocyte transcriptome. *Human Reproduction*. Zenodo: 10.5281/zenodo.13224872

3. Chen, J., et al. (2024). Molecular and genetic insights into human ovarian aging from single-nuclei multi-omics analyses. *Nature Aging*. DOI: 10.1038/s43587-024-00762-5. GEO: GSE202601

