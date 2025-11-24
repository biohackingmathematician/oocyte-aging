# Data Integration Plan: Sample Size Expansion

**Priority**: 1 (Critical for Forum Readiness)  
**Status**: Planning Phase  
**Target**: Expand from n=20 to n=50+ cells

---

## Current Status

### Existing Data

**Zenodo Dataset (Primary)**:
- **Samples**: 20 oocytes
- **Stages**: 6 GV, 14 MI, partial MII
- **Age Range**: 25-35 years (mapped from GEO metadata)
- **Data Type**: Single-cell RNA-seq (kallisto abundance.tsv files)
- **Status**: Fully integrated and analyzed

**GEO Datasets (Metadata Only)**:
- **GSE155179**: 12 samples with age metadata (6 <30 years, 6 ≥40 years)
- **GSE95477**: 32 samples with age metadata (estimated 25-35 years)
- **Status**: Age metadata extracted and mapped to Zenodo samples
- **Limitation**: Transcriptomic expression data not yet integrated

---

## Integration Requirements

### GSE155179 Integration

**Dataset Characteristics**:
- **Accession**: GSE155179
- **Reference**: Zhang, Y.L., et al. (2020). Vitamin C enhances the number of ovarian follicles and fertility in aged female mice. Aging, 12(13), 13018.
- **Samples**: 12 MII oocytes
- **Age Groups**: 6 young (<30 years), 6 old (≥40 years)
- **Value**: Provides age contrast (young vs old) not available in current dataset

**Integration Steps**:

1. **Data Download**:
   ```python
   # Download expression matrix from GEO
   import GEOparse
   gse = GEOparse.get_GEO("GSE155179", destdir="./geo_data/")
   # Extract expression data from supplementary files
   ```

2. **Data Processing**:
   - Normalize expression data (TPM or counts)
   - Map gene identifiers to match Zenodo dataset
   - Quality control (filter low-quality cells/genes)

3. **Batch Correction**:
   - Apply scVI or Harmony for batch correction
   - Account for protocol differences (sequencing platform, library prep)

4. **Integration**:
   - Combine with existing Zenodo data
   - Verify stage labels (should all be MII)
   - Map age information

**Expected Impact**: +12 cells → n = 32 total

---

### GSE95477 Integration

**Dataset Characteristics**:
- **Accession**: GSE95477
- **Reference**: Reyes, J.M., et al. (2017). Differing molecular response of young and advanced maternal age human oocytes to IVM. Human Reproduction, 32(11), 2199-2208.
- **Samples**: 32 oocytes
- **Stages**: GV and MII (based on publication)
- **Age Range**: Estimated 25-35 years (from metadata)

**Integration Steps**:

1. **Data Download**:
   ```python
   gse = GEOparse.get_GEO("GSE95477", destdir="./geo_data/")
   # Extract expression data
   ```

2. **Data Processing**:
   - Normalize and map gene identifiers
   - Quality control
   - Verify stage annotations

3. **Batch Correction**:
   - Apply batch correction methods
   - Harmonize with existing data

4. **Integration**:
   - Combine with Zenodo + GSE155179
   - Update stage distribution
   - Map age information

**Expected Impact**: +32 cells → n = 64 total (after GSE155179 integration)

---

## Technical Challenges

### 1. Data Format Differences

**Challenge**: Different GEO datasets may use different:
- Expression units (FPKM, TPM, counts, log-transformed)
- Gene identifiers (Ensembl, RefSeq, Gene Symbol)
- Sequencing platforms (Illumina, other)

**Solution**:
- Standardize to TPM or counts
- Use biomaRt or similar for ID mapping
- Document all transformations

### 2. Batch Effects

**Challenge**: Different studies = different:
- Sequencing protocols
- Laboratories
- Sample preparation methods
- Time periods

**Solution**:
- Apply scVI (if available) or Harmony for batch correction
- Verify biological signal preserved after correction
- Use batch as covariate in downstream analysis

### 3. Stage Annotation Consistency

**Challenge**: Stage labels may differ between studies:
- GV vs VG (germinal vesicle)
- MI vs MII annotation quality

**Solution**:
- Standardize stage labels
- Verify using marker genes
- Document any discrepancies

### 4. Age Data Availability

**Challenge**: Age information may be:
- Missing for some samples
- In different formats (categorical vs continuous)
- From different sources (donor age vs oocyte age)

**Solution**:
- Extract from GEO metadata
- Standardize to continuous years
- Document missing data

---

## Implementation Timeline

### Phase 1: GSE155179 Integration (2-3 days)

**Day 1**:
- Download and inspect GSE155179 data
- Extract expression matrices
- Map gene identifiers
- Quality control

**Day 2**:
- Normalize expression data
- Apply batch correction
- Integrate with Zenodo data
- Verify integration quality

**Day 3**:
- Update analysis pipeline
- Re-run trajectory inference
- Validate results
- Update documentation

### Phase 2: GSE95477 Integration (2-3 days)

**Day 1**:
- Download and inspect GSE95477 data
- Extract expression matrices
- Map gene identifiers
- Quality control

**Day 2**:
- Normalize expression data
- Apply batch correction
- Integrate with existing combined data
- Verify integration quality

**Day 3**:
- Update analysis pipeline
- Re-run all analyses
- Validate results
- Update documentation

### Phase 3: Validation (1 day)

- Cross-study validation
- Compare results before/after integration
- Update all metrics and visualizations
- Final documentation

**Total Timeline**: 5-7 days

---

## Expected Outcomes

### Sample Size

**Before**: n = 20 cells
- 6 GV, 14 MI, partial MII
- Age range: 25-35 years

**After**: n = 64 cells (target)
- GV: 6 + X from GSE95477
- MI: 14 (from Zenodo)
- MII: X from GSE155179 + X from GSE95477
- Age range: <30 to ≥40 years (with GSE155179)

### Statistical Power

**Before**:
- Limited power for correlation analyses
- High variance in cross-validation
- Small test sets

**After**:
- Improved statistical power
- More stable cross-validation
- Larger test sets for validation
- Ability to stratify by age groups

### Validation Capabilities

**Before**:
- Single study (Zenodo)
- No cross-study validation possible

**After**:
- Multiple studies (Zenodo, GSE155179, GSE95477)
- Leave-one-study-out cross-validation
- Assessment of generalizability
- Age-stratified analysis

---

## Current Limitations

### Why Full Integration Not Yet Completed

1. **Time Constraints**: Full integration requires 5-7 days of focused work
2. **Data Processing Complexity**: Requires careful normalization and batch correction
3. **Validation Requirements**: Need to verify integration quality before using in analysis
4. **Current Focus**: Prioritizing validation metrics and statistical rigor first

### Workaround for Forum Presentation

**Current Approach**:
- Acknowledge sample size limitation honestly
- Emphasize strong validation framework despite small n
- Highlight cross-validation results
- Note plans for expansion in future work

**Presentation Statement**:
> "Current analysis includes 20 oocytes from a single study (Zenodo). While this limits statistical power, we have implemented comprehensive validation metrics including train/test splits and cross-validation. We are actively working to integrate additional datasets (GSE155179: 12 MII oocytes, GSE95477: 32 oocytes) to expand to n=64 cells, which will enable cross-study validation and improved statistical power."

---

## Files to Create/Update

### New Files
- `scripts/integrate_geo_datasets.py` - Integration script
- `data/integrated_expression_matrix.h5ad` - Combined AnnData object
- `docs/INTEGRATION_VALIDATION.md` - Integration quality metrics

### Updated Files
- `ADSPROJECT_new.ipynb` - Update data loading sections
- `RESULTS.md` - Update sample size and add cross-study validation
- `METRICS.md` - Add cross-study validation metrics
- `README.md` - Update dataset descriptions

---

## Success Criteria

### Integration Quality

- [ ] Expression data successfully downloaded and processed
- [ ] Gene identifiers mapped correctly (>90% overlap)
- [ ] Batch effects corrected (verify with batch mixing metrics)
- [ ] Biological signal preserved (stage separation maintained)

### Sample Size Targets

- [ ] Total n ≥ 50 cells
- [ ] Age range expanded (include <30 and ≥40 groups)
- [ ] Multiple stages represented
- [ ] Multiple studies integrated

### Validation

- [ ] Leave-one-study-out CV implemented
- [ ] Cross-study correlation r > 0.70
- [ ] Batch mixing metrics acceptable (kBET, LISI)
- [ ] Results consistent across studies

---

## References

1. Zhang, Y.L., et al. (2020). Vitamin C enhances the number of ovarian follicles and fertility in aged female mice. Aging, 12(13), 13018. [GSE155179]

2. Reyes, J.M., et al. (2017). Differing molecular response of young and advanced maternal age human oocytes to IVM. Human Reproduction, 32(11), 2199-2208. [GSE95477]

3. Lopez, R., et al. (2018). Deep generative modeling for single-cell transcriptomics. Nature Methods, 15(12), 1053-1058. [scVI for batch correction]

---

**Last Updated**: November 2025  
**Status**: Planning - Ready for Implementation

