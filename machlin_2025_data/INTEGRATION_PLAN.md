# Machlin et al. 2025 Dataset Integration Plan

**Source**: Zenodo DOI: 10.5281/zenodo.13224872  
**Paper**: "Single-cell analysis comparing early-stage oocytes from fresh and slow-frozen/thawed human ovarian cortex reveals minimal impact of cryopreservation on the oocyte transcriptome"

## Dataset Overview

### Files Downloaded

1. **4311_counts.txt** (29.4 MB)
   - ~60,629 genes (rows)
   - 192 samples (columns: Sample_4311-JM-1 through Sample_4311-JM-192)
   - Format: Tab-separated, gene-level counts
   - Columns: `gene_id`, `entrezgene_id`, `external_gene_name`, `description`, then sample columns
   - Gene IDs: Ensembl IDs (ENSG...)

2. **4311_metadata.csv** (19.8 KB)
   - 192 rows (one per sample)
   - Sequencing metadata only (barcodes, quality scores, etc.)
   - Missing: Biological metadata (age, stage, condition)

3. **3799_counts.txt** (16.9 MB)
   - Similar structure to 4311
   - ~60,629 genes
   - 96 samples (based on metadata file)

4. **3799_metadata.csv** (9.5 KB)
   - 96 rows
   - Sequencing metadata only

### Total Samples
- **4311 dataset**: 192 samples
- **3799 dataset**: 96 samples
- **Combined**: 288 samples (but paper mentions 144 oocytes)

**Note**: The discrepancy between 288 samples and 144 oocytes suggests:
- Possibly multiple technical replicates per biological sample
- Or different experimental conditions (fresh vs frozen)
- Need to consult paper for clarification

## Integration Steps Required

### 1. Metadata Enrichment
- [ ] Consult original paper for sample metadata (age, stage, fresh/frozen status)
- [ ] Create comprehensive metadata file with:
  - Donor age
  - Oocyte stage (GV, MI, MII)
  - Condition (fresh vs frozen/thawed)
  - Donor ID
  - Batch information

### 2. Data Preprocessing
- [ ] Load count matrices
- [ ] Match gene IDs to current dataset (Ensembl IDs should match)
- [ ] Combine 4311 and 3799 datasets if appropriate
- [ ] Filter cells/samples (QC similar to current pipeline)

### 3. Integration with Existing Data
- [ ] Merge with current 20-oocyte dataset
- [ ] Ensure consistent gene ID format
- [ ] Handle batch effects between datasets
- [ ] Apply scVI batch correction to combined dataset

### 4. Re-run Analysis Pipeline
- [ ] Trajectory learning (scVI + GPLVM)
- [ ] Pathway scoring
- [ ] Risk stratification
- [ ] Performance metrics

### 5. Validation
- [ ] Compare findings with current n=20 results
- [ ] Assess improvement in statistical power
- [ ] Update RESULTS.md with new sample size

## Current Status

**Completed**:
- Downloaded all 4 files from Zenodo
- Examined file structure
- Identified that data is at gene level (not transcript level)

**Remaining**:
- Need biological metadata (age, stage, condition)
- Need to determine sample organization (144 oocytes vs 288 samples)
- Need to integrate into pipeline
- Need to handle potential batch effects

## Notes

- Current dataset uses transcript-level (Kallisto) data, aggregated to genes
- Machlin dataset appears to be gene-level from the start
- Both use Ensembl IDs, so gene matching should be straightforward
- May need to check for gene ID version differences (e.g., ENSG00000000003 vs ENSG00000000003.15)

## Next Steps

1. **Immediate**: Read Machlin et al. 2025 paper to understand:
   - Sample structure (144 vs 288)
   - Available metadata
   - Experimental design (fresh vs frozen comparison)

2. **Short-term**: Create comprehensive metadata file

3. **Medium-term**: Modify pipeline to handle combined dataset

4. **Long-term**: Re-run full analysis and compare results

