# Raw EDA Visualizations Documentation

**Date**: November 2025  
**Purpose**: Comprehensive exploratory data analysis of raw, untransformed expression data

---

## Overview

This document describes the raw exploratory data analysis (EDA) visualizations generated from untransformed expression data. These plots provide insights into data quality, distribution characteristics, and basic patterns before any normalization or transformation steps.

---

## Generated Visualizations

### 1. Library Size Distribution

**File**: `visualizations/raw_eda_library_size.png`

**Description**: 
- **Left panel**: Histogram of total library size (sum of all expression values) across all samples
- **Right panel**: Boxplot of library size stratified by developmental stage (GV, MI, MII)

**Interpretation**:
- Library size indicates total sequencing depth per sample
- Large variation may indicate technical artifacts or biological differences
- Stage-specific differences may reflect developmental stage characteristics

**Key Metrics**:
- Mean and median library sizes are displayed
- Sample counts per stage are annotated

---

### 2. Mitochondrial Gene Percentage

**File**: `visualizations/raw_eda_mitochondrial.png`

**Description**:
- **Left panel**: Histogram of mitochondrial gene expression percentage across all samples
- **Right panel**: Boxplot of mitochondrial percentage by developmental stage

**Interpretation**:
- High mitochondrial percentage (>10-20%) may indicate cell stress or low-quality samples
- Mitochondrial genes are identified using multiple patterns:
  - MT- prefix (human mitochondrial genes)
  - Common mitochondrial gene symbols (COX, ND, ATP, etc.)
  - Transcript IDs containing MT patterns

**Biological Relevance**:
- Mitochondrial function is critical for oocyte quality
- Changes in mitochondrial gene expression may reflect aging or quality decline

---

### 3. Top Expressed Genes by Stage

**File**: `visualizations/raw_eda_top_genes_by_stage.png`

**Description**:
- Horizontal bar charts showing the top 20 most highly expressed genes for each developmental stage
- Genes ranked by mean expression across samples within each stage

**Interpretation**:
- Identifies stage-specific highly expressed genes
- Reveals biological differences between GV, MI, and MII stages
- Helps identify potential marker genes for each stage

**Use Cases**:
- Quality control: Verify expected stage-specific markers are highly expressed
- Biological insight: Identify genes driving stage-specific biology

---

### 4. Gene Expression Heatmap

**File**: `visualizations/raw_eda_heatmap.png`

**Description**:
- Heatmap of the top 50 most variable genes across all samples
- Samples ordered by developmental stage (if metadata available)
- Expression values log-transformed for visualization (log(expression + 1))

**Interpretation**:
- Visualizes expression patterns across samples
- Identifies sample clusters and outliers
- Reveals stage-specific expression patterns
- Variable genes are selected based on coefficient of variation (CV)

**Technical Details**:
- Color scale: viridis colormap
- Log transformation applied only for visualization (raw data preserved)
- Gene and sample names truncated for readability

---

### 5. Gene-Gene Correlation Matrix

**File**: `visualizations/raw_eda_correlation_matrix.png`

**Description**:
- Correlation matrix for the top 30 most expressed genes
- Pearson correlation coefficients between all pairs of genes
- Color-coded heatmap (red = positive correlation, blue = negative correlation)

**Interpretation**:
- Identifies co-expressed gene modules
- Reveals functional relationships between genes
- High correlations may indicate:
  - Genes in the same pathway
  - Coordinated regulation
  - Technical artifacts (if correlations are unusually high)

**Use Cases**:
- Pathway analysis: Identify functionally related genes
- Quality control: Detect technical artifacts (unusually high correlations)
- Network analysis: Identify gene modules

---

### 6. Raw Counts Boxplots

**File**: `visualizations/raw_eda_boxplots.png`

**Description**:
- Boxplots of raw expression values for the top 20 most expressed genes
- Stratified by developmental stage (if metadata available)
- Each subplot shows one gene's expression distribution

**Interpretation**:
- Visualizes expression distribution for individual genes
- Reveals stage-specific expression differences
- Identifies highly variable genes
- Shows outliers and expression ranges

**Key Features**:
- Median, quartiles, and outliers displayed
- Stage-specific comparisons enable identification of differentially expressed genes
- Useful for quality control and biological interpretation

---

## Data Characteristics

### Expression Data

- **Source**: Kallisto abundance.tsv files (TPM values)
- **Total Genes**: 204,563 transcripts
- **Total Samples**: 20 oocytes
- **Expression Range**: 0 - 113,029 TPM
- **Data Type**: Transcript-level expression (TPM)

### Sample Distribution

- **GV (Germinal Vesicle)**: 6 samples
- **MI (Metaphase I)**: 14 samples
- **MII (Metaphase II)**: Variable representation

---

## Quality Control Insights

### Library Size

- **Purpose**: Assess sequencing depth and sample quality
- **Expected Range**: Varies by sequencing protocol
- **Interpretation**: 
  - Very low library sizes may indicate failed samples
  - Large variation may require normalization
  - Stage-specific differences may be biological

### Mitochondrial Percentage

- **Purpose**: Assess cell quality and stress
- **Expected Range**: <10-20% for healthy cells
- **High MT% Indicates**:
  - Cell stress or damage
  - Low-quality samples
  - Potential contamination

### Expression Distribution

- **Purpose**: Assess data quality and identify artifacts
- **Expected Patterns**:
  - Most genes have low expression (typical for RNA-seq)
  - Few highly expressed genes (housekeeping, stage-specific markers)
  - Right-skewed distribution

---

## Biological Insights

### Stage-Specific Patterns

1. **GV Stage**:
   - Expected high expression of early meiotic genes
   - Lower expression of maturation-specific genes

2. **MI Stage**:
   - Transition genes activated
   - Cell cycle and spindle assembly genes upregulated

3. **MII Stage**:
   - Maturation completion genes
   - Chromosome segregation genes

### Gene Expression Patterns

- **Highly Expressed Genes**: Likely include housekeeping genes and stage-specific markers
- **Variable Genes**: May represent biological variation or technical noise
- **Correlated Genes**: May indicate functional modules or pathways

---

## Technical Notes

### Data Loading

The script attempts to load data in the following order:
1. AnnData object (`.h5ad` file) with counts layer
2. Individual `abundance.tsv` files from kallisto directories

### Gene Identification

- **Transcript IDs**: Original data uses transcript IDs (e.g., ENST...)
- **Gene Mapping**: May require additional mapping to gene symbols for biological interpretation
- **Mitochondrial Detection**: Uses multiple patterns to identify MT genes from transcript IDs

### Visualization Parameters

- **Top Genes**: Configurable (default: 20-50 genes)
- **Color Schemes**: Optimized for publication
- **Figure Size**: Optimized for clarity and readability
- **DPI**: 300 for publication quality

---

## Usage

### Generate All Plots

```bash
python scripts/generate_raw_eda.py
```

### Output Location

All plots are saved to: `visualizations/raw_eda_*.png`

### Customization

Edit the script to modify:
- Number of genes displayed
- Color schemes
- Figure sizes
- Output directory

---

## Integration with Analysis Pipeline

These raw EDA plots should be generated:

1. **Before normalization**: To assess raw data quality
2. **After quality control**: To verify filtering effects
3. **For publication**: As supplementary figures showing data characteristics

### Recommended Workflow

1. Generate raw EDA plots (this script)
2. Perform quality control filtering
3. Generate post-QC EDA plots
4. Compare before/after to assess filtering impact

---

## References

- Kallisto: Bray, N.L., et al. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34(5), 525-527.
- TPM: Transcripts Per Million - standard normalization for RNA-seq data

---

**Last Updated**: November 2025  
**Status**: Complete

