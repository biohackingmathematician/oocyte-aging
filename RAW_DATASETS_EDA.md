# Raw Datasets EDA: GEO and Zenodo

## Overview
Exploratory data analysis of the raw datasets used in the project, showing sample distribution, stage breakdown, and age information.

**File**: `raw_datasets_eda.png`

---

## Panel A: Dataset Sources Overview

Shows the number of samples from each dataset source:

- **Zenodo (10.5281/zenodo.14163313)**: 20 samples
  - Single-cell RNA-seq data
  - Quantified with Kallisto
  - Mean sequencing depth: ~39.6M reads per sample

- **GSE155179 (GEO)**: 12 samples
  - Age metadata (30-40 years range)
  - Used for age integration

- **GSE95477 (GEO)**: 32 samples
  - Age metadata (25-35 years range)
  - Used for age integration

**Total**: 64 samples across all datasets

---

## Panel B: Stage Distribution (Zenodo)

Shows the distribution of developmental stages in the Zenodo dataset:

- **GV (Germinal Vesicle)**: 6 samples
- **MI (Metaphase I)**: 14 samples
- **MII (Metaphase II)**: 0 samples (partial representation)

**Note**: Stage information is available for the 20 Zenodo samples used in transcriptomic analysis.

---

## Panel C: Age Distribution (GEO Datasets)

Histogram showing chronological age distribution from GEO datasets:

- **Age Range**: 25-35 years (mean: 32.0 years)
- **Age Groups**:
  - Young (<30): 6 samples (30%)
  - Middle (30-35): 14 samples (70%)

**Source**: Age metadata extracted from GSE155179 and GSE95477 GEO datasets and mapped to Zenodo samples.

---

## Panel D: Sample Distribution by Stage Across Datasets

Grouped bar chart showing sample counts by stage and dataset:

- Visualizes how samples are distributed across developmental stages
- Shows which datasets contribute to each stage
- Helps identify data availability for different comparisons

---

## Panel E: Dataset Summary Statistics

Text summary panel with key statistics:

```
RAW DATASETS SUMMARY

ZENODO (10.5281/zenodo.14163313):
  • Total samples: 20
  • Stages: GV (6), MI (14)
  • Technology: Single-cell RNA-seq
  • Quantification: Kallisto

GEO DATASETS:
  • GSE155179: 12 samples
  • GSE95477: 32 samples
  • Data type: Age metadata

COMBINED:
  • Total samples: 64
  • Age range: 25-35 years
  • Mean age: 32.0 years
```

---

## Key Findings

### 1. Dataset Sources
- **Primary transcriptomic data**: Zenodo (20 oocytes)
- **Age metadata**: GEO datasets (44 total samples with age info)
- **Integration**: Age data mapped to 20/20 Zenodo samples

### 2. Stage Distribution
- **Imbalanced**: More MI (14) than GV (6) samples
- **No MII**: MII stage not represented in current dataset
- **Clinical relevance**: Focus on GV→MI transition (critical quality decline)

### 3. Age Distribution
- **Limited range**: 25-35 years (reproductive prime)
- **Missing cohorts**: No samples <25 or >35 years
- **Clinical context**: Represents typical fertility preservation age range

### 4. Sequencing Quality
- **Deep sequencing**: ~39.6M reads per sample (high quality)
- **Technology**: Single-cell RNA-seq with Kallisto quantification
- **Coverage**: Sufficient for transcriptomic analysis

---

## Data Integration Strategy

1. **Transcriptomic Data**: 20 oocytes from Zenodo
   - Raw counts from Kallisto quantification
   - 204,563 genes (reduced to 2,000 most variable)

2. **Age Data**: Integrated from GEO datasets
   - GSE155179: 12 samples (30-40 years)
   - GSE95477: 32 samples (25-35 years)
   - Age mapped to 20/20 Zenodo samples

3. **Stage Information**: From Zenodo metadata
   - GV: 6 samples
   - MI: 14 samples

---

## Limitations

1. **Sample Size**: 20 oocytes (small but sufficient for proof-of-concept)
2. **Stage Imbalance**: More MI than GV (14 vs 6)
3. **Age Range**: Limited to 25-35 years (missing extremes)
4. **Single Study**: All transcriptomic data from one Zenodo study
5. **No MII**: MII stage not represented (incomplete trajectory)

---

## How to Generate

```bash
python3 create_raw_datasets_eda.py
```

**Requirements**:
- `sample_metadata_with_age.csv` - Sample metadata
- `clinical_decision_framework_final.csv` - Age data (optional)
- `geo_data/` - GEO dataset files (optional, for metadata)
- `zenodo_data/final_code/kallisto/` - Zenodo kallisto output files

---

## References

1. **Zenodo Dataset**: 
   - DOI: 10.5281/zenodo.14163313
   - Llonch, S., et al. (2021). Single human oocyte transcriptome analysis reveals distinct maturation stage-dependent pathways impacted by age. Aging Cell, 20(5), e13360.

2. **GEO Datasets**:
   - **GSE155179**: Zhang, Y.L., et al. (2020). Vitamin C enhances the number of ovarian follicles and fertility in aged female mice. Aging, 12(13), 13018.
   - **GSE95477**: Reyes, J.M., et al. (2017). Differing molecular response of young and advanced maternal age human oocytes to IVM. Human Reproduction, 32(11), 2199-2208.

---

## Visualization Details

- **Format**: Multi-panel figure (2×3 GridSpec)
- **Resolution**: 300 DPI
- **Size**: 16×10 inches
- **File**: `raw_datasets_eda.png`
- **Color Scheme**: Consistent with project (GV=green, MI=red, MII=blue)

This EDA provides a comprehensive overview of the raw datasets before processing, helping understand data availability, quality, and limitations for the analysis pipeline.

