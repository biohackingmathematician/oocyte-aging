# Oocyte Aging Analysis: Bayesian GPLVM for Fertility Preservation Timing

**Authors**: Agna Chan, Aniqa Nayim, Rimjhim Singh  
**Course**: STAT GR5243 Applied Data Science  
**Institution**: Columbia University  
**Date**: November 2025

---

## Project Overview

This uses a multi-dimensional Bayesian generative model to quantify variability and uncertainty in oocyte aging trajectories from single-cell transcriptomics data. The goal is to identify intervention windows for fertility preservation by combining:

- **scVI** for batch correction across studies
- **Bayesian GPLVM** for uncertainty-aware trajectory analysis
- **Clinical risk classification** for personalized fertility preservation
- **AMH calibration** for population-level fertility predictions

---

## Repository Structure

```
.
├── ADSPROJECT_new.ipynb       # Main analysis notebook
├── ADSPROJECT_base.ipynb     # Base notebook (reference)
├── README.md                  # This file
├── RESULTS.md                 # Detailed results summary
├── METRICS.md                 # Validation metrics and performance evaluation
├── LICENSE                    # MIT License
├── environment.yml            # Conda environment file (for scVI support)
├── requirements.txt           # Python requirements (full installation)
├── requirements-minimal.txt   # Minimal requirements (no scVI/GPLVM)
├── setup_environment.sh       # Automated environment setup script
│
├── docs/                      # Documentation
│   ├── INSTALLATION_GUIDE.md
│   ├── QUICK_START.md
│   ├── METRICS_EVALUATION.md  # Comprehensive metrics evaluation
│   ├── RAW_DATASETS_EDA.md
│   ├── FORUM_EDA_VISUALIZATIONS.md
│   └── *.pdf                  # Reference papers and reports
│
├── scripts/                   # Analysis and visualization scripts
│   ├── calculate_metrics.py  # Calculate validation metrics
│   ├── create_combined_intervention_plot.py
│   ├── create_complete_results_summary.py
│   ├── create_forum_eda_visualizations.py
│   ├── create_forum_raw_eda.py
│   ├── create_raw_datasets_eda.py
│   ├── run_all_visualizations.sh
│   └── run_visualizations_from_csv.py
│
├── visualizations/            # Generated figures
│   ├── combined_intervention_plot.png
│   ├── complete_results_summary.png
│   ├── gplvm_trajectory_analysis.png
│   ├── risk_stratification.png
│   ├── metrics_roc_curve.png
│   ├── metrics_pr_curve.png
│   └── [other visualization files]
│
├── data/                      # Data files
│   ├── clinical_decision_framework_final.csv
│   ├── sample_metadata.csv
│   └── sample_metadata_with_age.csv
│
├── zenodo_data/              # Zenodo dataset (20 oocytes)
│   └── final_code/
│       ├── kallisto/         # Kallisto abundance files
│       ├── DeSeq2/           # Differential expression results
│       ├── EdgeR/            # EdgeR analysis results
│       └── [other analysis files]
│
└── geo_data/                 # GEO datasets
    ├── GSE155179_family.soft.gz
    └── GSE95477_family.soft.gz
```

---

## Requirements

### Python Version

**Recommended**: Python 3.10-3.13 (for full functionality)  
**Current**: Python 3.14 (with fallback methods)

**Note**: Some packages (scvi-tools, tensorflow, gpflow) require Python <3.14. The notebook includes fallback methods for Python 3.14 compatibility.

### Required Packages

#### Core Packages (Required)
```bash
pip install pandas numpy matplotlib scipy scikit-learn
```

#### For Full Functionality (Python 3.10-3.13)
```bash
pip install scvi-tools gpflow tensorflow scanpy anndata GEOparse
```

#### For Basic Functionality (Python 3.14)
```bash
pip install pandas numpy matplotlib scipy scikit-learn GEOparse
```

### Installation Instructions

#### Option 1: Automated Setup Script (Easiest)

```bash
# Run the setup script
./setup_environment.sh

# Activate the environment
conda activate oocyte_analysis

# Verify scVI installation
python -c "import scvi; print('scVI version:', scvi.__version__)"
```

#### Option 2: Using Conda with environment.yml (Recommended for scVI)

```bash
# Create environment from file (includes Python 3.11 and all dependencies)
conda env create -f environment.yml

# Activate environment
conda activate oocyte_analysis

# Verify scVI installation
python -c "import scvi; print('scVI version:', scvi.__version__)"
```

#### Option 3: Using Conda Manually

```bash
# Create environment with Python 3.11
conda create -n oocyte_analysis python=3.11
conda activate oocyte_analysis

# Install from requirements.txt
pip install -r requirements.txt
```

#### Option 4: Using pip with requirements.txt

```bash
# Requires Python 3.10-3.13 for full functionality
# Create virtual environment first
python3.11 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install all packages
pip install -r requirements.txt
```

### Verifying Installation

After setting up the environment, verify that scVI and other key packages are installed:

```bash
conda activate oocyte_analysis

# Check scVI
python -c "import scvi; print('scVI version:', scvi.__version__)"

# Check scanpy
python -c "import scanpy as sc; print('scanpy version:', sc.__version__)"

# Check tensorflow
python -c "import tensorflow as tf; print('TensorFlow version:', tf.__version__)"

# Check gpflow
python -c "import gpflow; print('GPflow version:', gpflow.__version__)"
```

### Troubleshooting Installation

#### Issue: "ImportError: cannot import name 'X' from 'scvi'"
**Solution**: Ensure you're using Python 3.10-3.13. Check with `python --version`. If using Python 3.14, use the minimal requirements or switch to the conda environment.

#### Issue: "No module named 'tensorflow'"
**Solution**: TensorFlow may have compatibility issues. Try installing specific version:
```bash
pip install tensorflow==2.13.0
```

#### Issue: Conda environment creation fails
**Solution**: Update conda and try again:
```bash
conda update conda
conda env create -f environment.yml
```

#### Option 5: Minimal Installation (Python 3.14, without scVI/GPLVM)

```bash
# For Python 3.14 with fallback methods only
pip install -r requirements-minimal.txt

# Note: scvi-tools, tensorflow, gpflow will not install on Python 3.14
# The notebook will use fallback methods automatically
```

---

## Data

### Datasets Used

1. **Zenodo Dataset** (Primary)
   - **Source**: Zenodo 14163313
   - **Full Dataset**: 72 oocytes across all maturation stages (18–43 years)
   - **Analysis Subset**: 20 oocytes (6 GV, 14 MI) with complete age + stage metadata
   - **Rationale for Subsetting**: 
     - Focus on the GV→MI transition (primary maturation checkpoint)
     - Keep the project computationally tractable
     - Ensure complete metadata availability for age integration and clinical analysis
   - **Raw Features**: 204,563 Kallisto transcripts (isoforms)
   - **Gene-level Matrix**: 126,966 gene symbols (after transcript-to-gene mapping)
   - **Detected per Cell**: ~4,256 ± 892 genes (after QC, ≥2 cells expressing)
   - **HVGs Used for Modeling**: 2,000 highly variable genes
   - **Format**: Kallisto abundance files (TPM)
   - **Location**: `zenodo_data/final_code/kallisto/`
   
   **Note**: The dataset contains transcript-level abundance estimates. We collapsed transcripts to gene symbols for downstream analysis, resulting in 126,966 unique genes from 204,563 transcripts.

2. **GEO Datasets** (Age Data)
   - **GSE155179**: 12 samples with age data (30-40 years)
   - **GSE95477**: 32 samples with age data
   - **Purpose**: Chronological age integration
   - **Location**: `geo_data/` (downloaded automatically)

### Data Download

The notebook automatically downloads the Zenodo dataset. GEO datasets are downloaded on first run.

---

## Usage

### Running the Notebook

1. **Start Jupyter Notebook**:
   ```bash
   jupyter notebook ADSPROJECT_new.ipynb
   ```

2. **Run Cells Sequentially**:
   - Cells 1-64: Data loading and preprocessing
   - Cell 65: Package installation check
   - Cell 66: Pre-flight validation
   - Cells 67-81: Upgrade sections (7 major analyses)

3. **Expected Runtime**:
   - Data loading: ~5-10 minutes
   - Age integration: ~2-3 minutes
   - scVI (if available): ~10-15 minutes
   - GPLVM: ~5-10 minutes
   - Risk stratification: ~1-2 minutes
   - **Total**: ~30-45 minutes

### Key Sections

#### Section 1: Age Data Integration (Cell 69)
- Parses GEO datasets for age information
- Maps age to oocyte samples
- Creates age groups

#### Section 2: scVI Batch Correction (Cell 71)
- Corrects for batch effects across studies
- Falls back to PCA if scVI unavailable

#### Section 3: Bayesian GPLVM (Cell 73)
- Learns 1D cellular age trajectory
- Quantifies uncertainty
- Falls back to simplified trajectory if tensorflow unavailable

#### Section 4: AMH Calibration (Cell 75)
- Maps cellular age to AMH predictions
- Requires gpflow (skipped if unavailable)

#### Section 5: Risk Stratification (Cell 77)
- Classifies oocytes into risk groups
- Low/Moderate/High risk categories

#### Section 6: Cross-Study Validation (Cell 79)
- Leave-one-study-out validation
- Requires multiple studies

#### Section 7: Final Results Integration (Cell 81)
- Compiles all results
- Generates summary visualizations
- Saves output files

---

## Output Files

### Generated Files

1. **`clinical_decision_framework_final.csv`**
   - Per-cell predictions and risk assessments
   - Columns: age, cellular_age_z, uncertainty, risk_group, risk_score

2. **`gplvm_trajectory_analysis.png`**
   - GPLVM trajectory visualization
   - Shows cellular age and uncertainty in UMAP space

3. **`risk_stratification.png`**
   - Risk group distribution
   - Features heatmap by risk group

4. **`complete_results_summary.png`**
   - Comprehensive summary figure
   - All analyses in one visualization

### Data Files

- Processed AnnData object (saved in notebook)
- Trajectory coordinates
- Gene correlation results

---

## Troubleshooting

### Common Issues

#### Issue: "ModuleNotFoundError: No module named 'scvi'"
**Solution**: 
- If using Python 3.14: This is expected. The notebook will use PCA fallback.
- If using Python 3.10-3.13: Install with `pip install scvi-tools`

#### Issue: "ModuleNotFoundError: No module named 'gpflow'"
**Solution**:
- If using Python 3.14: This is expected. AMH calibration will be skipped.
- If using Python 3.10-3.13: Install with `pip install gpflow tensorflow`

#### Issue: "NameError: name 'adata' is not defined"
**Solution**: Run data loading cells (1-64) before upgrade sections.

#### Issue: GEO datasets not downloading
**Solution**: Check internet connection. Files are cached in `geo_data/` after first download.

### Python 3.14 Compatibility

The notebook is designed to work with Python 3.14 using fallback methods:
- **scVI → PCA**: Uses PCA for dimensionality reduction
- **GPLVM → Simplified**: Uses simplified trajectory analysis
- **AMH → Skipped**: AMH calibration requires gpflow

For full functionality, use Python 3.10-3.13.

---

## Results Summary

### Key Findings

1. **Age Integration**: Successfully integrated age data from 2 GEO datasets (20/20 cells)
2. **Cellular Aging**: Computed cellular age trajectory with uncertainty (r=0.27 with chronological age)
3. **Risk Stratification**: Identified 3 risk groups:
   - 65% Low Risk (Resilient Agers)
   - 30% Moderate Risk
   - 5% High Risk (Accelerated Agers)
4. **Health Score**: Strong correlation with maturation trajectory (r=-0.79, p<0.001)

### Model Performance Metrics

**Risk Stratification Performance**:
- AUC-ROC: 1.000 (perfect separation of high-risk group)
- Precision-Recall AUC: 1.000

**Clinical Health Score Validation**:
- Discriminative ability (GV vs MI): AUC = 1.000
- Effect size (Cohen's d): 2.161 (large effect)
- Mann-Whitney U test: p < 0.001

**Correlation Analysis**:
- Cellular age Z vs Health Score: r = -0.837, p < 0.001 (strong negative correlation)
- Cellular age Z vs Chronological Age: r = 0.270, p = 0.249 (moderate, not significant)

**Latent Space Quality**:
- Silhouette score (GV vs MI separation): -0.139
- Davies-Bouldin index: 1.422
- Calinski-Harabasz score: 1.417

See `RESULTS.md` for detailed results and `METRICS_EVALUATION.md` for comprehensive metrics evaluation and recommendations.

---

## Citation

If you use this code, please cite:

```
Chan, A., Nayim, A., & Singh, R. (2025). Evaluating Oocyte Aging Uncertainty: 
A Multidimensional Bayesian Approach for Personalized Fertility Preservation Timing. 
STAT GR5243 Applied Data Science, Columbia University.
```

### Data Sources

1. **Zenodo**: Llonch, S., et al. (2021). Single human oocyte transcriptome analysis. 
   Aging Cell, 20(5), e13360. DOI: 10.5281/zenodo.14163313

2. **GSE155179**: Zhang, Y.L., et al. (2020). Vitamin C enhances ovarian follicles. 
   Aging, 12(13), 13018.

3. **GSE95477**: Reyes, J.M., et al. (2017). Differing molecular response of oocytes. 
   Human Reproduction, 32(11), 2199-2208.

---

## Contact

For questions or issues, please contact:
- **Agna Chan**: [email]
- **Aniqa Nayim**: [email]
- **Rimjhim Singh**: [email]

---

## License

This is for academic use only. Please respect data usage terms from Zenodo and GEO datasets.

---

## Acknowledgments

- Columbia University STAT GR5243 Applied Data Science course
- Zenodo and GEO for providing datasets
- scVI, GPflow, and scanpy developers for excellent tools

---

**Last Updated**: November 18, 2025

