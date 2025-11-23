# Installation Guide for Upgrade Sections

## Python Version Compatibility Issue

**Current Python Version**: 3.14.0

**Problem**: Some required packages don't support Python 3.14 yet:
- `scvi-tools` requires Python <3.14 (supports 3.10-3.13)
- `tensorflow` doesn't support Python 3.14
- `numba` (dependency) doesn't support Python 3.14

## Solutions

### Option 1: Use Python 3.10-3.13 (Recommended)

Create a conda environment with compatible Python version:

```bash
# Create conda environment with Python 3.11
conda create -n oocyte_analysis python=3.11
conda activate oocyte_analysis

# Install packages
pip install scvi-tools gpflow tensorflow scanpy scikit-learn pandas numpy matplotlib GEOparse
```

### Option 2: Use Fallback Methods (Current Setup)

The notebook is designed to work with fallback methods:
- **scVI → PCA**: If scVI fails, uses PCA for dimensionality reduction
- **GPLVM → Simplified**: If tensorflow/gpflow unavailable, uses simplified trajectory
- **All sections have error handling**: Missing packages won't crash the notebook

### Option 3: Use Docker/Container

Use a pre-configured environment:
```bash
# Example with Jupyter Docker
docker run -p 8888:8888 jupyter/scipy-notebook:latest
```

## Package Installation Commands

### For Python 3.10-3.13:
```bash
pip install scvi-tools gpflow tensorflow scanpy scikit-learn pandas numpy matplotlib GEOparse
```

### For Python 3.14 (with limitations):
```bash
# These may work:
pip install scanpy scikit-learn pandas numpy matplotlib GEOparse

# These will fail:
# pip install scvi-tools  # Requires Python <3.14
# pip install tensorflow   # Not available for 3.14
# pip install gpflow       # May have issues
```

## What Still Works Without Full Installation

Even without scVI and tensorflow, you can still run:
1.  **Section 1**: Age Data Integration (uses GEOparse)
2.  **Section 2**: Will use PCA instead of scVI
3.  **Section 3**: GPLVM will use simplified version
4.  **Section 4**: AMH calibration may be limited
5.  **Section 5**: Risk Stratification (uses sklearn)
6.  **Section 6**: Cross-validation (uses sklearn)
7.  **Section 7**: Results integration

## Recommended Approach

1. **For full functionality**: Use Python 3.11 or 3.12 with conda
2. **For testing**: Use current setup with fallbacks
3. **For production**: Set up proper environment with all dependencies

## Checking Your Setup

Run this in Python to check what's available:
```python
import sys
print(f"Python: {sys.version}")

packages = ['scvi', 'gpflow', 'tensorflow', 'scanpy', 'sklearn']
for pkg in packages:
    try:
        __import__(pkg)
        print(f" {pkg}")
    except ImportError:
        print(f" {pkg} - MISSING")
```

