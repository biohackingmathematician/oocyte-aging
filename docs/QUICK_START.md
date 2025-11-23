# Quick Start Guide

## Current Situation

You're running **Python 3.14**, which is too new for some packages:
-  `scvi-tools` - Not compatible (needs Python <3.14)
-  `tensorflow` - Not compatible (needs Python <3.14)
-  `gpflow` - May have issues

## What You Can Do Now

### Option 1: Run with Fallbacks (Easiest)

The notebook is designed to work even without scVI and tensorflow:

1. **Run Cell 65**: Package installation (will show warnings, that's OK)
2. **Run Cell 66**: Pre-flight checks
3. **Run Sections 1-7**: They will use fallback methods:
   - Section 2: Uses PCA instead of scVI
   - Section 3: Uses simplified GPLVM
   - Other sections work normally

### Option 2: Use Compatible Python Version (Best for Full Features)

```bash
# Install conda if you don't have it
# Then create environment:
conda create -n oocyte python=3.11
conda activate oocyte
pip install scvi-tools gpflow tensorflow scanpy
```

### Option 3: Test What Works

Run this to see what's available:
```python
import sys
print(f"Python: {sys.version_info.major}.{sys.version_info.minor}")

for pkg in ['scanpy', 'sklearn', 'pandas', 'numpy', 'GEOparse', 'scvi', 'gpflow', 'tensorflow']:
    try:
        __import__(pkg)
        print(f" {pkg}")
    except:
        print(f" {pkg}")
```

## Running the Notebook

1. **Start from Cell 65** (package check)
2. **Continue through all sections** - errors are handled gracefully
3. **Check output files** - even with fallbacks, you'll get results

## Expected Behavior

-  Sections 1, 5, 6, 7: Will work fully
-  Section 2: Will use PCA fallback (still functional)
-  Section 3: Will use simplified GPLVM (still functional)
-  Section 4: May have limited functionality

All sections save outputs and provide meaningful results even with fallbacks!

