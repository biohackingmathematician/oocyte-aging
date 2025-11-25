# Activating the Environment

## Virtual Environment (Currently Deployed)

The environment has been deployed using Python 3.11.9 in a virtual environment.

### Activate the Environment

```bash
source venv/bin/activate
```

### Verify Installation

```bash
# Check scVI
python -c "import scvi; print('scVI version:', scvi.__version__)"

# Check other key packages
python -c "import scanpy as sc; print('scanpy version:', sc.__version__)"
python -c "import tensorflow as tf; print('TensorFlow version:', tf.__version__)"
python -c "import gpflow; print('GPflow version:', gpflow.__version__)"
```

### Deactivate the Environment

```bash
deactivate
```

## Installed Packages

All packages from `requirements.txt` have been successfully installed:
- scVI-tools: 1.4.0.post1
- scanpy: 1.11.5
- TensorFlow: 2.15.1
- GPflow: 2.10.0
- All other dependencies

## Using with Jupyter Notebook

If you want to use this environment with Jupyter notebooks:

```bash
source venv/bin/activate
pip install ipykernel
python -m ipykernel install --user --name=oocyte_analysis --display-name "Python (oocyte_analysis)"
```

Then select the kernel "Python (oocyte_analysis)" in your Jupyter notebook.

