#!/bin/bash
# Setup script for oocyte analysis environment with scVI support

set -e  # Exit on error

echo "Setting up oocyte analysis environment with scVI support..."
echo ""

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed. Please install Miniconda or Anaconda first."
    echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if environment already exists
if conda env list | grep -q "^oocyte_analysis "; then
    echo "Environment 'oocyte_analysis' already exists."
    read -p "Do you want to remove it and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n oocyte_analysis
    else
        echo "Using existing environment. Activate it with: conda activate oocyte_analysis"
        exit 0
    fi
fi

# Create environment from yml file
echo "Creating conda environment from environment.yml..."
conda env create -f environment.yml

# Activate environment
echo ""
echo "Environment created successfully!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate oocyte_analysis"
echo ""
echo "To verify scVI installation, run:"
echo "  python -c 'import scvi; print(\"scVI version:\", scvi.__version__)'"
echo ""

