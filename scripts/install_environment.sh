#!/bin/bash

# BIO559R Environment Installation Script
# This script automates the installation of the conda environment for BIO559R

set -e  # Exit on any error

echo "=========================================="
echo "BIO559R Environment Installation Script"
echo "=========================================="

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH."
    echo "Please install Miniconda first following the tutorial instructions."
    exit 1
fi

echo "✓ conda found"

# Check if environment.yml exists
if [ ! -f "environment.yml" ]; then
    echo "Error: environment.yml file not found in current directory."
    echo "Please make sure you are in the correct directory."
    exit 1
fi

echo "✓ environment.yml found"

# Create the conda environment
echo "Creating conda environment 'bio559r'..."
echo "This may take several minutes..."

conda env create -f environment.yml

if [ $? -eq 0 ]; then
    echo "✓ Environment created successfully"
else
    echo "✗ Environment creation failed"
    exit 1
fi

# Activate the environment and register R kernel
echo "Activating environment and registering R kernel..."

# Source conda to make conda activate work in script
eval "$(conda shell.bash hook)"
conda activate bio559r

# Register R kernel with Jupyter
echo "Registering R kernel with Jupyter..."
R --slave -e "IRkernel::installspec(user = FALSE)"

if [ $? -eq 0 ]; then
    echo "✓ R kernel registered successfully"
else
    echo "✗ R kernel registration failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "Installation completed successfully!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
echo "  conda activate bio559r"
echo ""
echo "To start Jupyter notebook, run:"
echo "  jupyter notebook"
echo ""
echo "You should now see both Python 3 and R kernels available in Jupyter."

