# BIO559R: Single-Cell and Spatial Transcriptomics Analysis Tutorial

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)

A comprehensive tutorial for setting up and using computational tools for single-cell RNA-seq and spatial transcriptomics analysis in the BIO559R Introduction to Systems Biology course.

## 🎯 Overview

This tutorial provides step-by-step instructions for:

- Installing Ubuntu Linux for bioinformatics
- Setting up conda environments with Python and R
- Configuring Jupyter notebooks with R magic
- Installing and using the scverse ecosystem (scanpy, squidpy, etc.)
- Leveraging R packages for statistical analysis and visualization
- Integrating Python and R workflows seamlessly

## 📋 Prerequisites

- A computer with at least 8GB RAM (16GB recommended)
- 50GB+ free disk space
- Basic familiarity with command line interfaces
- No prior experience with Linux, Python, or R required

## 🚀 Quick Start

### Option 1: Automated Installation (Recommended)

1. **Clone this repository:**
   ```bash
   git clone https://github.com/your-username/BIO559R-Tutorial.git
   cd BIO559R-Tutorial
   ```

2. **Run the installation script:**
   ```bash
   chmod +x scripts/install_environment.sh
   ./scripts/install_environment.sh
   ```

3. **Test your installation:**
   ```bash
   conda activate bio559r
   python scripts/test_installation.py
   ```

4. **Start Jupyter notebook:**
   ```bash
   jupyter notebook
   ```

### Option 2: Manual Installation

Follow the detailed instructions in [`docs/BIO559R_tutorial.md`](docs/BIO559R_tutorial.md).

## 📁 Repository Structure

```
BIO559R-Tutorial/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── environment.yml              # Conda environment specification
├── .gitignore                   # Git ignore rules
├── CONTRIBUTING.md              # Contribution guidelines
├── docs/
│   └── BIO559R_tutorial.md     # Complete tutorial documentation
├── scripts/
│   ├── install_environment.sh  # Automated installation script
│   └── test_installation.py    # Installation verification script
├── notebooks/                   # Mathematical foundations notebooks
│   ├── README.md               # Notebook overview and instructions
│   ├── 01_ODEs_Mathematical_Foundations.ipynb
│   ├── 02_Linear_Algebra_Network_Theory.ipynb
│   ├── 03_Probability_Optimization_Steady_State.ipynb
│   └── example_notebook.ipynb  # Example analysis workflow
├── examples/
│   └── README.md               # Example datasets and analyses
└── assets/
    └── README.md               # Images and other assets
```

## 🛠 What's Included

### Python Packages (scverse ecosystem)
- **scanpy**: Single-cell analysis toolkit
- **anndata**: Annotated data structures
- **squidpy**: Spatial transcriptomics analysis
- **cellrank**: RNA velocity and trajectory inference
- **scvelo**: RNA velocity analysis
- And many more...

### R Packages
- **Seurat**: Single-cell genomics toolkit
- **ggplot2**: Advanced data visualization
- **dplyr**: Data manipulation
- **Bioconductor**: Bioinformatics packages
- **IRkernel**: R kernel for Jupyter

### Key Features
- ✅ Seamless Python-R integration via rpy2
- ✅ Jupyter notebook environment with both kernels
- ✅ Comprehensive quality control and preprocessing tools
- ✅ Advanced visualization capabilities
- ✅ Statistical analysis and modeling tools
- ✅ Spatial transcriptomics support

## 📖 Tutorial Contents

### Main Tutorial (docs/BIO559R_tutorial.md)
1. **Chapter 1**: Installing Ubuntu for Windows Users
2. **Chapter 2**: Setting Up Your Conda Environment
3. **Chapter 3**: Installing Python and R Packages
4. **Chapter 4**: Configuring Jupyter Notebooks for Python and R
5. **Chapter 5**: Introduction to scverse and Scanpy
6. **Chapter 6**: R for Visualization and Statistical Analysis
7. **Chapter 7**: Troubleshooting

### Mathematical Foundations Notebooks (notebooks/)
1. **01_ODEs_Mathematical_Foundations.ipynb**: Ordinary differential equations in systems biology
2. **02_Linear_Algebra_Network_Theory.ipynb**: Matrix operations and network analysis
3. **03_Probability_Optimization_Steady_State.ipynb**: Stochastic processes and optimization
4. **example_notebook.ipynb**: Integrated Python-R analysis workflow

## 🔬 Example Workflows

### Single-Cell Analysis (example_notebook.ipynb)
- Loading and preprocessing single-cell data with scanpy
- Quality control and filtering
- Dimensionality reduction (PCA, UMAP)
- Clustering with the Leiden algorithm
- Creating publication-quality plots with R/ggplot2
- Statistical analysis using R
- Seamless data transfer between Python and R

### Mathematical Foundations (notebooks/01-03)
- Modeling gene regulatory networks with ODEs
- Analyzing metabolic networks using linear algebra
- Stochastic simulations of molecular processes
- Parameter estimation and optimization
- Network topology analysis
- Steady-state and stability analysis

## 💻 System Requirements

### Minimum Requirements
- **OS**: Ubuntu 20.04+ (or Windows with WSL2)
- **RAM**: 8GB
- **Storage**: 50GB free space
- **CPU**: 2+ cores

### Recommended Requirements
- **OS**: Ubuntu 22.04+
- **RAM**: 16GB+
- **Storage**: 100GB+ free space
- **CPU**: 4+ cores
- **GPU**: Optional, for accelerated computing

## 🆘 Getting Help

### Common Issues
- Check the [Troubleshooting section](docs/BIO559R_tutorial.md#chapter-7-troubleshooting) in the tutorial
- Run the test script: `python scripts/test_installation.py`
- Ensure you're in the correct conda environment: `conda activate bio559r`

### Resources
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [scverse Ecosystem](https://scverse.org/)
- [Seurat Documentation](https://satijalab.org/seurat/)
- [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/)

## 👥 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Instructor**: Tommy W. Terooatea
- **Course**: BIO559R - Introduction to Systems Biology
- **scverse Community**: For developing excellent single-cell analysis tools
- **R/Bioconductor Community**: For statistical and visualization packages

## 📞 Contact

For questions about this tutorial:
- Course-related: Contact your instructor
- Technical issues: Open an issue on GitHub
- General questions: Check the documentation first

---

**Happy analyzing! 🧬📊**

