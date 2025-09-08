# Mathematical Foundations Notebooks - BIO559R

This directory contains comprehensive Jupyter notebooks covering the mathematical foundations of systems biology, designed for **BIO559R - Introduction to Systems Biology** taught by Tommy W. Terooatea.

## 📚 Notebook Overview

### 01_ODEs_Mathematical_Foundations.ipynb
**Ordinary Differential Equations in Systems Biology**

- **Duration**: 3-4 hours
- **Prerequisites**: Basic calculus, Python fundamentals
- **Topics Covered**:
  - First-order ODEs and biological applications
  - Systems of ODEs for multi-component models
  - Phase plane analysis and stability
  - Bistable switches and cellular decision-making
  - Oscillatory dynamics (repressilator, circadian rhythms)
  - Numerical methods and parameter sensitivity

### 02_Linear_Algebra_Network_Theory.ipynb
**Linear Algebra and Network Analysis**

- **Duration**: 3-4 hours
- **Prerequisites**: Basic linear algebra, matrix operations
- **Topics Covered**:
  - Matrix operations in biological contexts
  - Stoichiometric matrices and flux balance analysis
  - Principal Component Analysis for omics data
  - Network topology and graph theory
  - Centrality measures and hub identification
  - Network motifs and biological significance

### 03_Probability_Optimization_Steady_State.ipynb
**Probability, Optimization, and Steady-State Analysis**

- **Duration**: 3-4 hours
- **Prerequisites**: Basic probability, statistics
- **Topics Covered**:
  - Stochastic processes in biology
  - Gillespie algorithm for molecular simulations
  - Parameter estimation and optimization
  - Evolutionary trade-offs and Pareto optimality
  - Steady-state analysis and stability
  - Bistability and oscillatory dynamics

## 🛠 Technical Requirements

### Python Environment
All notebooks require the conda environment specified in `../environment.yml`:

```bash
# Create and activate environment
conda env create -f ../environment.yml
conda activate bio559r

# Start Jupyter
jupyter notebook
```

### Key Dependencies
- **Python**: NumPy, SciPy, Pandas, Matplotlib, Seaborn, Scikit-learn
- **R Integration**: rpy2 for R magic (`%%R` cells)
- **R Packages**: ggplot2, dplyr, tidyr, gridExtra
- **Specialized**: NetworkX, CVXPY (for optimization)

## 🎯 Learning Objectives

By completing these notebooks, students will:

1. **Model biological systems** using mathematical frameworks
2. **Implement numerical methods** for solving biological problems
3. **Analyze network structures** in biological systems
4. **Apply optimization techniques** to biological questions
5. **Understand stochastic effects** in molecular biology
6. **Visualize complex data** using Python and R
7. **Interpret mathematical results** in biological contexts

## 📊 Visualization Strategy

### Python + R Integration
- **Python**: Mathematical computations, data processing, simulations
- **R (via %%R magic)**: Advanced statistical plots using ggplot2
- **Benefits**: 
  - Python's computational power
  - R's superior statistical visualization
  - Seamless data transfer between languages

### Example Visualization Workflow
```python
# Python: Generate data
data = simulate_gene_expression()
df = pd.DataFrame(data)

# R: Create publication-quality plots
%%R -i df -w 12 -h 8
library(ggplot2)
ggplot(df, aes(x=time, y=expression)) + 
  geom_line() + theme_classic()
```

## 🧪 Practical Applications

### Real-World Examples
- **Gene regulatory networks**: Toggle switches, oscillators
- **Metabolic pathways**: Flux balance analysis
- **Population dynamics**: Predator-prey models
- **Drug responses**: Dose-response curves
- **Evolution**: Fitness landscapes and trade-offs

### Computational Skills
- Numerical integration of ODEs
- Matrix operations for network analysis
- Stochastic simulations
- Parameter estimation
- Data visualization

## 📝 Assessment and Exercises

### Interactive Problems
Each notebook includes:
- **Guided exercises**: Step-by-step implementations
- **Challenge problems**: Open-ended investigations
- **Your turn sections**: Independent problem-solving
- **Biological interpretation**: Connecting math to biology

### Example Exercise Types
- Implement the Gillespie algorithm
- Analyze bistability in gene circuits
- Perform flux balance analysis
- Identify network motifs
- Estimate kinetic parameters

## 🔗 Integration with Course Materials

### Connection to Other Modules
- **Module 1**: Introduction to systems thinking
- **Module 3**: Experimental techniques and data analysis
- **Module 4**: Computational tools (scanpy, scverse)
- **Module 5**: Applications in disease and drug discovery

### Preparation for Advanced Topics
- Single-cell RNA sequencing analysis
- Spatial transcriptomics
- Multi-omics integration
- Machine learning in biology

## 💡 Tips for Success

### Before Starting
1. **Review prerequisites**: Ensure comfort with calculus and basic programming
2. **Set up environment**: Test that all packages install correctly
3. **Allocate time**: Each notebook requires 3-4 hours of focused work

### During Exercises
1. **Run code sequentially**: Notebooks build on previous cells
2. **Experiment**: Modify parameters to see effects
3. **Visualize**: Use plots to understand mathematical concepts
4. **Connect to biology**: Always ask "what does this mean biologically?"

### Troubleshooting
- **R magic issues**: Ensure rpy2 is properly installed
- **Package conflicts**: Use the provided environment.yml
- **Slow simulations**: Reduce simulation parameters for testing
- **Plot rendering**: Restart kernel if visualizations don't appear

## 📚 Additional Resources

### Mathematical Background
- Strogatz, S. "Nonlinear Dynamics and Chaos"
- Edelstein-Keshet, L. "Mathematical Models in Biology"
- Alon, U. "An Introduction to Systems Biology"

### Computational Resources
- SciPy documentation: https://scipy.org/
- NetworkX tutorials: https://networkx.org/
- R ggplot2 reference: https://ggplot2.tidyverse.org/

### Biological Context
- Phillips, R. "Physical Biology of the Cell"
- Alberts, B. "Molecular Biology of the Cell"
- Systems biology review papers (provided in course materials)

---

**Happy modeling!** 🧬📊

*For questions or issues, contact the course instructor or teaching assistants.*

