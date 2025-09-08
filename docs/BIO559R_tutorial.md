# BIO559R: Single-Cell and Spatial Transcriptomics Analysis Tutorial

**Author:** Manus AI (on behalf of Tommy W. Terooatea)

**Course:** BIO559R - Introduction to Systems Biology

## Introduction

Welcome to the BIO559R tutorial on single-cell and spatial transcriptomics analysis. This guide will walk you through setting up a complete computational environment for analyzing complex biological data using cutting-edge tools. You will learn how to:

*   Install Ubuntu on your computer.
*   Set up a conda environment with Python and R.
*   Use Jupyter notebooks with R magic for seamless integration of both languages.
*   Install and use the scverse ecosystem (including scanpy) for single-cell RNA-seq analysis.
*   Leverage R packages for advanced data visualization and statistical analysis.

This tutorial is designed for students with varying computational backgrounds. We will provide detailed instructions to ensure everyone can follow along.










## Chapter 1: Installing Ubuntu for Windows Users

For students who are not using a Mac, installing Ubuntu (a popular Linux distribution) is a great way to create a powerful environment for bioinformatics. You have two main options: dual-booting Ubuntu alongside Windows, or using a virtual machine. For this course, we recommend dual-booting for better performance, but both methods are acceptable.

### Dual-Booting Ubuntu with Windows

Dual-booting allows you to choose between Windows and Ubuntu when you start your computer. This gives you the full power of your machine for Ubuntu, which is ideal for computationally intensive tasks.

**1. Pre-installation Steps in Windows:**

*   **Back up your data:** Before making any changes to your system, it is crucial to back up all your important files to an external hard drive or cloud storage.
*   **Create a bootable Ubuntu USB drive:**
    *   Download the latest version of Ubuntu Desktop from the official website [1].
    *   Download a tool like Rufus [2] or Balena Etcher [3] to create a bootable USB drive from the Ubuntu ISO file you downloaded.
    *   You will need a USB drive with at least 8 GB of storage.
*   **Create space for Ubuntu:**
    *   Open the "Disk Management" tool in Windows.
    *   Right-click on your main drive (usually C:) and select "Shrink Volume."
    *   Shrink the volume to create at least 50 GB of unallocated space for Ubuntu. We recommend 100 GB or more if you have enough space.

**2. Installing Ubuntu:**

*   **Boot from the USB drive:**
    *   Restart your computer and enter the BIOS/UEFI settings. The key to enter the BIOS is usually F2, F10, F12, or DEL, depending on your computer manufacturer.
    *   In the BIOS settings, change the boot order to prioritize the USB drive.
    *   Save the changes and exit the BIOS. Your computer will now boot from the Ubuntu USB drive.
*   **Follow the installation wizard:**
    *   You will be greeted with the Ubuntu installer. Choose "Install Ubuntu."
    *   Select your language and keyboard layout.
    *   When you get to the "Installation type" screen, choose "Install Ubuntu alongside Windows Boot Manager." This is the simplest option for dual-booting.
    *   The installer will automatically detect the unallocated space you created earlier and partition it for Ubuntu.
    *   Follow the remaining prompts to set your location, create a user account, and start the installation.

**3. Post-installation:**

*   Once the installation is complete, restart your computer. You will now see a menu (called GRUB) that allows you to choose between booting into Ubuntu or Windows.

### Using a Virtual Machine

A virtual machine (VM) allows you to run Ubuntu within a window on your Windows desktop. This is a less invasive option, but it comes with a performance overhead.

**1. Install VirtualBox:**

*   Download and install VirtualBox [4] for Windows.

**2. Create a new virtual machine:**

*   Open VirtualBox and click "New."
*   Give your VM a name (e.g., "Ubuntu for BIO559R").
*   Set the type to "Linux" and the version to "Ubuntu (64-bit)."
*   Allocate at least 4 GB of RAM (8 GB or more is recommended).
*   Create a virtual hard disk. A dynamically allocated disk of at least 50 GB is a good starting point.

**3. Install Ubuntu on the VM:**

*   Start the newly created VM.
*   You will be prompted to select a start-up disk. Choose the Ubuntu ISO file you downloaded earlier.
*   The VM will boot from the ISO, and you can follow the same Ubuntu installation steps as in the dual-booting guide. However, in the "Installation type" screen, you can choose "Erase disk and install Ubuntu," as this will only affect the virtual hard disk you created.

---

### References

[1] Ubuntu Desktop: [https://ubuntu.com/download/desktop](https://ubuntu.com/download/desktop)
[2] Rufus: [https://rufus.ie/](https://rufus.ie/)
[3] Balena Etcher: [https://www.balena.io/etcher/](https://www.balena.io/etcher/)
[4] VirtualBox: [https://www.virtualbox.org/](https://www.virtualbox.org/)




## Chapter 2: Setting Up Your Conda Environment

Conda is an open-source package and environment management system that runs on Windows, macOS, and Linux. It is an essential tool for bioinformaticians as it allows you to create isolated environments containing specific versions of packages and their dependencies. This ensures that your projects are reproducible and do not conflict with each other.

For this course, we will use Miniconda, a minimal installer for conda. It is a smaller, bootstrap version of Anaconda that includes only conda, Python, the packages they both depend on, and a small number of other useful packages.

### 1. Installing Miniconda

*   **Download the Miniconda installer:**
    *   Go to the Miniconda documentation page [5] and download the installer for your operating system (Linux, in our case).
    *   Choose the Python 3.x version.
*   **Run the installer:**
    *   Open a terminal in Ubuntu.
    *   Navigate to the directory where you downloaded the installer.
    *   Run the installer script using the following command (replace the filename with the one you downloaded):

        ```bash
        bash Miniconda3-latest-Linux-x86_64.sh
        ```

    *   Follow the prompts in the installer. It is recommended to accept the default settings.
    *   When asked "Do you wish the installer to initialize Miniconda3 by running conda init?", we recommend you select "yes".
*   **Restart your terminal:**
    *   After the installation is complete, close and reopen your terminal for the changes to take effect.
    *   You should now see `(base)` at the beginning of your terminal prompt, which indicates that the base conda environment is active.

### 2. Creating a Conda Environment

Now that you have Miniconda installed, you can create a dedicated environment for our course. This will keep all the packages we need for single-cell analysis separate from other projects.

*   **Create a new environment:**
    *   In your terminal, run the following command to create a new environment named `bio559r` with Python 3.9:

        ```bash
        conda create -n bio559r python=3.9
        ```

    *   Conda will show you a list of packages that will be installed and ask you to confirm. Type `y` and press Enter.
*   **Activate the environment:**
    *   To start using the new environment, you need to activate it:

        ```bash
        conda activate bio559r
        ```

    *   Your terminal prompt should now start with `(bio559r)`, indicating that you are in the correct environment.

From now on, whenever you want to work on this course's projects, you should always activate the `bio559r` environment first.

---

### References

[5] Miniconda: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)




## Chapter 3: Installing Python and R Packages

Now that we have our `bio559r` conda environment, we can install the necessary packages for our analysis. We will install packages from different "channels," which are locations where packages are stored. The main channels we will use are `conda-forge` and `bioconda`.

### 1. Installing Python Packages for Single-Cell Analysis

The core of our Python-based analysis will be the scverse ecosystem, which is a collection of tools for single-cell omics data analysis. The main package we will use is `scanpy`.

*   **Activate your conda environment:**
    *   If you haven't already, activate the `bio559r` environment:

        ```bash
        conda activate bio559r
        ```

*   **Install scverse:**
    *   The easiest way to install the complete scverse is to install the `scverse` package, which includes `scanpy`, `anndata`, and other core packages. We will install it from the `conda-forge` channel:

        ```bash
        conda install -c conda-forge scverse
        ```

    *   This command will install all the necessary Python packages for single-cell analysis.

### 2. Installing R and R Packages

We will also use R for some statistical analysis and for creating publication-quality figures. We can install R and R packages directly within our conda environment.

*   **Install R and essential packages:**
    *   We will install R from the `conda-forge` channel. We will also install `r-essentials`, which includes a collection of commonly used R packages, and `r-irkernel` to enable the use of R in Jupyter notebooks.

        ```bash
        conda install -c conda-forge r-essentials r-irkernel
        ```

*   **Install Seurat and other R packages:**
    *   Seurat is a popular R package for single-cell analysis. We can install it from the `conda-forge` channel as well. We will also install `ggplot2` for plotting.

        ```bash
        conda install -c conda-forge r-seurat r-ggplot2
        ```

### 3. Verifying the Installation

After installing all the packages, you can verify that they are installed correctly.

*   **Check Python packages:**
    *   Open a Python interpreter by typing `python` in your terminal.
    *   Try importing `scanpy`:

        ```python
        import scanpy as sc
        print(sc.__version__)
        ```

    *   If the installation was successful, you should see the version number of scanpy printed to the console.
*   **Check R packages:**
    *   Open an R session by typing `R` in your terminal.
    *   Try loading the `Seurat` library:

        ```R
        library(Seurat)
        print(packageVersion("Seurat"))
        ```

    *   If the installation was successful, you should see the version number of Seurat printed to the console.




## Chapter 4: Configuring Jupyter Notebooks for Python and R

Jupyter notebooks provide an interactive, web-based environment for data analysis, visualization, and collaboration. We will configure Jupyter to work with both our Python and R kernels, allowing us to seamlessly switch between the two languages in a single notebook.

### 1. Installing Jupyter Notebook

First, we need to install Jupyter Notebook within our `bio559r` conda environment.

*   **Activate your conda environment:**

    ```bash
    conda activate bio559r
    ```

*   **Install Jupyter:**

    ```bash
    conda install -c conda-forge jupyter
    ```

### 2. Registering the R Kernel

We have already installed the `r-irkernel` package, but we need to register it with Jupyter so that it appears as an available kernel.

*   **Open an R session:**

    ```bash
    R
    ```

*   **Register the kernel:**

    ```R
    IRkernel::installspec()
    ```

    You might be asked if you want to install it for the current user or for all users. Installing for the current user is usually sufficient.

*   **Quit R:**

    ```R
    q()
    ```

### 3. Using R Magic in Jupyter Notebooks

One of the most powerful features of using R with Jupyter is "R magic," which allows you to run R code within a Python notebook. This is enabled by the `rpy2` package, which should have been installed as a dependency of `scverse`.

*   **Start a Jupyter Notebook:**

    ```bash
    jupyter notebook
    ```

    This will open a new tab in your web browser with the Jupyter interface.

*   **Create a new Python 3 notebook:**

    From the Jupyter dashboard, click on "New" and select "Python 3".

*   **Load the R magic extension:**

    In a new cell, type the following and run it:

    ```python
    %load_ext rpy2.ipython
    ```

*   **Run R code:**

    Now you can run R code in any cell by starting the cell with `%%R`.

    For example, you can create an R data frame and display it:

    ```R
    %%R
    my_df <- data.frame(x = 1:5, y = c("a", "b", "c", "d", "e"))
    print(my_df)
    ```

*   **Passing variables between Python and R:**

    You can also pass variables between your Python and R environments.

    *   **From Python to R:** Use the `-i` flag to pass a Python variable to the R environment.

        ```python
        import pandas as pd
        python_df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
        ```

        ```R
        %%R -i python_df
        print(head(python_df))
        ```

    *   **From R to Python:** Use the `-o` flag to pass an R variable back to the Python environment.

        ```R
        %%R -o r_df
        r_df <- data.frame(a = c(5, 6, 7), b = c(8, 9, 10))
        ```

        ```python
        print(r_df)
        ```

This seamless integration of Python and R within a single Jupyter notebook is a powerful paradigm for bioinformatics analysis, allowing you to leverage the strengths of both languages in your workflow.




## Chapter 5: Introduction to scverse and Scanpy

The scverse is a decentralized ecosystem of interoperable Python packages for single-cell omics data analysis. It provides a comprehensive suite of tools for all stages of the analysis workflow, from data loading and preprocessing to visualization and advanced modeling. At the core of the scverse is **Scanpy**, a powerful and scalable toolkit for single-cell gene expression analysis.

### The AnnData Data Structure

Before diving into Scanpy, it's essential to understand the `AnnData` data structure, which is the central data object in the scverse. An `AnnData` object stores the gene expression matrix, along with annotations for both cells (observations) and genes (variables). It is a highly efficient and flexible data structure that can be easily manipulated and extended.

An `AnnData` object has several key components:

*   `X`: The main data matrix, where rows are cells and columns are genes.
*   `obs`: A pandas DataFrame containing cell annotations (e.g., cell type, batch information).
*   `var`: A pandas DataFrame containing gene annotations (e.g., gene names, feature types).
*   `obsm`: A dictionary for storing multi-dimensional cell annotations (e.g., PCA, UMAP coordinates).
*   `varm`: A dictionary for storing multi-dimensional gene annotations.
*   `uns`: A dictionary for storing unstructured annotations (e.g., color palettes for plotting).

### A Typical Scanpy Workflow

A typical single-cell analysis workflow with Scanpy involves the following steps:

1.  **Data Loading:** Reading in the count matrix from various formats (e.g., 10x Genomics, Loom).
2.  **Preprocessing and Quality Control:**
    *   Filtering out low-quality cells and genes.
    *   Normalizing the data to account for differences in sequencing depth.
    *   Identifying highly variable genes.
    *   Scaling the data to have zero mean and unit variance.
3.  **Dimensionality Reduction:**
    *   Principal Component Analysis (PCA) to reduce the dimensionality of the data.
    *   t-SNE or UMAP for visualizing the data in 2D or 3D.
4.  **Clustering:**
    *   Graph-based clustering (e.g., Leiden algorithm) to identify cell clusters.
5.  **Finding Marker Genes:**
    *   Identifying genes that are differentially expressed between clusters.
6.  **Cell Type Annotation:**
    *   Annotating clusters with known cell types based on marker gene expression.

### Example: A Minimal Scanpy Analysis

Here is a minimal example of a Scanpy analysis workflow:

```python
import scanpy as sc

# Load the 10x Genomics dataset
adata = sc.read_10x_mtx(
    'path/to/your/data/',  # the directory with the .mtx, .tsv, and .barcodes files
    var_names='gene_symbols',                # use gene symbols for the variable names (genes)
    cache=True)                              # write a cache file for faster loading

# Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# Dimensionality reduction and clustering
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Visualization
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
```

This is just a brief introduction to the vast capabilities of Scanpy and the scverse. We encourage you to explore the official tutorials and documentation [6] to learn more.

---

### References

[6] Scanpy Tutorials: [https://scanpy-tutorials.readthedocs.io/](https://scanpy-tutorials.readthedocs.io/)




## Chapter 6: R for Visualization and Statistical Analysis

While Python and the scverse provide a powerful environment for single-cell analysis, R remains a popular choice for statistical modeling and creating publication-quality visualizations. Thanks to our conda environment and the `rpy2` package, we can seamlessly integrate R into our workflow.

### Seurat: A Powerful R Toolkit for Single-Cell Genomics

Seurat is one of the most widely used R packages for single-cell analysis. It provides a comprehensive suite of tools for quality control, analysis, and exploration of single-cell RNA-seq data. While we are using a Python-centric workflow in this tutorial, it is beneficial to be familiar with Seurat, as you may encounter it in other contexts or want to use some of its unique features.

### ggplot2: The Grammar of Graphics

`ggplot2` is a powerful and flexible R package for creating a wide range of static and interactive visualizations. It is based on the "Grammar of Graphics," a concept that allows you to build complex plots by combining simple components. `ggplot2` is an excellent choice for creating customized, publication-quality figures.

### Example: Using R for Visualization

Let's say we have performed our analysis in Scanpy and have our data stored in an `AnnData` object. We can easily transfer this data to R and use `ggplot2` to create a UMAP plot.

First, we need to extract the UMAP coordinates and cluster information from our `AnnData` object into a pandas DataFrame:

```python
# Assuming 'adata' is our AnnData object

umap_coords = adata.obsm['X_umap']
clusters = adata.obs['leiden']

plot_df = pd.DataFrame(data=umap_coords, columns=['UMAP1', 'UMAP2'])
plot_df['cluster'] = clusters.values
```

Now, we can pass this DataFrame to R and use `ggplot2` to create a plot:

```R
%%R -i plot_df

library(ggplot2)

ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5) +
    theme_classic() +
    labs(title = "UMAP of Cells by Cluster")
```

This is just a simple example, but it demonstrates the power of combining Python and R in your analysis. You can use the extensive statistical modeling capabilities of R to perform differential expression analysis, trajectory inference, or other advanced analyses, and then bring the results back into your Python environment.

We encourage you to explore the vast ecosystem of R packages for bioinformatics and statistical analysis. The Bioconductor project [7] is an excellent resource for finding R packages for genomics.

---

### References

[7] Bioconductor: [https://www.bioconductor.org/](https://www.bioconductor.org/)




## Chapter 7: Troubleshooting

This section provides solutions to common problems you might encounter while following this tutorial.

### Conda and Environment Issues

*   **`conda: command not found`**: This usually means that conda is not in your system's PATH. If you selected "yes" to initialize Miniconda during installation, restarting your terminal should fix this. If not, you may need to manually add the Miniconda `bin` directory to your PATH.

*   **Package installation conflicts**: Sometimes, conda may have trouble resolving dependencies, especially when installing packages from different channels. To minimize this, it is good practice to specify the channel for each package and to create a new environment for each project.

### Jupyter Notebook Issues

*   **R kernel not showing up**: If the R kernel is not available in Jupyter, make sure you have run `IRkernel::installspec()` in an R session within your `bio559r` conda environment.

*   **`rpy2` errors**: If you encounter errors related to `rpy2`, ensure that it is installed correctly in your conda environment. You can try reinstalling it with `conda install -c conda-forge rpy2`.

### General Tips

*   **Check your environment**: Always make sure you have the correct conda environment activated before installing packages or running your analysis.

*   **Read the error messages**: Error messages often provide valuable clues about what went wrong. Carefully read the message and try to understand the problem.

*   **Consult the documentation**: The documentation for the tools we are using (Scanpy, Seurat, conda, etc.) is an excellent resource for troubleshooting and learning more about their features.


