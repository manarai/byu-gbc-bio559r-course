#!/usr/bin/env python3
"""
BIO559R Installation Test Script

This script tests whether all required packages are installed correctly
in the bio559r conda environment.
"""

import sys
import importlib
import subprocess

def test_python_packages():
    """Test Python package imports"""
    print("Testing Python packages...")
    
    required_packages = [
        'numpy',
        'pandas', 
        'scipy',
        'matplotlib',
        'seaborn',
        'sklearn',
        'scanpy',
        'anndata',
        'squidpy',
        'umap',
        'rpy2'
    ]
    
    failed_imports = []
    
    for package in required_packages:
        try:
            if package == 'sklearn':
                importlib.import_module('sklearn')
            elif package == 'umap':
                importlib.import_module('umap')
            else:
                importlib.import_module(package)
            print(f"  ✓ {package}")
        except ImportError as e:
            print(f"  ✗ {package} - {e}")
            failed_imports.append(package)
    
    return failed_imports

def test_r_packages():
    """Test R package availability through rpy2"""
    print("\nTesting R packages...")
    
    try:
        from rpy2.robjects.packages import importr
        from rpy2.robjects import r
        
        required_r_packages = [
            'base',
            'stats',
            'ggplot2',
            'dplyr',
            'Seurat'
        ]
        
        failed_r_imports = []
        
        for package in required_r_packages:
            try:
                importr(package)
                print(f"  ✓ {package}")
            except Exception as e:
                print(f"  ✗ {package} - {e}")
                failed_r_imports.append(package)
        
        return failed_r_imports
        
    except ImportError:
        print("  ✗ rpy2 not available - cannot test R packages")
        return ['rpy2']

def test_jupyter_kernels():
    """Test available Jupyter kernels"""
    print("\nTesting Jupyter kernels...")
    
    try:
        result = subprocess.run(['jupyter', 'kernelspec', 'list'], 
                              capture_output=True, text=True)
        
        if result.returncode == 0:
            output = result.stdout
            print("Available kernels:")
            print(output)
            
            # Check for Python and R kernels
            has_python = 'python3' in output
            has_r = 'ir' in output
            
            if has_python:
                print("  ✓ Python kernel available")
            else:
                print("  ✗ Python kernel not found")
                
            if has_r:
                print("  ✓ R kernel available")
            else:
                print("  ✗ R kernel not found")
                
            return not (has_python and has_r)
        else:
            print("  ✗ Failed to list Jupyter kernels")
            return True
            
    except FileNotFoundError:
        print("  ✗ Jupyter not found")
        return True

def test_scanpy_basic():
    """Test basic scanpy functionality"""
    print("\nTesting scanpy basic functionality...")
    
    try:
        import scanpy as sc
        import numpy as np
        
        # Create a simple test dataset
        n_obs, n_vars = 100, 50
        X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
        
        # Create AnnData object
        adata = sc.AnnData(X)
        adata.var_names = [f'Gene_{i}' for i in range(n_vars)]
        adata.obs_names = [f'Cell_{i}' for i in range(n_obs)]
        
        # Test basic preprocessing
        sc.pp.filter_cells(adata, min_genes=1)
        sc.pp.filter_genes(adata, min_cells=1)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        print("  ✓ Basic scanpy operations successful")
        return False
        
    except Exception as e:
        print(f"  ✗ Scanpy test failed: {e}")
        return True

def main():
    """Main test function"""
    print("=" * 50)
    print("BIO559R Installation Test")
    print("=" * 50)
    
    # Test Python packages
    failed_python = test_python_packages()
    
    # Test R packages
    failed_r = test_r_packages()
    
    # Test Jupyter kernels
    jupyter_failed = test_jupyter_kernels()
    
    # Test scanpy functionality
    scanpy_failed = test_scanpy_basic()
    
    # Summary
    print("\n" + "=" * 50)
    print("Test Summary")
    print("=" * 50)
    
    total_failures = len(failed_python) + len(failed_r) + jupyter_failed + scanpy_failed
    
    if total_failures == 0:
        print("🎉 All tests passed! Your installation is ready for BIO559R.")
    else:
        print(f"⚠️  {total_failures} issues found:")
        
        if failed_python:
            print(f"  - Failed Python packages: {', '.join(failed_python)}")
        
        if failed_r:
            print(f"  - Failed R packages: {', '.join(failed_r)}")
            
        if jupyter_failed:
            print("  - Jupyter kernel issues detected")
            
        if scanpy_failed:
            print("  - Scanpy functionality test failed")
        
        print("\nPlease check the installation instructions and try again.")
    
    return total_failures

if __name__ == "__main__":
    sys.exit(main())

