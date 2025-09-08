#!/usr/bin/env python3
"""
BIO559R Installation Test Script

- Verifies core Python packages (scanpy stack) and rpy2
- Ensures R is discoverable (sets R_HOME from `R RHOME` if missing)
- Confirms Jupyter kernels for Python and R
- Runs a minimal Scanpy smoke test
"""

import os
import sys
import importlib
import subprocess

def ensure_r_home():
    """Ensure rpy2 will use THIS environment's R."""
    if os.environ.get("R_HOME"):
        return
    try:
        r_home = subprocess.check_output(["R", "RHOME"], text=True).strip()
        if r_home:
            os.environ["R_HOME"] = r_home
            print(f"  ✓ R_HOME -> {r_home}")
    except Exception as e:
        print(f"  ⚠️  Could not determine R_HOME automatically: {e}")

def test_python_packages():
    print("Testing Python packages...")
    required = [
        "numpy", "pandas", "scipy", "matplotlib", "seaborn",
        "sklearn", "scanpy", "anndata", "squidpy", "umap", "rpy2"
    ]
    failed = []
    for name in required:
        try:
            mod = "sklearn" if name == "sklearn" else name
            importlib.import_module(mod)
            print(f"  ✓ {name}")
        except ImportError as e:
            print(f"  ✗ {name} - {e}")
            failed.append(name)
    return failed

def test_r_packages():
    print("\nTesting R packages...")
    ensure_r_home()
    try:
        from rpy2.robjects.packages import importr
        required_r = ["base", "methods", "stats", "ggplot2", "dplyr", "Seurat"]
        failed_r = []
        for pkg in required_r:
            try:
                importr(pkg)
                print(f"  ✓ {pkg}")
            except Exception as e:
                print(f"  ✗ {pkg} - {e}")
                failed_r.append(pkg)
        return failed_r
    except ImportError:
        print("  ✗ rpy2 not available - cannot test R packages")
        return ["rpy2"]

def test_jupyter_kernels():
    print("\nTesting Jupyter kernels...")
    cmds = [["jupyter", "kernelspec", "list"], [sys.executable, "-m", "jupyter", "kernelspec", "list"]]
    for cmd in cmds:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                out = r.stdout
                print("Available kernels:")
                print(out)
                has_py = ("python3" in out) or ("Python (bio559r)" in out)
                has_r  = ("ir" in out) or ("R " in out)
                if has_py: print("  ✓ Python kernel available")
                else:      print("  ✗ Python kernel not found")
                if has_r:  print("  ✓ R kernel available")
                else:      print("  ✗ R kernel not found")
                return not (has_py and has_r)
        except Exception:
            continue
    print("  ✗ Failed to list Jupyter kernels")
    return True

def test_scanpy_basic():
    print("\nTesting scanpy basic functionality...")
    try:
        import scanpy as sc
        import numpy as np
        n_obs, n_vars = 100, 50
        X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
        adata = sc.AnnData(X)
        adata.var_names = [f"Gene_{i}" for i in range(n_vars)]
        adata.obs_names = [f"Cell_{i}" for i in range(n_obs)]
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
    print("=" * 50); print("BIO559R Installation Test"); print("=" * 50)
    failed_py = test_python_packages()
    failed_r  = test_r_packages()
    jup_fail  = test_jupyter_kernels()
    sc_fail   = test_scanpy_basic()

    print("\n" + "=" * 50); print("Test Summary"); print("=" * 50)
    total = len(failed_py) + len(failed_r) + (1 if jup_fail else 0) + (1 if sc_fail else 0)
    if total == 0:
        print("🎉 All tests passed! Your installation is ready for BIO559R.")
    else:
        print(f"⚠️  {total} issues found:")
        if failed_py: print("  - Failed Python packages:", ", ".join(failed_py))
        if failed_r:  print("  - Failed R packages:", ", ".join(failed_r))
        if jup_fail:  print("  - Jupyter kernel issues detected")
        if sc_fail:   print("  - Scanpy functionality test failed")
        print("\nIf R packages failed with 'methods.dylib', ensure R in this env is used:\n"
              "  conda activate bio559r && export R_HOME=\"$(R RHOME)\"")
    return total

if __name__ == "__main__":
    sys.exit(main())

