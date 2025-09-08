#!/usr/bin/env bash
# BIO559R Environment Installation Script (robust + Apple Silicon–friendly)

set -euo pipefail

banner() { printf "\n==========================================\n%s\n==========================================\n" "$1"; }
ok()    { printf "✓ %s\n" "$1"; }
err()   { printf "✗ %s\n" "$1" >&2; }
info()  { printf "→ %s\n" "$1"; }

banner "BIO559R Environment Installation Script"

# 0) Ensure conda exists and initialize shell integration for THIS script
if ! command -v conda >/dev/null 2>&1; then
  err "conda is not installed or not on PATH. Install Miniconda first."
  exit 1
fi
# Make 'conda activate' work in scripts
eval "$(conda shell.bash hook)"
ok "conda found"

# 1) Verify env file
if [[ ! -f "environment.yml" ]]; then
  err "environment.yml not found in $(pwd). Run this script from the repo root."
  exit 1
fi
ok "environment.yml found"

# 2) Use libmamba solver + strict channel priority; prefer conda-forge/bioconda only
info "Configuring solver and channels..."
conda install -n base -c conda-forge -y conda-libmamba-solver mamba >/dev/null
conda config --set solver libmamba
# keep only forge/bioconda to avoid defaults/repo.anaconda.com conflicts
conda config --remove-key channels >/dev/null 2>&1 || true
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
ok "libmamba set; channels = [conda-forge, bioconda] (strict)"

# 3) Create env
ENV_NAME="bio559r"
banner "Creating conda environment '${ENV_NAME}' (this may take a few minutes)..."
set +e
mamba env create -f environment.yml
CREATE_RC=$?
set -e
if [[ $CREATE_RC -ne 0 ]]; then
  err "Environment creation failed."
  echo "Tip: On Apple Silicon, ensure packages exist for osx-arm64. Avoid osx-64-only bioconda recipes."
  exit 1
fi
ok "Environment '${ENV_NAME}' created"

# 4) Activate env
info "Activating environment '${ENV_NAME}'..."

# guard flaky activate hooks that read unset vars  <<< YOUR REQUESTED LINES
export GFORTRAN="${GFORTRAN:-}"
export FC_FOR_BUILD="${FC_FOR_BUILD:-}"
export CC_FOR_BUILD="${CC_FOR_BUILD:-}"
export CXX_FOR_BUILD="${CXX_FOR_BUILD:-}"

conda activate "${ENV_NAME}"

# 5) Point rpy2 to THIS env's R (warning-proof)
if command -v R >/dev/null 2>&1; then
  # Use only the last non-empty line from `R RHOME` (avoids 'WARNING: ignoring environment value of R_HOME')
  export R_HOME="$(R RHOME | tail -n 1)"
  if [[ -f "$R_HOME/lib/libR.dylib" || -f "$R_HOME/lib/libR.so" ]]; then
    ok "R_HOME set to: ${R_HOME}"
  else
    info "Could not validate libR under R_HOME (${R_HOME}); continuing but rpy2 may fail if R is misconfigured."
  fi
else
  info "R not found in environment; skipping R_HOME setup."
fi

# 6) Register Jupyter kernels (user scope; avoids admin perms)
banner "Registering Jupyter kernels"
if command -v python >/dev/null 2>&1; then
  python -m ipykernel install --user --name "${ENV_NAME}" --display-name "Python (${ENV_NAME})" >/dev/null
  ok "Python kernel registered"
else
  err "Python not found in the environment (unexpected)."
fi

if command -v R >/dev/null 2>&1; then
  # Install IRkernel spec if package is present; otherwise instruct the user
  R --slave -e 'if (!requireNamespace("IRkernel", quietly=TRUE)) { message("IRkernel not installed; add `r-irkernel` to environment.yml."); quit(status=1) } else { IRkernel::installspec(user = TRUE); message("R kernel registered (user scope).") }'
  if [[ $? -eq 0 ]]; then ok "R kernel registered"; else err "R kernel registration failed"; fi
else
  info "R not found; skipping R kernel registration."
fi

banner "Installation completed!"
echo "To activate:  conda activate ${ENV_NAME}"
echo "To start:     jupyter lab"
echo "You should see kernels: Python (${ENV_NAME}) and R"

