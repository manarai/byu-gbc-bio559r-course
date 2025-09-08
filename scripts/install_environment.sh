#!/usr/bin/env bash
# BIO559R Environment Installation Script (robust, Apple Silicon–friendly)

set -euo pipefail

banner() { printf "\n==========================================\n%s\n==========================================\n" "$1"; }
ok()    { printf "✓ %s\n" "$1"; }
err()   { printf "✗ %s\n" "$1" >&2; }
info()  { printf "→ %s\n" "$1"; }

banner "BIO559R Environment Installation Script"

# --- 0) Locate conda and initialize shell integration
if ! command -v conda >/dev/null 2>&1; then
  err "conda is not installed or not in PATH. Install Miniconda/Mambaforge first."
  exit 1
fi
ok "conda found"

# Ensure 'conda activate' works inside scripts
# (do this BEFORE any 'conda config' or 'conda install' calls)
eval "$(conda shell.bash hook)"

# --- 1) Verify environment.yml presence
if [[ ! -f "environment.yml" ]]; then
  err "environment.yml not found in $(pwd). Run this script from the repo root."
  exit 1
fi
ok "environment.yml found"

# --- 2) Prefer libmamba solver + mamba for reliability/speed
info "Ensuring libmamba solver + mamba are available..."
conda install -n base -c conda-forge -y conda-libmamba-solver mamba >/dev/null
conda config --set solver libmamba >/dev/null
ok "libmamba solver set; mamba available"

# --- 3) Set channels + strict priority for stable solves
info "Configuring channels and priority..."
conda config --add channels conda-forge >/dev/null 2>&1 || true
conda config --add channels bioconda    >/dev/null 2>&1 || true
conda config --add channels defaults    >/dev/null 2>&1 || true
conda config --set channel_priority strict >/dev/null
ok "Channels set (conda-forge, bioconda, defaults) with strict priority"

# --- 4) Show platform (helps debug Apple Silicon vs Rosetta)
PLAT="$(conda info | awk -F': ' '/platform/{print $2}')"
info "Conda platform: ${PLAT}"

# --- 5) Create env with mamba (retry once with extra hints on Apple Silicon)
ENV_NAME="bio559r"
banner "Creating conda environment '${ENV_NAME}' (this may take a few minutes)..."

set +e
mamba env create -f environment.yml
CREATE_RC=$?
set -e

if [[ $CREATE_RC -ne 0 ]]; then
  err "Environment creation failed on first attempt."

  if [[ "${PLAT}" == "osx-arm64" ]]; then
    cat <<'EOS'

Tips for Apple Silicon (osx-arm64):
  • Ensure packages in environment.yml exist for osx-arm64 (some bioconda recipes are osx-64 only).
  • Prefer conda-forge builds; avoid hard pins to osx-64-only packages.
  • As a last resort for osx-64-only packages, you can build an x86_64 env via Rosetta:
        # new terminal under Rosetta OR prefix commands with: arch -x86_64
        CONDA_SUBDIR=osx-64 conda create -n bio559r_x86 python=3.11
        conda activate bio559r_x86
        conda config --env --set subdir osx-64
        conda install -c conda-forge -c bioconda <needed-packages>

EOS
  fi
  exit 1
fi

ok "Environment '${ENV_NAME}' created successfully"

# --- 6) Activate env
info "Activating environment '${ENV_NAME}'..."
conda activate "${ENV_NAME}"

# --- 7) Register kernels (Python + R) in user scope
banner "Registering Jupyter kernels"

# Python kernel (explicit name/display-name for clarity)
if command -v python >/dev/null 2>&1; then
  python -m ipykernel install --user --name "${ENV_NAME}" --display-name "Python (${ENV_NAME})" >/dev/null
  ok "Python kernel registered"
else
  err "Python not found inside the environment (unexpected)."
fi

# R kernel via IRkernel (if R is present)
if command -v R >/dev/null 2>&1; then
  info "Registering R kernel (IRkernel) for current user..."
  R --slave -e 'if (!requireNamespace("IRkernel", quietly=TRUE)) { message("IRkernel not installed; install via conda: r-irkernel"); quit(status=1) } else { IRkernel::installspec(user = TRUE); message("IRkernel registered (user scope).") }'
  if [[ $? -eq 0 ]]; then
    ok "R kernel registered"
  else
    err "R kernel registration skipped/failed. Ensure r-irkernel is in environment.yml."
  fi
else
  info "R not found in environment; skipping R kernel registration."
fi

banner "Installation completed successfully!"
printf "\nTo activate the environment, run:\n  conda activate %s\n" "${ENV_NAME}"
printf "\nTo start Jupyter Lab (recommended), run:\n  jupyter lab\n"
printf "\nYou should now see **Python (%s)** and **R** kernels in Jupyter.\n\n" "${ENV_NAME}"
