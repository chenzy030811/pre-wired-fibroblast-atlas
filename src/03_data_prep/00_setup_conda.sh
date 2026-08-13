#!/bin/bash
# One-time conda environment setup on Eddie
# Creates 'scvi' env with scvi-tools + scanpy + all deps
# Run on the login node (safe, no compute); ~15 min.

set -e

module load anaconda
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

ENV_NAME=scvi

if conda env list | grep -q "^${ENV_NAME}\s"; then
  echo ">>> Env '${ENV_NAME}' already exists. Activating..."
else
  echo ">>> Creating env '${ENV_NAME}' (this will take ~15 min)..."
  conda create -y -n ${ENV_NAME} python=3.10
fi

conda activate ${ENV_NAME}

echo ">>> Installing scvi-tools + friends..."
pip install --quiet \
  scvi-tools \
  scanpy \
  anndata \
  anndata2ri \
  pandas \
  numpy \
  scipy \
  matplotlib \
  seaborn \
  leidenalg \
  igraph \
  rpy2

echo ">>> DONE. To use:"
echo "  module load anaconda"
echo "  conda activate ${ENV_NAME}"
