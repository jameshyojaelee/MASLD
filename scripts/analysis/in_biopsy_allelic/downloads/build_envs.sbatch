#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
# Build the two user-approved environments (2026-09-23) off the login node.
#   allelic_qtl : Models A and B (tensorQTL, decoupler, NumPyro, pysam)
#   he_molecular: Model C (OpenSlide with DICOM support, timm, IDC download client)
# Usage: sbatch --job-name=<n> --output=<log>/%j.out build_envs.sbatch <log_dir>
set -euo pipefail
log=$1
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
export PYTHONNOUSERSITE=1

have_env() { "$mm" env list | awk -v e="$1" '$1 == e {f = 1} END {exit !f}'; }
if ! have_env allelic_qtl; then
  "$mm" create -y -n allelic_qtl -c conda-forge -c bioconda \
    python=3.11 numpy scipy pandas pyarrow statsmodels scikit-learn pysam pgenlib \
    bcftools=1.21 htslib=1.21 openjdk=21
  "$mm" run -n allelic_qtl python -m pip install --no-cache-dir \
    --index-url https://download.pytorch.org/whl/cpu torch
  "$mm" run -n allelic_qtl python -m pip install --no-cache-dir tensorqtl decoupler "jax[cpu]" numpyro arviz
fi
# tensorqtl imports pandas_plink unconditionally (genotypeio.py); pip did not pull it in.
"$mm" install -y -n allelic_qtl -c conda-forge -c bioconda pandas-plink
"$mm" run -n allelic_qtl python -c "import tensorqtl, decoupler, numpyro, pysam, pgenlib; print('allelic_qtl ok', tensorqtl.__version__, numpyro.__version__)"
"$mm" run -n allelic_qtl python -m pip freeze > "$log/allelic_qtl.pip_freeze.txt"

have_env he_molecular || "$mm" create -y -n he_molecular -c conda-forge python=3.11 numpy scipy pandas scikit-learn scikit-image pyarrow
"$mm" run -n he_molecular python -m pip install --no-cache-dir \
  --index-url https://download.pytorch.org/whl/cu124 torch torchvision
"$mm" run -n he_molecular python -m pip install --no-cache-dir \
  openslide-python openslide-bin timm huggingface_hub idc-index wsidicom opencv-python-headless
"$mm" run -n he_molecular python -c "import openslide, timm, idc_index; print('he_molecular ok', openslide.__library_version__, timm.__version__)"
"$mm" run -n he_molecular python -m pip freeze > "$log/he_molecular.pip_freeze.txt"
echo envs done
