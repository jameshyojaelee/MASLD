#!/usr/bin/env bash
#SBATCH --job-name=music_install
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=RNA-seq/deconvolution/logs/%x_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/%x_%j.err

set -euo pipefail

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}
export R_LIBS_USER=
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

if [[ ! -x "${MICROMAMBA}" ]]; then
  echo "micromamba not found at ${MICROMAMBA}" >&2
  exit 1
fi

# Ensure env exists (setup_env.sh will create if missing)
if [[ ! -d "${ENV_PREFIX}" ]]; then
  RNA-seq/deconvolution/scripts/setup_env.sh
fi

# Install MuSiC if missing
set +e
"${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript -e 'if (!requireNamespace("MuSiC", quietly=TRUE)) { quit(status=2) }'
status=$?
set -e
if [[ $status -ne 0 ]]; then
  echo "MuSiC not found; installing..." >&2
  "${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript RNA-seq/deconvolution/scripts/install_r_packages.R
else
  echo "MuSiC already installed." >&2
fi

"${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript -e 'if (!requireNamespace("MuSiC", quietly=TRUE)) { quit(status=3) }'
