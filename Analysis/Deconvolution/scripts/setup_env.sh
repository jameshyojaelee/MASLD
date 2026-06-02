#!/usr/bin/env bash
set -euo pipefail

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}
ENV_FILE=${ENV_FILE:-RNA-seq/deconvolution/environment.yml}

if [[ ! -x "${MICROMAMBA}" ]]; then
  echo "micromamba not found at ${MICROMAMBA}" >&2
  exit 1
fi

"${MICROMAMBA}" create -y -p "${ENV_PREFIX}" -f "${ENV_FILE}"
"${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript RNA-seq/deconvolution/scripts/install_r_packages.R
