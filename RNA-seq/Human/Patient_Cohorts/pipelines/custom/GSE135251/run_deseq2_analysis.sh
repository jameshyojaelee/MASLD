#!/bin/bash

# Lightweight wrapper to launch the standard DESeq2 batch job.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SBATCH_SCRIPT="${RESULTS_DIR}/run_deseq2.sh"

if [[ ! -f "${SBATCH_SCRIPT}" ]]; then
  echo "ERROR: Expected submission script not found at ${SBATCH_SCRIPT}" >&2
  exit 1
fi

echo "Submitting DESeq2 job via ${SBATCH_SCRIPT}"

if command -v sbatch >/dev/null 2>&1; then
  sbatch "${SBATCH_SCRIPT}"
else
  echo "sbatch not found; running script locally."
  bash "${SBATCH_SCRIPT}"
fi
