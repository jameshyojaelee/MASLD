#!/bin/bash

# Summarize all available DESeq2 logs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

LOG_DIR="${RESULTS_DIR}/logs/deseq2"

if [[ ! -d "${LOG_DIR}" ]]; then
  echo "No DESeq2 log directory found at ${LOG_DIR}"
  exit 0
fi

echo "Available DESeq2 logs under ${LOG_DIR}:"
ls -lh "${LOG_DIR}"
