#!/bin/bash

# Show the tail of the most recent DESeq2 log.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

LOG_DIR="${RESULTS_DIR}/logs/deseq2"

if [[ ! -d "${LOG_DIR}" ]]; then
  echo "No DESeq2 log directory found at ${LOG_DIR}"
  exit 0
fi

LATEST_LOG="$(ls -1t "${LOG_DIR}" 2>/dev/null | head -n 1 || true)"

if [[ -z "${LATEST_LOG}" ]]; then
  echo "DESeq2 log directory exists but is empty."
  exit 0
fi

echo "Displaying last 40 lines of ${LOG_DIR}/${LATEST_LOG}"
tail -n 40 "${LOG_DIR}/${LATEST_LOG}"
