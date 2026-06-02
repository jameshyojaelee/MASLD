#!/bin/bash

# Provide a concise view of the results directory layout.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

if [[ ! -d "${RESULTS_DIR}" ]]; then
  echo "Results directory not found: ${RESULTS_DIR}" >&2
  exit 1
fi

echo "Results structure for $(basename "${RESULTS_DIR}")"
echo "================================================"
find "${RESULTS_DIR}" -maxdepth 2 -mindepth 1 -type d | sort
