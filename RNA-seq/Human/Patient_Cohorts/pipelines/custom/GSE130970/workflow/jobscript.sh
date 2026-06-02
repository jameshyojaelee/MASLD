#!/bin/bash
# properties = {properties}
set -euo pipefail
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
if [ -x "$MICROMAMBA" ]; then
  eval "$($MICROMAMBA shell hook --shell bash)"
  micromamba activate rnaseq >/dev/null 2>&1 || micromamba activate rnaseq
else
  echo "micromamba executable not found at $MICROMAMBA" >&2
  exit 1
fi
{exec_job}
