#!/usr/bin/env bash
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail
cd "${PROJECT}"
Rscript scripts/figures/figS_tier1_contrast_summary.R
