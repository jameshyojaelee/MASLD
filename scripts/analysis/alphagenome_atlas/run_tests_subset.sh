#!/usr/bin/env bash
set -euo pipefail
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$MASLD_PROJECT_ROOT/scripts/analysis/alphagenome_atlas"
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python -m unittest "$@" 2>&1 | tail -25
