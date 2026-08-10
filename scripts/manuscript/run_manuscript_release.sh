#!/usr/bin/env bash
set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
RELEASE_ID="${MANUSCRIPT_RELEASE_ID:-2026-07-15-r2}"
R_BIN="${RNASEQ_RSCRIPT:-/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript}"
PY_BIN="${SPATIAL_PYTHON:-/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python}"

export MASLD_PROJECT_ROOT="$BASE"
export MANUSCRIPT_RELEASE_ID="$RELEASE_ID"
export PYTHONNOUSERSITE=1

cd "$BASE"

"$R_BIN" scripts/manuscript/build_evidence_class_release.R
"$PY_BIN" scripts/manuscript/build_strict_drug_benchmark.py
"$PY_BIN" Analysis/Spatial/scripts/15l_gsmap_pcc_narrowing.py
"$PY_BIN" scripts/manuscript/build_gsmap_class_validation.py
"$R_BIN" scripts/manuscript/build_coding_architecture_release.R
"$R_BIN" scripts/figures/fig1_fig4_evidence_classes.R
"$PY_BIN" scripts/manuscript/check_release_consistency.py

printf 'Manuscript release %s complete.\n' "$RELEASE_ID"
