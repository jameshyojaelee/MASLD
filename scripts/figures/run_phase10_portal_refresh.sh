#!/bin/bash
# =============================================================================
# ⛔ RETIRED / DEAD — DO NOT RUN.
# This script points at a non-existent `masld-atlas-portal/` directory and
# invokes `python -m data.build_parquet`, a module that no longer exists.
# It is SUPERSEDED by the single web-data orchestrator:
#     scripts/portal/rebuild_web_data.sbatch
# (which rebuilds the atlas + all portal data under one 7-phase DAG). For a
# one-off atlas.parquet rebuild, run preprocess_atlas_data.py directly under
# the `spatial` env. Kept only for provenance; retired 2026-07-08.
# =============================================================================
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --job-name=phase10_portal
#SBATCH --output=scripts/figures/logs/phase10_portal_%j.out
#SBATCH --error=scripts/figures/logs/phase10_portal_%j.err

# Phase 10 Team C: rebuild atlas.parquet for masld-atlas-portal so the live
# Streamlit/Next.js portal serves the PolyFun-backed canonical SuSiE-COLOC
# values from the refreshed multi_evidence_atlas.csv.

set -eo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}/masld-atlas-portal"

# set -u disabled because conda/micromamba activate scripts reference unbound vars.
source ~/.bashrc
micromamba activate rnaseq

echo "============================================================"
echo "Phase 10 portal refresh"
echo "Started: $(date)"
echo "============================================================"

python -m data.build_parquet

# Mirror the parquet to the Next.js public/data dir so the deployed site
# serves the new atlas (build_parquet.py header notes this is the served copy).
NEXT_PUBLIC_DATA="${PROJECT_ROOT}/masld-atlas-v2/public/data"
if [[ -d "${NEXT_PUBLIC_DATA}" ]]; then
  echo "Mirroring atlas.parquet -> ${NEXT_PUBLIC_DATA}/atlas.parquet"
  cp -v atlas.parquet "${NEXT_PUBLIC_DATA}/atlas.parquet"
fi

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
