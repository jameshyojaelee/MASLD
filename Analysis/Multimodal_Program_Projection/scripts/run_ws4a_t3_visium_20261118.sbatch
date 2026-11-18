#!/usr/bin/env bash
#SBATCH --job-name=ws4a_t3_visium
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/ws4a_t3_visium_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/ws4a_t3_visium_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"

"$PY_BIN" Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection_20261118.py
"$PY_BIN" Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/31_run_full_composition_sensitivity_20261118.py
"$PY_BIN" Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN4/20261118/an4_spatial_20261118.py
"$PY_BIN" scripts/analysis/hypotheses_2026q4/V4_spatial_replication/v4_spatial_replication_20261118.py
