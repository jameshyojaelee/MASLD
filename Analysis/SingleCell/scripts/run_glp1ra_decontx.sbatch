#!/bin/bash
#SBATCH --job-name=decontX
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/decontX_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/decontX_%j.log

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCR=$ROOT/Analysis/SingleCell/scripts
mkdir -p "$SCR/logs"

echo "[sbatch] host=$(hostname) job=$SLURM_JOB_ID start=$(date)"

echo "[sbatch] STEP 1/2: export raw Hep+Endo counts (env: spatial)"
micromamba run -n spatial python "$SCR/glp1ra_decontx_export.py"

echo "[sbatch] STEP 2/2: decontX run + pool (env: decontx)"
micromamba run -n decontx Rscript "$SCR/glp1ra_decontx_run.R"

echo "[sbatch] DONE $(date)"
