#!/bin/bash
#SBATCH --job-name=block_b_diag
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/block_b_diag_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/block_b_diag_%j.err

# Block B diagnostics: bootstrap (B1) + permutation null (B2) + leverage (B3).
# Run sequentially in one job so they share the loaded LR data via filesystem
# cache and we stay under the qos=interactive 4-job ceiling.

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3

echo "[$(date)] BLOCK B start"
echo "[$(date)] B1 bootstrap (B=500)"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/v2/350_v3_bootstrap.R
echo "[$(date)] B2 permutation null (B=500)"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/v2/350b_permutation_null_v3.R
echo "[$(date)] B3 leverage diagnostics"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/v2/351_leverage_diagnostics_v3.R
echo "[$(date)] BLOCK B done"
