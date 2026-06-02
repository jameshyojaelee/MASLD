#!/bin/bash
#SBATCH --job-name=figM_ccc_main
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/figM_ccc_main_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/figM_ccc_main_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting figM_stage_ccc_main on $(hostname)"
micromamba run -n rnaseq Rscript scripts/figures/figM_stage_ccc_main.R
echo "[$(date)] DONE"
