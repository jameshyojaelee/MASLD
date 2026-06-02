#!/bin/bash
#SBATCH --job-name=07c_endocrine
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/07c_endocrine_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/07c_endocrine_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting 07c_endocrine_lr_filter on $(hostname)"
micromamba run -n rnaseq Rscript \
    Analysis/SingleCell/scripts/v2/07c_endocrine_lr_filter.R
echo "[$(date)] DONE"
