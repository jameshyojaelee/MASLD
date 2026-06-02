#!/bin/bash
#SBATCH --job-name=343k_audit
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343k_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343k_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

micromamba run -n rnaseq Rscript \
  Analysis/SingleCell/scripts/343k_gated_fstage_audit.R
