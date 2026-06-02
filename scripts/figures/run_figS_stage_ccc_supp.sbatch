#!/bin/bash
#SBATCH --job-name=figS_ccc_supp
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/figS_ccc_supp_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/figS_ccc_supp_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting figS_stage_ccc_supp on $(hostname)"
micromamba run -n rnaseq Rscript scripts/figures/figS_stage_ccc_supp.R
echo "[$(date)] DONE"
