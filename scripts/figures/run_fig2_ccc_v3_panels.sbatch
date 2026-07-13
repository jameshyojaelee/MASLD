#!/bin/bash
#SBATCH --job-name=fig2_ccc_v3
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/fig2_ccc_v3_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/fig2_ccc_v3_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting ccc_v3_panels on $(hostname)"
micromamba run -n rnaseq Rscript scripts/figures/ccc_v3_panels.R
echo "[$(date)] DONE"
