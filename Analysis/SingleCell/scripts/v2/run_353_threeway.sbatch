#!/bin/bash
#SBATCH --job-name=353_threeway
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/353_threeway_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/353_threeway_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting 353_three_way_concordance"
micromamba run -n rnaseq Rscript \
    Analysis/SingleCell/scripts/v2/353_three_way_concordance.R
echo "[$(date)] DONE"
