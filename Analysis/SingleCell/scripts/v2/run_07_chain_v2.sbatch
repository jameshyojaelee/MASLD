#!/bin/bash
#SBATCH --job-name=07_chain_v2
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --time=12:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/07_chain_v2_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/07_chain_v2_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2

echo "[$(date)] starting 07_chain_346_349_v2 on $(hostname)"
micromamba run -n rnaseq Rscript \
    Analysis/SingleCell/scripts/v2/07_chain_346_349_v2.R
echo "[$(date)] DONE"
