#!/bin/bash
#SBATCH --job-name=07b_chain_v3
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/07b_chain_v3_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/07b_chain_v3_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3

echo "[$(date)] starting 07b_chain_v3_hardened on $(hostname)"
micromamba run -n rnaseq Rscript \
    Analysis/SingleCell/scripts/v2/07b_chain_v3_hardened.R
echo "[$(date)] DONE"
