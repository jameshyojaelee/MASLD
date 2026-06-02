#!/bin/bash
#SBATCH --job-name=04_fstage_v2
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/04_fstage_v2_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/04_fstage_v2_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2

echo "[$(date)] starting 04_fstage_pipeline_v2 on $(hostname)"
micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/v2/04_fstage_pipeline_v2.py
echo "[$(date)] DONE"
