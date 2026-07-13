#!/bin/bash
#SBATCH --job-name=fwdval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/perturb_forward_validation/fwdval_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/perturb_forward_validation/fwdval_%j.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
echo "[$(date)] starting forward-validation TEST 1 on $(hostname)"
$PY Analysis/SingleCell/scripts/perturb_forward_validation.py
echo "[$(date)] python analysis done"
echo "[$(date)] ALL DONE"
