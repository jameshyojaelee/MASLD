#!/bin/bash
#SBATCH --job-name=352a_scan
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/352a_scanorama_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/352a_scanorama_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting 352a_scanorama_integration on $(hostname)"
micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/v2/352a_scanorama_integration.py
echo "[$(date)] DONE"
