#!/bin/bash
#SBATCH --job-name=07f_canonical_audit
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/07f_canonical_audit_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/07f_canonical_audit_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting 07f_canonical_lr_audit on $(hostname)"
micromamba run -n rnaseq Rscript \
    Analysis/SingleCell/scripts/v2/07f_canonical_lr_audit.R
echo "[$(date)] DONE"
