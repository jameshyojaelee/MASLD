#!/bin/bash
#SBATCH --job-name=05_donor_meta_v2
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/05_donor_meta_v2_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/05_donor_meta_v2_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs

echo "[$(date)] starting 05_donor_metadata_v2 on $(hostname)"
micromamba run -n rnaseq Rscript \
    Analysis/SingleCell/scripts/v2/05_donor_metadata_v2.R
echo "[$(date)] DONE"
