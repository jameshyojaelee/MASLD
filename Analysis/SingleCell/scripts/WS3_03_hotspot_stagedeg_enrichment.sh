#!/bin/bash
#SBATCH --job-name=hotspot
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/WS3_03_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/WS3_03_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "Host: $(hostname)  Job: ${SLURM_JOB_ID:-NA}  Start: $(date)"
micromamba run -n rnaseq Rscript \
  Analysis/SingleCell/scripts/WS3_03_hotspot_stagedeg_enrichment.R
echo "End: $(date)"
