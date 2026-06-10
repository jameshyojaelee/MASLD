#!/bin/bash
#SBATCH --job-name=saturn
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --output=logs/dtu_%x_%j.out
#SBATCH --error=logs/dtu_%x_%j.err
# Usage: sbatch --export=ALL,SPECIES=human,CONTRAST=group 04_dtu.sh
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
SPECIES=${SPECIES:-human}; CONTRAST=${CONTRAST:-group}
echo "[$(date)] DTU $SPECIES $CONTRAST"
micromamba run -n dtu Rscript 04_dtu.R "$SPECIES" "$CONTRAST"
echo "[$(date)] DONE"
