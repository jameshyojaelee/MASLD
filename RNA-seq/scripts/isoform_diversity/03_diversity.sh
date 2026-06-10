#!/bin/bash
#SBATCH --job-name=diversity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/diversity_%x_%j.out
#SBATCH --error=logs/diversity_%x_%j.err
# Usage: sbatch --export=ALL,SPECIES=human 03_diversity.sh
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
SPECIES=${SPECIES:-human}
echo "[$(date)] diversity $SPECIES"
micromamba run -n dtu Rscript 03_diversity.R "$SPECIES"
echo "[$(date)] DONE"
