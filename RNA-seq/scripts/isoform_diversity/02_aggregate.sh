#!/bin/bash
#SBATCH --job-name=tximport
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=logs/aggregate_%x_%j.out
#SBATCH --error=logs/aggregate_%x_%j.err
# Usage: sbatch --export=ALL,SPECIES=human 02_aggregate.sh
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
SPECIES=${SPECIES:-human}
echo "[$(date)] aggregate $SPECIES"
micromamba run -n dtu Rscript 02_aggregate.R "$SPECIES"
echo "[$(date)] DONE"
