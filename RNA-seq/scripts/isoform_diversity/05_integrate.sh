#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/integrate_%x_%j.out
#SBATCH --error=logs/integrate_%x_%j.err
# sbatch --export=ALL,SPECIES=human 05_integrate.sh
set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"
micromamba run -n dtu Rscript 05_integrate.R "${SPECIES:-human}"
