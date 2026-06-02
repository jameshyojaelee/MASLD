#!/bin/bash
#SBATCH --job-name=gtex_coloc_ukbb
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00

# GTEx v8 Liver × UKBB ALT replication COLOC
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd "${BASE}"
Rscript RNA-seq/35x_gtex_coloc_ukbb.R
