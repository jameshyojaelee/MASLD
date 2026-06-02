#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=72:00:00

# SuSiE-COLOC wrapper — run via run_susie_coloc.sh which sets GWAS_NAME + log paths
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

module load PLINK/2.0a5.13

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd "${BASE}"
Rscript RNA-seq/35s_susie_coloc_broadaway.R
