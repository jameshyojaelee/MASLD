#!/bin/bash
#SBATCH --job-name=plink2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/ld_rebuild/verify_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/ld_rebuild/verify_%j.err
set -eo pipefail
set +u
export PATH="$HOME/.local/bin:$PATH"; eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
set -u
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/verify_rebuild.R
