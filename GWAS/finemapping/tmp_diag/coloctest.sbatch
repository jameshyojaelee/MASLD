#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/ld_rebuild/coloctest_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/ld_rebuild/coloctest_%j.err
set -eo pipefail
set +u
export PATH="$HOME/.local/bin:$PATH"; eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
set -u
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/test_coloc_scale.R
