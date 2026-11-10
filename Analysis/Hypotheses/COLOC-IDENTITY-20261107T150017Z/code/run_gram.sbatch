#!/bin/bash
#SBATCH --job-name=mprip-genotype-LD
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/gram_%j.out
#SBATCH --error=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/gram_%j.err
set -euo pipefail
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
work=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z
sha256sum -c "$work/gram_code.sha256"
Rscript "$work/code/map_gram_variants.R"
python "$work/code/build_genotype_gram.py"
Rscript "$work/code/refit_gram.R"
sha256sum -c "$work/gram_code.sha256"
