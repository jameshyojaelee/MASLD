#!/bin/bash
#SBATCH --job-name=mprip-identity-refit
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/refit_%j.out
#SBATCH --error=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/refit_%j.err
set -euo pipefail
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export LD_PANEL=polyfun
unset UKBB_LD_DIR
work=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z
sha256sum -c "$work/refit_code.sha256"
python "$work/code/complete_identity_audit.py"
python "$work/code/prepare_mprip_refit.py"
Rscript "$work/code/refit_mprip.R"
sha256sum -c "$work/refit_code.sha256"
