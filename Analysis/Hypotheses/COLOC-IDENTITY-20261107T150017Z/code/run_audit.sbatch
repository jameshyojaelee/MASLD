#!/bin/bash
#SBATCH --job-name=coloc-identity-audit
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/audit_%j.out
#SBATCH --error=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/audit_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
task=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z
sha256sum -c "$task/SPEC.sha256"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$task/code/audit_metadata.py"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$task/code/audit_eqtl.R"
