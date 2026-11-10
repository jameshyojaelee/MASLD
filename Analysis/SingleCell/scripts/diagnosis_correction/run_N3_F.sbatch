#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/n3_F_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/n3_F_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
cp Analysis/SingleCell/scripts/diagnosis_correction/refit_N3_documented_F.R "$TASK/job_scripts/refit_N3_F_${SLURM_JOB_ID}.R"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$TASK/job_scripts/refit_N3_F_${SLURM_JOB_ID}.R"
