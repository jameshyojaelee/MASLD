#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/programs_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/programs_%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/PREPARED
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
mkdir -p "$TASK/job_scripts"
SCRIPT="$TASK/job_scripts/refit_programs_${SLURM_JOB_ID}.R"
cp Analysis/SingleCell/scripts/diagnosis_correction/refit_programs.R "$SCRIPT"
sha256sum "$SCRIPT" > "$SCRIPT.sha256"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla "$SCRIPT"
sha256sum --check "$SCRIPT.sha256"
