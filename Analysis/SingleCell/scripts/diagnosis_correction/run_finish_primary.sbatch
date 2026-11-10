#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/finish_primary_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/finish_primary_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TASK/job_scripts"
cp Analysis/SingleCell/scripts/diagnosis_correction/primary_diagnostics.R "$TASK/job_scripts/primary_diagnostics_${SLURM_JOB_ID}.R"
cp Analysis/SingleCell/scripts/diagnosis_correction/render_results.R "$TASK/job_scripts/render_results_${SLURM_JOB_ID}.R"
"$RS" "$TASK/job_scripts/primary_diagnostics_${SLURM_JOB_ID}.R"
"$RS" "$TASK/job_scripts/render_results_${SLURM_JOB_ID}.R"
touch "$TASK/FIGURES_COMPLETED"
