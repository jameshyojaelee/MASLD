#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/verify_models_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/verify_models_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
cp Analysis/SingleCell/scripts/diagnosis_correction/verify_program_nulls.R "$TASK/job_scripts/verify_program_nulls_${SLURM_JOB_ID}.R"
"$RS" "$TASK/job_scripts/verify_program_nulls_${SLURM_JOB_ID}.R"
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/verify_results.py --phase models
touch "$TASK/MODEL_CHECKS_COMPLETED"
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/collect_corrected_results.py
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/prepare_manuscript.py
