#!/bin/bash
#SBATCH --job-name=qualified_programs
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/qualified_programs_%j.out
#SBATCH --error=logs/qualified_programs_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export DIAGNOSIS_REFIT_ROOT="$PWD/figures/candidates/source-diagnosis-programs-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
export QUALIFIED_DONOR_METADATA="$PWD/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
mkdir -p "$DIAGNOSIS_REFIT_ROOT/source"
printf '%s\n' "$DIAGNOSIS_REFIT_ROOT" > "logs/qualified_programs_${SLURM_JOB_ID}.root"
cp scripts/analysis/bulk_refit_followup/refit_qualified_program_diagnosis.R scripts/analysis/bulk_refit_followup/refit_qualified_ambient.py scripts/analysis/bulk_refit_followup/run_qualified_program_diagnosis.sbatch "$DIAGNOSIS_REFIT_ROOT/source/"
"$bin/Rscript" --vanilla "$DIAGNOSIS_REFIT_ROOT/source/refit_qualified_program_diagnosis.R"
"$bin/python" "$DIAGNOSIS_REFIT_ROOT/source/refit_qualified_ambient.py"
