#!/bin/bash
#SBATCH --job-name=qualified_program_fig
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/qualified_program_fig_%j.out
#SBATCH --error=logs/qualified_program_fig_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export DIAGNOSIS_REFIT_ROOT="$(cat logs/qualified_programs_22237556.root)"
export FIGURE_CANDIDATE_ROOT="$DIAGNOSIS_REFIT_ROOT"
export PROGRAM_RELEASE_ROOT="$DIAGNOSIS_REFIT_ROOT/hotspot"
export SINGLECELL_STAGE_METADATA_FILE="$PWD/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected.tsv"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
export RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/bulk_refit_followup/render_verify_qualified_programs.py
