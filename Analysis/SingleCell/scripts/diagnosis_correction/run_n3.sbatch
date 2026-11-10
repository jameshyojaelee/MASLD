#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/n3_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/n3_%j.err
set -euo pipefail
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f "$TASK/REFITS_PREPARED"
"$RS" "$TASK/refit_scripts/n3_model.R" "$TASK/N3/raw" results 1000 8
"$RS" "$TASK/refit_scripts/n3_model.R" "$TASK/N3/corrected" results_corrected 1000 8
"$RS" "$TASK/refit_scripts/n3_summary.R" "$TASK/N3/corrected" "$TASK/N3/raw"
touch "$TASK/N3_COMPLETED"
