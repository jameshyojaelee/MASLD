#!/bin/bash
#SBATCH --job-name=N3source
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N3-source-diagnosis-20261014/logs/%x-%j.out
set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N3-source-diagnosis-20261014
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum "$TASK_ROOT"/scripts/*.R "$TASK_ROOT"/raw/donor_meta.tsv "$TASK_ROOT"/corrected/donor_meta.tsv > "$TASK_ROOT/input_script_sha256.txt"
"$RS" "$TASK_ROOT/scripts/03_donor_model.R" "$TASK_ROOT/raw" results 1000 8
"$RS" "$TASK_ROOT/scripts/03_donor_model.R" "$TASK_ROOT/corrected" results_corrected 1000 8
"$RS" "$TASK_ROOT/scripts/04_ambient_summary.R" "$TASK_ROOT/corrected" "$TASK_ROOT/raw"
