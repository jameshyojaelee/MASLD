#!/bin/bash
#SBATCH --job-name=class_negcontrol
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/class_negcontrol_%j.out
#SBATCH --error=logs/class_negcontrol_%j.err
set -euo pipefail
export CLASS_COMPARISON_ROOT="$PWD/figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference"
export CLASS_NEGATIVE_CONTROL_ROOT="$CLASS_COMPARISON_ROOT/negative-control-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/bulk_refit_followup/class_pair_label_negative_control.py
cp scripts/analysis/bulk_refit_followup/class_pair_label_negative_control.py scripts/analysis/bulk_refit_followup/run_class_pair_label_negative_control.sbatch "$CLASS_NEGATIVE_CONTROL_ROOT/"
printf '%s\n' "$CLASS_NEGATIVE_CONTROL_ROOT" > "logs/class_negcontrol_${SLURM_JOB_ID}.root"
