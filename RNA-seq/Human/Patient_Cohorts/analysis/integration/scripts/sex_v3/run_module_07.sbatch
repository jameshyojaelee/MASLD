#!/bin/bash
#SBATCH --job-name=s_v3_M07
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M07_%j.out
#SBATCH --error=logs/M07_%j.err
# ---------------------------------------------------------------------------
# Module 07 — sex_v3 calibration (synthetic-truth + LOCO posterior predictive).
# Wraps: scripts/sex_v3/07_calibration.R
#
# Requires:
#   - intermediates/sex_v3_input.rds (Module 01)
#   - intermediates/sva_factors.rds  (Module 02)
#   - intermediates/age_mi.rds       (Module 03)
#   - dream_M2_v3.rds                (Module 04)
#   - intermediates/mashr_fit.rds    (Module 05)
#   - mashr_posterior_summary.csv    (Module 05)
#
# Compute notes:
#   - LOCO is the heavy step (5 dream M2 fits sequentially, each ~3-6h).
#   - Synthetic simulation requires ONE additional dream fit on a ~2,000-gene
#     panel — fast (~30-60 min) because of the small N_genes.
#   - 96 G memory budget matches Module 04.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

# Saturate BLAS threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "=== Module 07 (calibration) start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 96G  Walltime: 48h"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo

Rscript 07_calibration.R

echo
echo "=== Module 07 end: $(date) ==="
