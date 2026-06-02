#!/bin/bash
#SBATCH --job-name=s_v3_M08
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/M08_%j.out
#SBATCH --error=logs/M08_%j.err
# ---------------------------------------------------------------------------
# Module 08 — sex_v3 validation (v2 vs v3 + biology sanity + cross-pipeline).
# Wraps: scripts/sex_v3/08_validation.R
#
# Requires:
#   - sex_v3/sex_deg_classification_v3.csv      (Module 06)
#   - sex_v3/dream_M2_v3.rds                    (Module 04)
#   - sex_v3/intermediates/calibration_simulation.csv (Module 07, optional)
#   - results/integration/sex_deg_classification.csv  (v2)
#   - audit_sensitivity/sex_power_matched/per_gene_class_freq.csv
#   - audit_sensitivity/sex_composition_adjusted/composition_adjusted_sex_classification.csv
#
# Compute notes:
#   - Aggregation + figure rendering only. Light job.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== Module 08 (validation) start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 32G  Walltime: 8h"
echo

Rscript 08_validation.R

echo
echo "=== Module 08 end: $(date) ==="
