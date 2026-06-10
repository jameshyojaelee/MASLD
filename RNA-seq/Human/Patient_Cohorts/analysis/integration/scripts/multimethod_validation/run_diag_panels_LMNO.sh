#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --job-name=varpart
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/diagnostics/diag_LMNO_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/diagnostics/diag_LMNO_%j.err
# ===========================================================================
# Cross-method DE diagnostics (panels L/M/N/O) compute job.
#   Fresh DESeq2 + edgeR-QLF fits, raw varpart on 846, residual PVCA,
#   lambda/pi1/QQ tables. Reads canonical/dream/metafor p-vectors from disk.
# Writes results/integration/multimethod_validation/diagnostics/.
# Figures rendered separately by scripts/figures/figS_multimethod_diagnostics.R
# (light; submit afterok or run by hand once this completes).
# ===========================================================================
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/diagnostics
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
mkdir -p "${LOG_DIR}"

cd "${PROJECT}"
export MASLD_PROJECT_ROOT="${PROJECT}"

echo "=== DE diagnostics LMNO compute ==="
echo "Started: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-<interactive>}  CPUS: ${SLURM_CPUS_PER_TASK:-1}"
${RBIN} "${SCRIPTS}/diag_panels_LMNO_compute.R"
echo "Finished: $(date)"
