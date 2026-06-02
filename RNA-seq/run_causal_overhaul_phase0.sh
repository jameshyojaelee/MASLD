#!/bin/bash
# =============================================================================
# Phase 0: Infrastructure for MR/TWAS/COLOC Pipeline Overhaul
#
# Submits all Phase 0 jobs in parallel:
#   0B: R package installation
#   0C: Python environment setup
#   0D: GWAS liftover to hg19
#
# Note: 0A (SuSiE-COLOC LD fix) is a code change, not a SLURM job.
#
# Usage:
#   bash RNA-seq/run_causal_overhaul_phase0.sh
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

mkdir -p RNA-seq/logs

echo "=== Phase 0: Infrastructure Setup ==="
echo "Start: $(date)"

# 0B: Install R packages
echo ""
echo "--- Submitting 0B: R package installation ---"
JOB_0B=$(sbatch --parsable RNA-seq/install_mr_packages.sbatch)
echo "  Job ID: ${JOB_0B}"

# 0C: Python environment setup
echo ""
echo "--- Submitting 0C: Python environment ---"
JOB_0C=$(sbatch --parsable RNA-seq/setup_causal_inference_env.sbatch)
echo "  Job ID: ${JOB_0C}"

# 0D: GWAS liftover (needs rnaseq env for R + rtracklayer)
echo ""
echo "--- Submitting 0D: GWAS liftover ---"
JOB_0D=$(sbatch --parsable \
  --job-name=liftover_gwas \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=64G \
  --time=48:00:00 \
  --output=RNA-seq/logs/liftover_gwas_%j.out \
  --error=RNA-seq/logs/liftover_gwas_%j.err \
  --wrap="eval \"\$(micromamba shell hook -s bash)\" && micromamba activate rnaseq && Rscript RNA-seq/60_liftover_gwas_hg19.R")
echo "  Job ID: ${JOB_0D}"

echo ""
echo "=== All Phase 0 jobs submitted ==="
echo "Job IDs: 0B=${JOB_0B}, 0C=${JOB_0C}, 0D=${JOB_0D}"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB_0B},${JOB_0C},${JOB_0D} --format=JobID,JobName,State,Elapsed"
echo ""
echo "Phase 1/2/3 can start after these complete."
echo "Save these job IDs for dependency chains:"
echo "  export PHASE0_JOBS=\"${JOB_0B}:${JOB_0C}:${JOB_0D}\""
