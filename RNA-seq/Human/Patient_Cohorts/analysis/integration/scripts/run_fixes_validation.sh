#!/bin/bash
#SBATCH --job-name=liver_validate_fixes
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/validate_fixes_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/validate_fixes_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=4:00:00

# Targeted re-run to validate Fixes 1-6 on the existing 5 cohorts.
# Run list:
#   04: Variance Partition (Fix 4 & 6)
#   05: Dream Mega-Analysis (Fix 2 - Sex Covariate)
#   07: Consensus DEGs (Updates based on 05)
#   08: Pathway Analysis (Fix 5 - GSEA t-stat)
#   10: Volcano Plots (Updates based on 05/07)
#   11: Ortholog Mapping (Fix 1 - Unified Mouse)
#   12: Library Intersection (Updates based on 07/11)

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Fix Validation Run Started: $(date) ==="
echo "SLURM_JOB_ID: $SLURM_JOB_ID"

run_step() {
  local name="$1"
  echo ""
  echo "============================================================"
  echo "  Running: $name"
  echo "============================================================"
  Rscript "${SCRIPTS}/${name}" 2>&1
  if [ $? -ne 0 ]; then
    echo "FAILED: $name at $(date)"
    exit 1
  fi
  echo "  Completed: $(date)"
}

# 04: Variance Partition (Independent)
run_step "04_variance_partition.R" &

# 05: Dream (Key dependency for 07)
run_step "05_dream_mega_analysis.R"

# Wait for 04 (optional, but good governance)
wait

# 07: Consensus (Needs 05 and existing 06)
# Note: 06 (Meta-analysis) is not re-run as per-study DE (02) didn't change.
run_step "07_consensus_degs.R"

# Downstream Parallel Block
run_step "08_pathway_analysis.R" &
run_step "10_volcano_plots.R" &

wait

# 12: Library (Needs 11)
run_step "12_library_intersection.R"

echo ""
echo "============================================================"
echo "  VALIDATION RUN COMPLETE: $(date)"
echo "============================================================"
