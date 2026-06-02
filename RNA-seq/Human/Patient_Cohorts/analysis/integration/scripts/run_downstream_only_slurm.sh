#!/bin/bash
#SBATCH --job-name=liver_downstream
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/downstream_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/downstream_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00

# Quick re-run of downstream scripts (07-12) only.
# Use this when changing padj/LFC thresholds — no need to re-run
# the expensive steps 02-05 (voom, dream).
# Those produce full unfiltered results; thresholds are applied here.

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Downstream re-run started: $(date) ==="

run_step() {
  local num="$1"
  local name="$2"
  echo ""
  echo "=== Step $num: $name ==="
  Rscript "${SCRIPTS}/${name}" 2>&1
  if [ $? -ne 0 ]; then
    echo "FAILED: $name at $(date)"
    exit 1
  fi
  echo "  Done: $(date)"
}

run_step "07" "07_consensus_degs.R"
run_step "08" "08_pathway_analysis.R"
run_step "09" "09_batch_correction_umap.R"
run_step "10" "10_volcano_plots.R"
run_step "12" "12_library_intersection.R"

echo ""
echo "=== ALL DOWNSTREAM COMPLETE: $(date) ==="
