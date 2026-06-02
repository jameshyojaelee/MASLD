#!/bin/bash
#SBATCH --job-name=p1_discovery
#SBATCH --output=logs/p1_discovery_%j.log
#SBATCH --error=logs/p1_discovery_%j.log
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Phase 1 Discovery Pipeline ==="
echo "Host: $(hostname)"
echo "Started: $(date)"

run_step() {
  local script=$1
  local label=$2
  echo ""
  echo "======================================================================"
  echo "[$label] Running scripts/$script at $(date)"
  echo "======================================================================"
  Rscript scripts/$script
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "ERROR: $script exited with code $rc"
    exit $rc
  fi
}

run_step 241_f3_substate_discovery.R              "1.1 sub-state discovery"
run_step 242_continuous_fibrosis_index.R           "1.2 continuous index"
run_step 243_multi_anchor_signature_concordance.R  "1.3 signature concordance"
run_step 244_transition_effect_sizes.R             "1.4 effect sizes"
run_step 245_phase1_gates_check.R                  "1.7 gate evaluation"

echo ""
echo "=== All Phase 1 scripts complete: $(date) ==="
echo "Inspect outputs: RNA-seq/results/granular_staging/"
echo "Decision file: RNA-seq/results/granular_staging/phase1_gates.csv"
