#!/bin/bash
#SBATCH --job-name=two_transition_substate
#SBATCH --output=logs/two_transition_substate_%j.log
#SBATCH --error=logs/two_transition_substate_%j.log
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

echo "=== Two-transition + B-bulk ==="
echo "Started: $(date)"

run_step() {
  local script=$1; local label=$2
  echo ""
  echo "============================================"
  echo "[$label] Running scripts/$script at $(date)"
  echo "============================================"
  Rscript scripts/$script
  local rc=$?
  if [ $rc -ne 0 ]; then echo "ERROR: $script exit $rc"; exit $rc; fi
}

run_step 250_two_transition_decomposition.R     "A two-transition decomposition"
run_step 251_f3_continuous_mixing_axis.R        "B-bulk F3a-F3b continuous axis"

echo ""
echo "=== Pillars A + B-bulk complete: $(date) ==="
