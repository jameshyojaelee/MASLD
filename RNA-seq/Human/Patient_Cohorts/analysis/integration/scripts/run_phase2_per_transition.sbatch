#!/bin/bash
#SBATCH --job-name=p2_transitions
#SBATCH --output=logs/p2_transitions_%j.log
#SBATCH --error=logs/p2_transitions_%j.log
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

# Override Phase 1 NO-GO gate; we want F1→F2 vs F2→F3 stratification to run
# regardless of F3 sub-state outcome.
export FORCE_PHASE2=1

echo "=== Phase 2 per-transition stratification ==="
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

run_step 246_substate_genetic_architecture.R    "246 COLOC × transition Fisher"
run_step 247_substate_drug_repurposing.R        "247 LINCS × drug class"
run_step 248_olink_substate_signature.R         "248 Olink panel signature"

echo ""
echo "=== Phase 2 per-transition complete: $(date) ==="
