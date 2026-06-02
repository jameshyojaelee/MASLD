#!/bin/bash
# Submit an arm's model-runner array + dependent consensus job.
#
# Usage:
#   bash run_arm_driver.sh d1_mechanism            # submit D1 with default settings
#   MODALITY=zero_shot CONTEXT=all bash run_arm_driver.sh d3_synergy
#
# This wraps the P0-X1 fix: instead of running consensus inside the last array
# task (which loses the chain if any runner fails with set -e), we submit
# consensus as a separate job with --dependency=afterany:$ARRAY_JOB so it
# always runs, even if some runners failed.

set -euo pipefail

ARM_DIR="${1:?Usage: $0 <arm_dir>  (e.g., d1_mechanism, d2_reversal, ...)}"
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
ARM_PATH="$PROJECT_ROOT/Analysis/Perturbation/scripts/$ARM_DIR"

cd "$ARM_PATH"

# 1. Submit model-runner array
echo "[+] Submitting $ARM_DIR model-runner array ..."
ARRAY_JOB=$(sbatch --parsable run_arm.sbatch)
echo "    Array job: $ARRAY_JOB"

# 2. Submit consensus + backtester + atlas integrator as dependent job
#    (runs after ALL array tasks finish, regardless of success/failure)
echo "[+] Submitting $ARM_DIR consensus/backtest/atlas (dependency=afterany:$ARRAY_JOB) ..."
CONSENSUS_JOB=$(sbatch --parsable --dependency=afterany:$ARRAY_JOB run_arm_consensus.sbatch)
echo "    Consensus job: $CONSENSUS_JOB"

echo ""
echo "[done] Submitted $ARM_DIR pipeline:"
echo "  Array : $ARRAY_JOB (model runners)"
echo "  Then  : $CONSENSUS_JOB (consensus + backtester + atlas)"
echo ""
echo "Monitor:"
echo "  squeue -u jameslee --jobs=$ARRAY_JOB,$CONSENSUS_JOB"
