#!/bin/bash
# Master submission script — chains all four pillars + aggregation.
# Run from RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/
#
# Usage:
#   bash run_robustness_master.sh
#
# Order:
#   1. Submit array jobs for Pillar A (CPSS, bootstrap, kfold), Pillar C (perm)
#      and stand-alone jobs for Pillar D (PVCA, SVA), Pillar B (transferability).
#   2. Capture all job IDs.
#   3. Submit aggregation job with afterok dependency on every upstream job.
#
# Environment: must be invoked from a node with sbatch in PATH; activates
# rnaseq inside each job rather than here.

set -euo pipefail

SUBMIT_LOG="logs/master_submit_$(date +%Y%m%d_%H%M%S).log"
mkdir -p logs

submit_array () {
  # Submit and capture the JobID
  local script=$1
  local jid
  jid=$(sbatch --parsable "$script")
  echo "$(date +%H:%M:%S) submit ${script} -> ${jid}" | tee -a "${SUBMIT_LOG}"
  echo "$jid"
}

CPSS_JID=$(submit_array  run_pillar_A_cpss.sbatch)
BOOT_JID=$(submit_array  run_pillar_A_bootstrap.sbatch)
KFOLD_JID=$(submit_array run_pillar_A_kfold.sbatch)
PERM_JID=$(submit_array  run_pillar_C_permutation.sbatch)
PVCA_JID=$(submit_array  run_pillar_D_residual_pvca.sbatch)
SVA_JID=$(submit_array   run_pillar_D_sva.sbatch)
PB_JID=$(submit_array    run_pillar_B_transferability.sbatch)

AGG_JID=$(sbatch --parsable \
  --dependency=afterok:${CPSS_JID}:${BOOT_JID}:${KFOLD_JID}:${PERM_JID}:${PVCA_JID}:${SVA_JID}:${PB_JID} \
  run_aggregate_robustness.sbatch)
echo "$(date +%H:%M:%S) submit aggregator -> ${AGG_JID}" | tee -a "${SUBMIT_LOG}"

cat <<EOF | tee -a "${SUBMIT_LOG}"

=== Robustness pipeline submitted ===
  CPSS pairs    : ${CPSS_JID}
  Bootstrap     : ${BOOT_JID}
  K-fold        : ${KFOLD_JID}
  Permutation   : ${PERM_JID}
  Residual PVCA : ${PVCA_JID}
  SVA           : ${SVA_JID}
  Pillar B      : ${PB_JID}
  Aggregator    : ${AGG_JID}  (afterok on all above)

Monitor with:
  squeue -u \$USER
  sacct -j ${AGG_JID}
EOF
