#!/bin/bash
#SBATCH --job-name=s_v3_driver
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/driver_%j.out
#SBATCH --error=logs/driver_%j.err
# ---------------------------------------------------------------------------
# sex_v3 driver — chains Modules 01 -> 02 -> 03 -> 04 via afterok dependency.
#
# Each child submission keeps its own SLURM resource block (cpu/16/96G/48h);
# this driver itself is a lightweight orchestrator that just dispatches the
# chain and prints job IDs. Run with `sbatch driver_v3.sbatch` from the
# sex_v3 directory.
#
# If a module step needs a re-run, manually submit run_module_0X.sbatch with
# the appropriate --dependency=afterok:<prior_id>.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

echo "=== sex_v3 driver start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "Submitting Modules 01 -> 02 -> 03 -> 04 (afterok chain)"
echo

# Submit Module 01
J01=$(sbatch --parsable run_module_01.sbatch)
echo "Submitted M01: ${J01}"

# Submit Module 02 (depends on M01)
J02=$(sbatch --parsable --dependency=afterok:${J01} run_module_02.sbatch)
echo "Submitted M02: ${J02}  (depends on M01=${J01})"

# Submit Module 03 (depends on M01 — independent of M02; could run in parallel)
J03=$(sbatch --parsable --dependency=afterok:${J01} run_module_03.sbatch)
echo "Submitted M03: ${J03}  (depends on M01=${J01})"

# Submit Module 04 (depends on M02 AND M03)
J04=$(sbatch --parsable --dependency=afterok:${J02}:${J03} run_module_04.sbatch)
echo "Submitted M04: ${J04}  (depends on M02=${J02}, M03=${J03})"

echo
echo "Chain summary:"
echo "  M01: ${J01}"
echo "  M02: ${J02}  (after M01)"
echo "  M03: ${J03}  (after M01, parallel to M02)"
echo "  M04: ${J04}  (after M02 & M03)"
echo
echo "Verify within 1-2 min:"
echo "  squeue -j ${J01},${J02},${J03},${J04}"
echo
echo "=== sex_v3 driver end: $(date) ==="
