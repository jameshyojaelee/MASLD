#!/bin/bash
#SBATCH --job-name=scploidy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-18%6
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ploidy/logs/scploidy_%A_%a.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ploidy/logs/scploidy_%A_%a.err

# Run scPloidy on hepatocyte fragments per donor (array task = donor index).
# Submit with: sbatch Analysis/SingleCell/scripts/run_340_scploidy.sbatch
# After completion, summary is at Analysis/SingleCell/results_gpu_v2/ploidy/scploidy_donor_summary.csv

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

# Source env from user profile so micromamba is on PATH
source ~/.bashrc 2>/dev/null || true

# Donor list (must match Analysis/ATAC/Human_Multiome/results/fragments/)
DONORS=(D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 D14 D15 D16 D17 D18)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
export DONOR_ID="${DONORS[$IDX]}"

echo "[$(date)] Array task $SLURM_ARRAY_TASK_ID -> donor $DONOR_ID"
echo "Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-1}  Mem: ${SLURM_MEM_PER_NODE:-?}M"

# Use the env that has scPloidy + nimble + GenomicRanges/Rsamtools installed.
ENV_NAME="${SCPLOIDY_ENV:-scploidy}"
echo "Using env: $ENV_NAME"

micromamba run -n "$ENV_NAME" Rscript Analysis/SingleCell/scripts/340_scploidy_hepatocyte.R
echo "[$(date)] Done $DONOR_ID"
