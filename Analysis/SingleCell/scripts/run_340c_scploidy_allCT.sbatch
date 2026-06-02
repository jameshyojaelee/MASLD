#!/bin/bash
#SBATCH --job-name=scpl_allCT
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-18%6
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ploidy/logs/scploidy_allCT_%A_%a.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ploidy/logs/scploidy_allCT_%A_%a.err
set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
source ~/.bashrc 2>/dev/null || true
DONORS=(D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 D14 D15 D16 D17 D18)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
export DONOR_ID="${DONORS[$IDX]}"
echo "[$(date)] Array task $SLURM_ARRAY_TASK_ID -> donor $DONOR_ID"
micromamba run -n scploidy Rscript Analysis/SingleCell/scripts/340c_scploidy_all_celltypes.R
