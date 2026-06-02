#!/bin/bash
#SBATCH --job-name=s_v3_09
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sex_v3_09_%j.out
#SBATCH --error=logs/sex_v3_09_%j.err

set -eo pipefail

# Module 09 — sex_v3 atlas refresh driver (light orchestrator).
# Submits sbatch chain: 27a -> 75 -> 217 with --dependency=afterok.
# Preflight asserts sex_deg_classification_v3.csv exists and 27a v3 fallback
# patch is present. Set DRY_RUN=FALSE to actually submit; default is DRY_RUN=TRUE.

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SDIR="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3"

cd "${SDIR}"
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== sex_v3 Module 09 — atlas refresh driver ==="
echo "Host:    $(hostname)"
echo "Started: $(date)"
echo "CPUs:    ${SLURM_CPUS_PER_TASK}"
echo "DRY_RUN: ${DRY_RUN:-TRUE}"
echo "Script:  ${SDIR}/09_atlas_refresh.R"

Rscript 09_atlas_refresh.R

echo "Finished: $(date)"
