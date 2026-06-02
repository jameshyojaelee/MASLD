#!/bin/bash
#SBATCH --job-name=ggplot
#SBATCH --qos=nslab
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_s2_2026-05-29/figS_loo_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_s2_2026-05-29/figS_loo_%j.err
# =====================================================================
# run_figS_loo_mash_masl.sh  (2026-05-29)
# Builds figS_lfc_sensitivity_loo_mash_masl.R from the freshly rebuilt
# STAR -s 2 LOO folds. Chained afterok on the LOO array job.
# =====================================================================
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u
cd "${PROJECT}"

echo "=== figS LOO mash/masl ==="
echo "host=$(hostname)  date=$(date)  JOB=${SLURM_JOB_ID}"

Rscript "${PROJECT}/scripts/figures/figS_lfc_sensitivity_loo_mash_masl.R"
rc=$?

echo "=== figure finished rc=${rc}  date=$(date) ==="
exit ${rc}
