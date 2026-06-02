#!/bin/bash
#SBATCH --job-name=stage_ccc_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/stage_ccc_regen_%j.out
#SBATCH --error=scripts/figures/logs/stage_ccc_regen_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOG_DIR="${BASE}/scripts/figures/logs"
mkdir -p "${LOG_DIR}"

micromamba run -n rnaseq Rscript "${BASE}/scripts/figures/figM_stage_ccc_main.R"
echo "[done] figM_stage_ccc_main.pdf"

micromamba run -n rnaseq Rscript "${BASE}/scripts/figures/figS_stage_ccc_supp.R"
echo "[done] figS_stage_ccc_supp.pdf"

micromamba run -n rnaseq Rscript "${BASE}/Analysis/SingleCell/scripts/v2/07e_figures_v3.R"
echo "[done] figS_stage_ccc_hardened.pdf"

micromamba run -n rnaseq Rscript "${BASE}/Analysis/SingleCell/scripts/352_lr_stage_trajectories.R"
echo "[done] figS_stage_ccc_LR_trajectories.pdf"

echo "[all done] stage_ccc figures regenerated"
