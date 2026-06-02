#!/bin/bash
#SBATCH --job-name=mcp_dialogue_mvp
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --array=0-5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/430_mvp_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/430_mvp_%A_%a.err

# MVP DIALOGUE runs: 2 phenotypes x 3 k = 6 combinations (array 0-5)
set -euo pipefail

PHENOTYPES=(disease_stage_numeric condition_binary_num)
K_LIST=(3 5 7)

idx=${SLURM_ARRAY_TASK_ID}
pheno_idx=$((idx / ${#K_LIST[@]}))
k_idx=$((idx % ${#K_LIST[@]}))
PHENO=${PHENOTYPES[$pheno_idx]}
K=${K_LIST[$k_idx]}

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate celltype_bio

echo "[430.mvp] array=${idx} pheno=${PHENO} k=${K}"
Rscript 430_dialogue_run.R --phenotype "${PHENO}" --k "${K}" \
  || echo "[430.mvp] phenotype ${PHENO} k=${K} failed (see logs)"

echo "[430.mvp] DONE."
