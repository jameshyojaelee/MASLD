#!/bin/bash
#SBATCH --job-name=mcp_dialogue_1khvg
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/430b_dialogue_1khvg_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/430b_dialogue_1khvg_%j.err

# DIALOGUE rehab rerun on 1K HVG/CT pre-filtered pseudobulk.
# Run only the MASLD-relevant phenotype (disease_stage_numeric) x k=5 to test convergence.
set -euo pipefail

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
set +u
micromamba activate rnaseq
set -u

echo "[430b] host=$(hostname) step 1 -- HVG pre-filter"
python 412b_hvg_filter_dialogue.py

echo "[430b] step 2 -- switch to celltype_bio and run DIALOGUE against 1K HVG dir"
set +u
micromamba deactivate
micromamba activate celltype_bio
set -u

# Override in_dir via a tiny wrapper: temporarily point the R script at the 1khvg dir.
# We pass --outdir and let the script read from the in_dir hardcoded to dialogue_pseudobulk.
# Since the R script hardcodes in_dir, we symlink 1khvg over it at runtime using a bash trick:
# copy the R to a temp location and sed in the new in_dir path.

TMP_R="/tmp/430_dialogue_run_1khvg_${SLURM_JOB_ID}.R"
sed 's|dialogue_pseudobulk|dialogue_pseudobulk_1khvg|' 430_dialogue_run.R > "${TMP_R}"

PHENOTYPES=(disease_stage_numeric condition_binary_num NAS fibrosis_stage)
K_LIST=(3 5 7)

for PHENO in "${PHENOTYPES[@]}"; do
  for K in "${K_LIST[@]}"; do
    OUTDIR="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/dialogue/rerun_1khvg/${PHENO}_k${K}"
    mkdir -p "${OUTDIR}"
    echo "[430b] phenotype=${PHENO} k=${K} outdir=${OUTDIR}"
    timeout 4h Rscript "${TMP_R}" --phenotype "${PHENO}" --k "${K}" --outdir "${OUTDIR}" \
      || echo "[430b] WARN phenotype=${PHENO} k=${K} failed or timed out (4h cap per combo)"
  done
done

rm -f "${TMP_R}"
echo "[430b] DONE."
