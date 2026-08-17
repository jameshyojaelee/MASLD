#!/bin/bash
#SBATCH --job-name=susieR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/refit_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/refit_%A_%a.err
#
# Targeted re-fit of the 204 eGenes that stopped converging at max_iter=500 when
# the eQTL SuSiE fits were regenerated with LD_PANEL pinned to polyfun.
#
# WHY IT MATTERS OUT OF PROPORTION TO 204: each eQTL fit is reused by ALL 50
# GWAS studies, so one non-converged gene costs ~50 gene-GWAS pairs. Measured:
# 204 genes -> 10,197 bail events (49.98x amplification), 2,659 lost SuSiE
# results, 171 of them above PP.H4.susie 0.5, and a 736 -> 712 headline drop.
#
# SAFETY: writes to results/eqtl_susie_refit/ (staging). The production fits in
# results/eqtl_susie_polyfun/ are NOT touched; merging is a separate, verified
# step. 08 skips genes whose RDS already exists, and the staging dir starts
# empty, so every requested gene is genuinely re-fitted.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export LD_PANEL=polyfun                    # identical to how eqtl_susie_polyfun was built
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_refit"
export REFIT_GENES="${FM}/tmp_diag/nonconverged_204_genes.txt"
export SUSIE_MAX_ITER="${SUSIE_MAX_ITER:-5000}"
export CHR=${SLURM_ARRAY_TASK_ID}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
case "${EQTL_SUSIE_DIR}" in
  *_refit) : ;;
  *) echo "FATAL: EQTL_SUSIE_DIR='${EQTL_SUSIE_DIR}' is not the staging dir; refusing" >&2; exit 1 ;;
esac
echo "[$(date -u +%FT%TZ)] refit chr${CHR}  max_iter=${SUSIE_MAX_ITER}  -> ${EQTL_SUSIE_DIR}"
cd "${FM}"
Rscript "${FM}/src/perf/08_refit_nonconverged.R"
