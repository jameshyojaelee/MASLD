#!/bin/bash
#SBATCH --job-name=coloc_r2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_eligibility/coloc_r2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_eligibility/coloc_r2_%j.err

# Corrected COLOC release r2 (review item 2, Phase 2 step 5). Submit after all
# five 09b eligibility batches succeed:
#
#   sbatch --dependency=afterok:<b1>:<b2>:<b3>:<b4>:<b5> \
#     --export=ALL,R2_ROOT=<FM>/results/susie_coloc_r2_<UTC ts> \
#     src/perf/11b_assemble_coloc_r2.sbatch
#
# 1 assemble   refuses if the 09b logs show an eligibility mismatch, an R error or an
#              unscanned resubmission (set BATCH_JOBS=<id,id,...> to include
#              resubmitted batch jobs); then adopted master + rerun SuSiE fields ->
#              R2_ROOT/susie_coloc_all_gwas.csv (fails on any rerun/offline
#              disagreement; see qc/, incl. eligibility_rule_comparison.tsv and
#              call_changes.tsv)
# 2 reproduce  the master-in combine path on the ADOPTED master must rebuild the
#              adopted gene_level_coloc.csv (checked in step 5)
# 3 combine    07_combine_rerun.R on the r2 master -> R2_ROOT/gene_level_coloc.csv
# 4 tier-1/2   rebuild_tier12.py -> R2_ROOT/tier12/, copied to R2_ROOT/ so the
#              release has the adopted layout (master, gene-level, tier-1/2)
# 5 summarize  old-vs-new tables -> R2_ROOT/summary/ and R2_ROOT/r2_summary.json
set -eo pipefail
: "${R2_ROOT:?Set R2_ROOT to a new results/susie_coloc_r2_<UTC ts> directory}"
test ! -e "${R2_ROOT}"

set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
ADOPTED=${FM}/results/susie_coloc
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${FM}"

python src/perf/11_assemble_coloc_r2.py assemble --release-root "${R2_ROOT}" \
  ${BATCH_JOBS:+--batch-jobs "${BATCH_JOBS}"}

REPRO=${R2_ROOT}/qc/adopted_gene_level_rerun
mkdir "${REPRO}"
cp "${ADOPTED}/ld_contamination_clusters.csv" "${REPRO}/"
COLOC_DIR="${REPRO}" COLOC_MASTER_IN="${ADOPTED}/susie_coloc_all_gwas.csv" \
  Rscript src/perf/07_combine_rerun.R > "${REPRO}/combine.log" 2>&1

COLOC_DIR="${R2_ROOT}" COLOC_MASTER_IN="${R2_ROOT}/susie_coloc_all_gwas.csv" \
  Rscript src/perf/07_combine_rerun.R > "${R2_ROOT}/qc/combine.log" 2>&1

python src/perf/rebuild_tier12.py --out-dir "${R2_ROOT}/tier12" \
  --master "${R2_ROOT}/susie_coloc_all_gwas.csv" > "${R2_ROOT}/qc/rebuild_tier12.log" 2>&1
cp "${R2_ROOT}/tier12/gene_level_coloc_tier12.csv" "${R2_ROOT}/gene_level_coloc_tier12.csv"

python src/perf/11_assemble_coloc_r2.py summarize --release-root "${R2_ROOT}"

Rscript -e 'sessionInfo()' > "${R2_ROOT}/sessionInfo.txt"
python -m pip freeze > "${R2_ROOT}/pip_freeze.txt"
sha256sum "${R2_ROOT}"/susie_coloc_all_gwas.csv "${R2_ROOT}"/gene_level_coloc.csv \
  "${R2_ROOT}"/gene_level_coloc_tier12.csv "${R2_ROOT}"/susie_coloc_signal_pairs.tsv.gz \
  "${R2_ROOT}"/r2_summary.json > "${R2_ROOT}/output_sha256.tsv"
echo "COLOC_R2_COMPLETE ${R2_ROOT}"
