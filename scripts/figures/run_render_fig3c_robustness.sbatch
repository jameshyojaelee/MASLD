#!/bin/bash
#SBATCH --job-name=fig3c_restore
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3c_render_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3c_render_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set a new candidate output directory}"
: "${STAGE_RELEASE_ROOT:?Set the validated five-cohort release directory}"
: "${STAGE_EXTENSION_ROOT:?Set the synchronized stage-extension directory}"
: "${ROBUSTNESS_SOURCE:?Set the validated Figure S3C source table}"
test ! -e "${FIGURE_CANDIDATE_ROOT}"
test -d "${STAGE_RELEASE_ROOT}"
test -d "${STAGE_EXTENSION_ROOT}"
test -f "${ROBUSTNESS_SOURCE}"
test "$(awk -F '\t' -v p="${STAGE_RELEASE_ROOT}/VALIDATED.json" '$2 == p {print $2}' \
  "${STAGE_EXTENSION_ROOT}/input_manifest.tsv")" = "${STAGE_RELEASE_ROOT}/VALIDATED.json"
mkdir -p "${FIGURE_CANDIDATE_ROOT}/source_tables"
cp "${ROBUSTNESS_SOURCE}" \
  "${FIGURE_CANDIDATE_ROOT}/source_tables/figs3c_cohort_robustness.tsv"
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export STAGE_EXTENSION_ROOT
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
"${rnaseq_env}/bin/Rscript" --vanilla \
  "${MASLD_PROJECT_ROOT}/scripts/figures/nas_fib_grid.R"
"${rnaseq_env}/bin/Rscript" --vanilla \
  "${MASLD_PROJECT_ROOT}/scripts/figures/render_fig3c_robustness.R"

main_pdf="${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3c_cohort_alluvial.pdf"
supp_pdf="${FIGURE_CANDIDATE_ROOT}/supplementary/figureS3/panels/figs3c_cohort_robustness.pdf"
for pdf in "${main_pdf}" "${supp_pdf}"; do
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  if strings "${pdf}" | grep -q '/Subtype /Type3'; then
    echo "Type 3 font detected: ${pdf}" >&2
    exit 1
  fi
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
done
test "$(pdfinfo "${main_pdf}" | awk '/^Page size:/ {print $3 " " $5}')" = "364 285"

"${rnaseq_env}/bin/Rscript" --vanilla - <<'RS'
library(data.table)
root <- Sys.getenv("FIGURE_CANDIDATE_ROOT")
alluvial <- fread(file.path(root, "source_tables", "fig3c_cohort_deg_replication.tsv"))
robustness <- fread(file.path(root, "source_tables", "figs3c_cohort_robustness.tsv"))
stopifnot(nrow(alluvial) == 25L, uniqueN(alluvial$dataset) == 5L,
          identical(sort(unique(alluvial$n_cohorts)), 1:5), all(alluvial$N > 0L),
          nrow(robustness) == 5L, all(robustness$canonical_genes == 1347L))
cat("FIG3C_RESTORE_CARDINALITY PASS\n")
RS

mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/fig3c_cohort_alluvial.png" "${main_pdf}"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/figs3c_cohort_robustness.png" "${supp_pdf}"
mkdir -p "${FIGURE_CANDIDATE_ROOT}/provenance"
{
  printf 'source_id\tpath\tsha256\n'
  printf 'validated_five_cohort_release\t%s\t%s\n' \
    "${STAGE_RELEASE_ROOT}/VALIDATED.json" \
    "$(sha256sum "${STAGE_RELEASE_ROOT}/VALIDATED.json" | cut -d' ' -f1)"
  printf 'synchronized_stage_extension\t%s\t%s\n' \
    "${STAGE_EXTENSION_ROOT}/cohort_disease_all_gene_results.tsv" \
    "$(sha256sum "${STAGE_EXTENSION_ROOT}/cohort_disease_all_gene_results.tsv" | cut -d' ' -f1)"
  printf 'cohort_robustness_source\t%s\t%s\n' \
    "${ROBUSTNESS_SOURCE}" "$(sha256sum "${ROBUSTNESS_SOURCE}" | cut -d' ' -f1)"
} > "${FIGURE_CANDIDATE_ROOT}/provenance/fig3c_input_manifest.tsv"
"${rnaseq_env}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "provenance", "sessionInfo.txt"))'
find "${FIGURE_CANDIDATE_ROOT}" -type f ! -name input_output_manifest.tsv -print0 | sort -z | \
  xargs -0 sha256sum > "${FIGURE_CANDIDATE_ROOT}/input_output_manifest.tsv"
python3 "${MASLD_PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
