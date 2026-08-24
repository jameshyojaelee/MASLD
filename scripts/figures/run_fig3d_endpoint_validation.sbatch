#!/bin/bash
#SBATCH --job-name=fig3d_endpoint
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_endpoint_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig3d_endpoint_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set a fresh candidate output directory}"
: "${STAGE_RELEASE_ROOT:?Set the validated five-cohort release directory}"
: "${STAGE_EXTENSION_ROOT:?Set the synchronized stage-extension directory}"
test -d "${FIGURE_CANDIDATE_ROOT}"
test -d "${STAGE_RELEASE_ROOT}"
test -d "${STAGE_EXTENSION_ROOT}"
test ! -e "${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf"

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${project_root}"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/figures/fig3d_endpoint_validation.R"

analysis="${FIGURE_CANDIDATE_ROOT}/analysis/endpoint_validation"
panel="${FIGURE_CANDIDATE_ROOT}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf"
supp="${FIGURE_CANDIDATE_ROOT}/supplementary/figureS3/panels/figs3_endpoint_deg_concordance.pdf"
all_sample_supp="${FIGURE_CANDIDATE_ROOT}/supplementary/figureS3/panels/figs3d_pca_all_samples.pdf"
figs3c="${FIGURE_CANDIDATE_ROOT}/supplementary/figureS3/panels/figs3c_cohort_robustness.pdf"

# Re-index the immutable all-sample PCA candidate as supplementary. The endpoint
# script rebuilds supplementary Figure S3C from synchronized per-cohort effects
# plus fresh leave-one-cohort-out fits.
all_sample_root="${project_root}/figures/candidates/fig3d-all-samples-2026-08-14-v2"
cp "${all_sample_root}/figure3/panels/fig3d_pca_fibrosis_gradient.pdf" \
  "${all_sample_supp}"
cp "${all_sample_root}/source_tables/fig3d_pca_coordinates.tsv" \
  "${FIGURE_CANDIDATE_ROOT}/source_tables/figs3d_all_sample_pca_coordinates.tsv"
cp "${all_sample_root}/source_tables/fig3d_pca_census.tsv" \
  "${FIGURE_CANDIDATE_ROOT}/source_tables/figs3d_all_sample_pca_census.tsv"
test "$(awk -F '\t' -v p="${STAGE_RELEASE_ROOT}/VALIDATED.json" '$2 == p {print $2}' \
  "${STAGE_EXTENSION_ROOT}/input_manifest.tsv")" = "${STAGE_RELEASE_ROOT}/VALIDATED.json"

"${rnaseq_env}/bin/Rscript" --vanilla - <<'RS'
library(data.table)
root <- Sys.getenv("FIGURE_CANDIDATE_ROOT")
roster <- fread(file.path(root, "analysis", "endpoint_validation", "endpoint_roster.tsv"))
audit <- fread(file.path(root, "analysis", "endpoint_validation", "endpoint_model_audit.tsv"))
de <- fread(file.path(root, "analysis", "endpoint_validation", "endpoint_de_results.tsv.gz"))
stopifnot(nrow(roster) == 418L, uniqueN(roster$dataset) == 4L,
          sum(roster$endpoint_code == "Advanced") == 303L,
          sum(roster$endpoint_code == "StrictControl") == 115L,
          nrow(de) == 23370L, audit$bh_family_size == 23370L,
          identical(audit$canonical_write, FALSE))
fig3c <- fread(file.path(root, "source_tables", "figs3c_cohort_robustness.tsv"))
stopifnot(nrow(fig3c) == 5L, all(fig3c$canonical_genes == 1347L),
          all(is.finite(fig3c$same_direction_fraction)),
          all(is.finite(fig3c$all_gene_spearman)))
cat("FIG3D_ENDPOINT_CONTRACT PASS: 418 participants, 4 matched cohorts, 23,370 genes\n")
RS

for pdf in "${panel}" "${supp}" "${all_sample_supp}" "${figs3c}"; do
  test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
  if strings "${pdf}" | grep -q '/Subtype /Type3'; then
    echo "Type 3 font detected: ${pdf}" >&2
    exit 1
  fi
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
done

mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/fig3d_endpoint_pca.png" "${panel}"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/figs3_endpoint_deg_concordance.png" "${supp}"

find "${FIGURE_CANDIDATE_ROOT}" -type f ! -name input_output_manifest.tsv -print0 | sort -z | \
  xargs -0 sha256sum > "${FIGURE_CANDIDATE_ROOT}/input_output_manifest.tsv"
python3 "${project_root}/scripts/manuscript/validate_resource_scope.py"
