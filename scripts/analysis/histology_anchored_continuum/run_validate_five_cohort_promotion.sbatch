#!/usr/bin/env bash
#SBATCH --job-name=hac_validate_5cohort
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
RENDERER="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/64_render_five_cohort_context.R"
PROMOTER="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/promote_five_cohort_context.py"
CANDIDATE="${PROJECT_ROOT}/figures/candidates/histology-continuum-five-cohort-context-20260824T175647Z"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

python3 -m json.tool "${PROJECT_ROOT}/config/method_exemptions.json" >/dev/null
python3 -m py_compile "${PROMOTER}"
"${R_BIN}" --vanilla -e "parse(file='${RENDERER}')" >/dev/null

for manifest in \
  "${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/manifests/continuum_five_cohort_revision_20260824.tsv" \
  "${PROJECT_ROOT}/figures/supplementary/figS03_bulk_transcriptomics/manifests/continuum_five_cohort_revision_20260824.tsv"; do
  while IFS=$'\t' read -r callout destination source old_sha new_sha descriptive inference; do
    [[ "${callout}" == "callout" ]] && continue
    pdf="${PROJECT_ROOT}/${destination}"
    [[ -f "${pdf}" && -f "${PROJECT_ROOT}/${source}" ]]
    [[ "$(sha256sum "${pdf}" | awk '{print $1}')" == "${new_sha}" ]]
    [[ "$(sha256sum "${PROJECT_ROOT}/${source}" | awk '{print $1}')" == "${new_sha}" ]]
    [[ "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" == "1" ]]
    gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
    ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${pdf}"
  done < "${manifest}"
done

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("MASLD_PROJECT_ROOT")
  candidate <- file.path(root, "figures/candidates/histology-continuum-five-cohort-context-20260824T175647Z")
  m <- fread(file.path(candidate, "figure_manifest.tsv"))
  stopifnot(nrow(m) == 13L, all(file.exists(m$path)))
  stopifnot(all(m[callout %in% c("4F", "S3P", "S3T", "S4F", "S4H", "S4I", "S4M", "S4O"), descriptive_cohorts] == 5L))
  stopifnot(all(m[callout %in% c("S3I", "S3K", "S3U", "S4G", "S4N"), descriptive_cohorts] == 4L))
  stopifnot(m[callout == "4F", primary_inference_cohorts] == 2L)
  expected_five <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
  for (name in c(
    "focal_hotspot_five_cohort_windows.tsv", "all10_nmf_five_cohort_windows.tsv",
    "six_hallmark_five_cohort_windows.tsv", "four_tf_five_cohort_windows.tsv",
    "five_gene_five_cohort_windows.tsv"
  )) {
    x <- fread(file.path(candidate, "source_tables", name))
    stopifnot(setequal(unique(x$dataset), expected_five), uniqueN(x$window_id) == 9L)
  }
  h <- fread(file.path(candidate, "source_tables", "all117_hotspot_five_cohort_windows.tsv.gz"))
  stopifnot(nrow(h) == 117L * 5L * 9L, uniqueN(h$program_uid) == 117L,
            setequal(unique(h$dataset), expected_five), all(h$visualization_only))
'

grep -Fq '64_render_five_cohort_context.R' \
  "${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/CANDIDATE_PANEL_INDEX.tsv"
grep -Fq '64_render_five_cohort_context.R' \
  "${PROJECT_ROOT}/figures/supplementary/figS03_bulk_transcriptomics/CANDIDATE_PANEL_INDEX.tsv"
python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
echo 'FIVE_COHORT_PROMOTION_VALIDATION PASS'
