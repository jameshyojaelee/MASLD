#!/usr/bin/env bash
#SBATCH --job-name=hac_5cohort_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ANALYSIS_CANDIDATE:?HAC_ANALYSIS_CANDIDATE is unset}"
: "${HAC_FIVE_COHORT_CANDIDATE:?HAC_FIVE_COHORT_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HAC_SOURCE_ROOT="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT_DIR}/64_render_five_cohort_context.R')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT_DIR}/64_render_five_cohort_context.R" \
  "${HAC_ANALYSIS_CANDIDATE}" "${HAC_FIVE_COHORT_CANDIDATE}"

while IFS=$'\t' read -r callout filename path rest; do
  [[ "${callout}" == "callout" ]] && continue
  [[ "$(pdfinfo "${path}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${path}"
  if grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${path}"; then
    echo "Type 3 font detected: ${path}" >&2
    exit 1
  fi
done < "${HAC_FIVE_COHORT_CANDIDATE}/figure_manifest.tsv"

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_FIVE_COHORT_CANDIDATE")
  m <- fread(file.path(root, "figure_manifest.tsv"))
  stopifnot(nrow(m) == 13L, all(file.exists(m$path)))
  five <- m[callout %in% c("4F", "S3P", "S3T", "S4F", "S4H", "S4I", "S4M", "S4O")]
  stage <- m[callout %in% c("S3I", "S3K", "S3U", "S4G", "S4N")]
  stopifnot(all(five$descriptive_cohorts == 5L), all(stage$descriptive_cohorts == 4L),
            m[callout == "4F", primary_inference_cohorts] == 2L,
            !m[callout == "4F", embedded_explanatory_text],
            !m[callout == "4F", forest_panel_displayed],
            !m[callout == "4F", participant_size_legend_displayed])
  figure4f_text <- system2("pdftotext", c(m[callout == "4F", path], "-"), stdout = TRUE)
  stopifnot(!any(grepl("Source-overlap|Cohort PC1|Fixed projection|Participants", figure4f_text)))
  h <- fread(file.path(root, "source_tables", "all117_hotspot_five_cohort_windows.tsv.gz"))
  stopifnot(uniqueN(h$program_uid) == 117L, uniqueN(h$dataset) == 5L,
            uniqueN(h$window_id) == 9L, nrow(h) == 117L * 5L * 9L,
            all(h$visualization_only))
  five_tables <- c(
    "focal_hotspot_five_cohort_windows.tsv",
    "all10_nmf_five_cohort_windows.tsv",
    "six_hallmark_five_cohort_windows.tsv",
    "four_tf_five_cohort_windows.tsv",
    "five_gene_five_cohort_windows.tsv"
  )
  for (name in five_tables) {
    x <- fread(file.path(root, "source_tables", name))
    stopifnot(uniqueN(x$dataset) == 5L, setequal(unique(x$dataset), c(
      "GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"
    )))
  }
'

printf 'candidate_ready\n' > "${HAC_FIVE_COHORT_CANDIDATE}/READY_FOR_REVIEW.ok"

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
