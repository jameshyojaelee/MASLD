#!/usr/bin/env bash
#SBATCH --job-name=hac_pathway_overview
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/65_render_pathway_effect_overviews.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_PATHWAY_COMPARISON_SOURCE:?HAC_PATHWAY_COMPARISON_SOURCE is unset}"
: "${HAC_PATHWAY_OVERVIEW_CANDIDATE:?HAC_PATHWAY_OVERVIEW_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${HAC_PATHWAY_COMPARISON_SOURCE}" "${HAC_PATHWAY_OVERVIEW_CANDIDATE}"

while IFS=$'\t' read -r figure_id path rest; do
  [[ "${figure_id}" == "figure_id" ]] && continue
  [[ "$(pdfinfo "${path}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${path}"
  ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${path}"
done < "${HAC_PATHWAY_OVERVIEW_CANDIDATE}/figure_manifest.tsv"

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_PATHWAY_OVERVIEW_CANDIDATE")
  summary <- fread(file.path(root, "source_tables", "pathway_family_summary.tsv"))
  expected_family <- c(50L, 658L, 1787L, 7583L, 1855L, 1042L)
  expected_estimable <- c(44L, 260L, 1027L, 2811L, 720L, 467L)
  stopifnot(nrow(summary) == 6L,
            identical(summary$complete_family, expected_family),
            identical(summary$coordinate_estimable, expected_estimable),
            all(summary$coordinate_estimable + summary$coordinate_unestimable == summary$complete_family))
  figures <- fread(file.path(root, "figure_manifest.tsv"))
  stopifnot(nrow(figures) == 7L, all(file.exists(figures$path)),
            !any(figures$result_adaptive_feature_selection))
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
