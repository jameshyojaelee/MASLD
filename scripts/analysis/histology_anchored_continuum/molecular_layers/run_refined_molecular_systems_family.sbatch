#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_polish
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/69_render_refined_molecular_systems_family.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
ATLAS_SOURCE="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824"
POOLED_SOURCE="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_pooled_20260824_v2"
: "${HAC_REFINED_SYSTEMS_CANDIDATE:?HAC_REFINED_SYSTEMS_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${ATLAS_SOURCE}" "${POOLED_SOURCE}" "${HAC_REFINED_SYSTEMS_CANDIDATE}"

MANIFEST="${HAC_REFINED_SYSTEMS_CANDIDATE}/figure_manifest.tsv"
while IFS=$'\t' read -r figure_id path rest; do
  [[ "${figure_id}" == "figure_id" ]] && continue
  [[ "$(pdfinfo "${path}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${path}"
  ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${path}"
done < "${MANIFEST}"

for panel in \
  "${HAC_REFINED_SYSTEMS_CANDIDATE}/panels/s4_molecular_systems_five_cohort_sliding_windows.pdf" \
  "${HAC_REFINED_SYSTEMS_CANDIDATE}/panels/s4_molecular_systems_pooled_continuum_windows.pdf"; do
  for label in $(seq -w 1 43); do
    pdftotext "${panel}" - | grep -q "S${label}"
  done
done

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_REFINED_SYSTEMS_CANDIDATE")
  order <- fread(file.path(root, "source_tables", "pooled_system_late_minus_early_order.tsv"))
  display <- fread(file.path(root, "source_tables", "pooled_system_heatmap_display_values.tsv"))
  audit <- fread(file.path(root, "provenance", "audit.tsv"))
  figures <- fread(file.path(root, "figure_manifest.tsv"))
  stopifnot(
    nrow(order) == 43L,
    uniqueN(order$community_id) == 43L,
    identical(order$display_rank, 1:43),
    all(diff(order$late_minus_early) <= 0),
    nrow(display) == 43L * 10L,
    uniqueN(display$community_id) == 43L,
    uniqueN(display$display_position) == 10L,
    nrow(figures) == 6L,
    sum(figures$result_adaptive_display_order) == 2L,
    !any(figures$data_values_changed),
    !any(figures$geometry_changed),
    !any(figures$inferential_use_of_display_order),
    !any(figures$embedded_explanatory_text),
    all(file.exists(figures$path)),
    all(audit$passed)
  )
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
printf '%s\n' "candidate_ready_for_visual_review" > \
  "${HAC_REFINED_SYSTEMS_CANDIDATE}/READY_FOR_REVIEW.ok"
