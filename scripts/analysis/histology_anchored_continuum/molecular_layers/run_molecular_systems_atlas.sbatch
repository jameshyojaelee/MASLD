#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_atlas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ANALYSIS_CANDIDATE:?HAC_ANALYSIS_CANDIDATE is unset}"
: "${HAC_COMPARISON_CANDIDATE:?HAC_COMPARISON_CANDIDATE is unset}"
: "${HAC_SYSTEMS_ATLAS_CANDIDATE:?HAC_SYSTEMS_ATLAS_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HAC_SOURCE_ROOT="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

for script in \
  "${SCRIPT_DIR}/66_freeze_molecular_systems_geometry.R" \
  "${SCRIPT_DIR}/67_overlay_molecular_systems_atlas.R"; do
  "${R_BIN}" --vanilla -e "parse(file='${script}')" >/dev/null
done

"${R_BIN}" --vanilla "${SCRIPT_DIR}/66_freeze_molecular_systems_geometry.R" \
  "${HAC_ANALYSIS_CANDIDATE}" "${HAC_SYSTEMS_ATLAS_CANDIDATE}"

"${R_BIN}" --vanilla "${SCRIPT_DIR}/67_overlay_molecular_systems_atlas.R" \
  "${HAC_ANALYSIS_CANDIDATE}" "${HAC_COMPARISON_CANDIDATE}" \
  "${HAC_SYSTEMS_ATLAS_CANDIDATE}"

while IFS=$'\t' read -r figure_id path rest; do
  [[ "${figure_id}" == "figure_id" ]] && continue
  [[ "$(pdfinfo "${path}" | awk '/^Pages:/ {print $2}')" == "1" ]]
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${path}"
  ! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${path}"
done < "${HAC_SYSTEMS_ATLAS_CANDIDATE}/figure_manifest.tsv"

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_SYSTEMS_ATLAS_CANDIDATE")
  ledger <- fread(file.path(root, "source_tables", "complete_feature_testability_ledger.tsv.gz"))
  nodes <- fread(file.path(root, "geometry", "latent_coordinates.tsv.gz"))
  communities <- fread(file.path(root, "geometry", "community_registry.tsv"))
  components <- fread(file.path(root, "geometry", "graph_component_registry.tsv"))
  representatives <- fread(file.path(root, "geometry", "community_top10_representatives.tsv.gz"))
  sensitivity <- fread(file.path(root, "geometry", "topology_sensitivity.tsv"))
  windows <- fread(file.path(root, "overlays", "all_node_five_cohort_fixed_windows.tsv.gz"))
  summary <- fread(file.path(root, "overlays", "atlas_summary.tsv"))
  effect_sensitivity <- fread(file.path(root, "overlays", "effect_aggregation_method_summary.tsv"))
  window_sensitivity <- fread(file.path(root, "overlays", "window_aggregation_method_summary.tsv"))
  label_manifest <- fread(file.path(root, "overlays", "rendered_label_manifest.tsv"))
  figures <- fread(file.path(root, "figure_manifest.tsv"))
  geometry_audit <- fread(file.path(root, "provenance", "outcome_blind_audit.tsv"))
  overlay_audit <- fread(file.path(root, "provenance", "overlay_audit.tsv"))
  stopifnot(
    nrow(ledger) == 13092L,
    sum(ledger$geometry_eligible) == 5445L,
    nrow(nodes) == 5445L,
    uniqueN(nodes$feature_id) == 5445L,
    uniqueN(nodes[nodes$collection == "hotspot", feature_id]) == 116L,
    uniqueN(nodes[nodes$collection != "hotspot", feature_id]) == 5329L,
    nrow(communities) == 43L,
    nrow(components) == 4L,
    components[component_role == "main_system_component", n_features] == 5439L,
    sum(components[component_role == "excluded_microcomponent", n_features]) == 6L,
    sum(nodes$component_role == "main_system_component") == 5439L,
    sum(nodes$component_role == "excluded_microcomponent") == 6L,
    all(is.na(nodes[component_role == "excluded_microcomponent", community_id])),
    setequal(representatives$community_id, communities$community_id),
    all(representatives[, .N, by = community_id]$N <= 10L),
    nrow(sensitivity) == 5L,
    all(is.finite(sensitivity$adjusted_rand_vs_baseline)),
    !any(sensitivity$outcome_values_read),
    nrow(windows) == 5445L * 5L * 9L,
    uniqueN(windows$dataset) == 5L,
    uniqueN(windows$window_id) == 9L,
    nrow(figures) == 5L,
    all(file.exists(figures$path)),
    all(figures$geometry_outcome_blind),
    !any(figures$result_adaptive_geometry),
    !any(figures$promoted),
    nrow(summary) == 15L,
    setequal(effect_sensitivity$aggregation_method,
             c("collection_balanced", "feature_balanced_median")),
    setequal(window_sensitivity$aggregation_method,
             c("collection_balanced", "feature_balanced_median")),
    all(label_manifest$intended_and_rendered),
    all(figures$rendered_label_count ==
          label_manifest[, uniqueN(community_id), by = figure_id][
            match(figures$figure_id, figure_id), V1]),
    all(geometry_audit$passed),
    all(overlay_audit$passed)
  )
  for (current_figure_id in unique(label_manifest$figure_id)) {
    path <- figures[figure_id == current_figure_id, path]
    pdf_text <- system2("pdftotext", c(path, "-"), stdout = TRUE, stderr = TRUE)
    expected_labels <- label_manifest[figure_id == current_figure_id, rendered_label]
    stopifnot(all(vapply(
      expected_labels,
      function(label) any(grepl(label, pdf_text, fixed = TRUE)),
      logical(1)
    )))
  }
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
printf '%s\n' "candidate_ready_for_visual_review" > \
  "${HAC_SYSTEMS_ATLAS_CANDIDATE}/READY_FOR_REVIEW.ok"
