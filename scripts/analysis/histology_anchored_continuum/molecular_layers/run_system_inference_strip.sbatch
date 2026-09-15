#!/usr/bin/env bash
#SBATCH --job-name=hac_system_sigfig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/71_render_system_inference_strip.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
REFINED_SOURCE="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_pooled_narrow_20260824"
INFERENCE_CANDIDATE="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/hac-system-inference-20260824T232135Z"
: "${HAC_SYSTEM_INFERENCE_FIGURE_CANDIDATE:?HAC_SYSTEM_INFERENCE_FIGURE_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${REFINED_SOURCE}" "${INFERENCE_CANDIDATE}" \
  "${HAC_SYSTEM_INFERENCE_FIGURE_CANDIDATE}"

PANEL="${HAC_SYSTEM_INFERENCE_FIGURE_CANDIDATE}/panels/s4_molecular_systems_pooled_continuum_windows_with_inference.pdf"
[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"
! pdftotext "${PANEL}" - | grep -q 'maxT'
! pdftotext "${PANEL}" - | grep -q 'BH 33/43'
! pdftotext "${PANEL}" - | grep -q 'Delta'
for label in $(seq -w 1 43); do
  pdftotext "${PANEL}" - | grep -q "S${label}"
done

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_SYSTEM_INFERENCE_FIGURE_CANDIDATE")
  strip <- fread(file.path(root, "source_tables", "system_inference_strip.tsv"))
  audit <- fread(file.path(root, "provenance", "audit.tsv"))
  figures <- fread(file.path(root, "figure_manifest.tsv"))
  stopifnot(
    nrow(strip) == 86L,
    uniqueN(strip$community_id) == 43L,
    strip[inference_column == "BH", sum(supported)] == 33L,
    strip[inference_column == "maxT", sum(supported)] == 18L,
    strip[inference_column == "BH", sum(rendered_in_panel)] == 0L,
    strip[inference_column == "maxT", sum(rendered_in_panel)] == 18L,
    strip[inference_column == "maxT" & rendered_in_panel, all(render_symbol == "*")],
    nrow(figures) == 1L,
    !figures$window_values_changed,
    !figures$display_order_changed,
    !figures$delta_tile_displayed,
    figures$delta_order_retained,
    !figures$embedded_explanatory_text,
    !figures$inferential_model_changed,
    !figures$windows_inferentially_tested,
    !figures$promoted,
    all(audit$passed)
  )
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
printf '%s\n' "candidate_ready_for_visual_review" > \
  "${HAC_SYSTEM_INFERENCE_FIGURE_CANDIDATE}/READY_FOR_REVIEW.ok"
