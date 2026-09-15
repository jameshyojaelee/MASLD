#!/usr/bin/env bash
#SBATCH --job-name=hac_fig4h_systems
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/72_render_main_system_stage_continuum.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ATLAS_SOURCE:?HAC_ATLAS_SOURCE is unset}"
: "${HAC_SYSTEM_SUPPORT_SOURCE:?HAC_SYSTEM_SUPPORT_SOURCE is unset}"
: "${HAC_SYSTEM_STRESS_SOURCE:?HAC_SYSTEM_STRESS_SOURCE is unset}"
: "${HAC_FIG4H_CANDIDATE:?HAC_FIG4H_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${HAC_ATLAS_SOURCE}" "${HAC_SYSTEM_SUPPORT_SOURCE}" \
  "${HAC_SYSTEM_STRESS_SOURCE}" "${HAC_FIG4H_CANDIDATE}"

PANEL="${HAC_FIG4H_CANDIDATE}/panels/fig4h_molecular_systems_stage_vs_continuum.pdf"
[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_FIG4H_CANDIDATE")
  m <- fread(file.path(root, "figure_manifest.tsv"))
  x <- fread(file.path(root, "source_tables", "fig4h_system_stage_continuum.tsv"))
  a <- fread(file.path(root, "provenance", "audit.tsv"))
  stopifnot(
    nrow(m) == 1L, m$callout == "4H", m$complete_system_family == "43/43",
    m$marked_systems == "1/43", m$connected_null_supported == "0/43",
    m$joint_strict_supported == "35/43", m$four_test_maxt_supported == "18/43",
    m$same_direction == "36/43",
    !m$geometry_used, !m$window_values_used, nrow(x) == 43L,
    sum(x$strict_maxt_supported) == 18L, sum(x$competitive_marked) == 1L,
    sum(x$connected_null_supported) == 0L,
    x[competitive_marked == TRUE, community_id] == "system_01",
    sum(x$fixed_outcome_blind_label) == 15L, all(a$passed)
  )
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
