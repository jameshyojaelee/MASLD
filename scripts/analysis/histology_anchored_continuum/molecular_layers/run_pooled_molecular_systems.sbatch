#!/usr/bin/env bash
#SBATCH --job-name=hac_systems_pool
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/68_render_pooled_molecular_systems.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
ATLAS_SOURCE="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824"
: "${HAC_POOLED_SYSTEMS_CANDIDATE:?HAC_POOLED_SYSTEMS_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" "${ATLAS_SOURCE}" "${HAC_POOLED_SYSTEMS_CANDIDATE}"

MANIFEST="${HAC_POOLED_SYSTEMS_CANDIDATE}/figure_manifest.tsv"
PANEL="$(awk -F $'\t' 'NR==2 {print $2}' "${MANIFEST}")"
[[ "$(pdfinfo "${PANEL}" | awk '/^Pages:/ {print $2}')" == "1" ]]
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PANEL}"
! grep -aEq '/Subtype[[:space:]]*/Type3|/FontType[[:space:]]*3' "${PANEL}"

for label in $(seq -w 1 43); do
  pdftotext "${PANEL}" - | grep -q "system_${label}"
done

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_POOLED_SYSTEMS_CANDIDATE")
  pooled <- fread(file.path(root, "source_tables", "pooled_molecular_system_windows.tsv"))
  audit <- fread(file.path(root, "provenance", "audit.tsv"))
  figures <- fread(file.path(root, "figure_manifest.tsv"))
  stopifnot(
    nrow(pooled) == 43L * 9L,
    uniqueN(pooled$community_id) == 43L,
    uniqueN(pooled$window_id) == 9L,
    all(pooled$n_cohorts == 5L),
    all(pooled$visualization_only),
    all(audit$passed),
    nrow(figures) == 1L,
    figures$figure_id == "pooled_continuum_windows",
    figures$pooling_weight == "equal_cohort",
    file.exists(figures$path)
  )
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
printf '%s\n' "candidate_ready_for_visual_review" > \
  "${HAC_POOLED_SYSTEMS_CANDIDATE}/READY_FOR_REVIEW.ok"
