#!/usr/bin/env bash
#SBATCH --job-name=hac_system_inference
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/70_infer_molecular_systems.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
ANALYSIS_CANDIDATE="${HAC_MOLECULAR_ANALYSIS_CANDIDATE:-${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-20260818T173348Z}"
ATLAS_SOURCE="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824"
: "${HAC_SYSTEM_INFERENCE_CANDIDATE:?HAC_SYSTEM_INFERENCE_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
if [[ -n "${HAC_SYSTEM_CONTEXT_ROSTER:-}" ]]; then
  : "${HAC_MOLECULAR_ANALYSIS_CANDIDATE:?Corrected inference requires recomputed molecular scores}"
  : "${HAC_SOURCE_ROOT:?Corrected inference requires recomputed fixed-loading axes}"
  : "${HAC_CORRECTED_STAGE_ROOT:?Corrected inference requires verified corrected model inputs}"
else
  export HAC_SOURCE_ROOT="${HAC_SOURCE_ROOT:-${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z}"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p "${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_systems_inference"
"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${ANALYSIS_CANDIDATE}" "${ATLAS_SOURCE}" "${HAC_SYSTEM_INFERENCE_CANDIDATE}"

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_SYSTEM_INFERENCE_CANDIDATE")
  scores <- fread(file.path(root, "tables", "participant_system_scores.tsv.gz"))
  cohort <- fread(file.path(root, "tables", "system_continuum_cohort_models.tsv"))
  meta <- fread(file.path(root, "tables", "system_continuum_meta_analysis.tsv"))
  stage <- fread(file.path(root, "tables", "system_fibrosis_stage_meta_analysis.tsv"))
  permutation <- fread(file.path(root, "tables", "system_within_stage_permutations.tsv"))
  membership <- fread(file.path(root, "tables", "system_inference_membership.tsv"))
  audit <- fread(file.path(root, "validation", "audit.tsv"))
  expected_permutations <- if (identical(Sys.getenv("HAC_ML_TEST_MODE"), "1")) {
    as.integer(Sys.getenv("HAC_ML_PERMUTATIONS"))
  } else {
    10000L
  }
  summary <- fread(file.path(root, "tables", "system_inference_summary.tsv"))
  expected_units <- as.integer(summary[metric == "participants_all_five_cohorts", value])
  if (!nzchar(Sys.getenv("HAC_SYSTEM_CONTEXT_ROSTER"))) stopifnot(expected_units == 844L)
  stopifnot(
    nrow(scores) == expected_units * 43L * 2L,
    uniqueN(scores$analysis_unit_id) == expected_units,
    uniqueN(scores$community_id) == 43L,
    nrow(cohort) == 2L * 2L * 2L * 43L,
    nrow(meta) == 2L * 2L * 43L,
    nrow(stage) == 2L * 43L,
    nrow(permutation) == 2L * 2L * 43L,
    all(permutation$permutation_replicates == expected_permutations),
    all(permutation$n_finite_permutations == expected_permutations),
    nrow(membership) == 43L,
    all(audit$passed)
  )
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
