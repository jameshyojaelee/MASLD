#!/usr/bin/env bash
#SBATCH --job-name=hac_system_stress
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers/73_stress_test_molecular_systems.R"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
ANALYSIS_CANDIDATE="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-20260818T173348Z"
ATLAS_SOURCE="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824"
INFERENCE_CANDIDATE="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/hac-system-inference-20260824T232135Z"
: "${HAC_SYSTEM_STRESS_CANDIDATE:?HAC_SYSTEM_STRESS_CANDIDATE is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HAC_SOURCE_ROOT="${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-20260818T024923Z"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export HAC_ML_WORKERS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${PROJECT_ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_systems_stress"
"${R_BIN}" --vanilla -e "parse(file='${SCRIPT}')" >/dev/null
"${R_BIN}" --vanilla "${SCRIPT}" \
  "${ANALYSIS_CANDIDATE}" "${ATLAS_SOURCE}" "${INFERENCE_CANDIDATE}" \
  "${HAC_SYSTEM_STRESS_CANDIDATE}"

"${R_BIN}" --vanilla -e '
  library(data.table)
  root <- Sys.getenv("HAC_SYSTEM_STRESS_CANDIDATE")
  joint <- fread(file.path(root, "tables", "system_joint_permutations.tsv"))
  competitive <- fread(file.path(root, "tables", "system_competitive_null.tsv"))
  draws <- fread(file.path(root, "tables", "system_competitive_null_draws.tsv.gz"))
  attenuation <- fread(file.path(root, "tables", "stage_continuum_attenuation.tsv"))
  redundancy <- fread(file.path(root, "tables", "system_redundancy.tsv"))
  membership <- fread(file.path(root, "tables", "system_stress_membership.tsv"))
  audit <- fread(file.path(root, "validation", "audit.tsv"))
  test_mode <- identical(Sys.getenv("HAC_ML_TEST_MODE"), "1")
  expected_permutations <- if (test_mode) {
    as.integer(Sys.getenv("HAC_ML_PERMUTATIONS", unset = "200"))
  } else {
    10000L
  }
  expected_draws <- if (test_mode) {
    as.integer(Sys.getenv("HAC_ML_COMPETITIVE_DRAWS", unset = "50"))
  } else {
    5000L
  }
  stopifnot(
    nrow(joint) == 43L * 2L,
    all(joint$permutation_replicates == expected_permutations),
    nrow(competitive) == 43L * 2L * 2L,
    all(competitive$requested_draws == expected_draws),
    all(competitive$finite_draws > 0L),
    nrow(draws) == sum(competitive$finite_draws),
    nrow(attenuation) == 43L * 2L * 2L,
    nrow(redundancy) == 3L,
    nrow(membership) == 43L,
    all(audit$passed)
  )
  cat("post-run assertions passed\n")
'

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
