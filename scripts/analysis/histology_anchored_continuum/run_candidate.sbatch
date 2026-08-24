#!/usr/bin/env bash
#SBATCH --job-name=hac_continuum
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/slurm-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

if [[ -z "${HAC_OUT_ROOT:-}" ]]; then
  echo "HAC_OUT_ROOT is unset" >&2
  exit 1
fi
if [[ ! -d "${HAC_OUT_ROOT}/external/kamzolas_v1" ]]; then
  echo "Frozen Kamzolas v1 checkout is missing from HAC_OUT_ROOT" >&2
  exit 1
fi

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export KAMZOLAS_V1_ROOT="${HAC_OUT_ROOT}/external/kamzolas_v1"
export HAC_RESOURCE_GENE_IDS_PATH="${HAC_OUT_ROOT}/inputs/resource_gene_ids.tsv"
export HAC_SCRIPT_DIR="${SCRIPT_DIR}"
export HAC_DGE_PATH="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/arms/F_five/results/integration/merged_dge.rds"
export HAC_MANIFEST_PATH="${PROJECT_ROOT}/figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis/stage_extensions/five_cohort_sample_manifest.tsv"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "${HAC_OUT_ROOT}/tests"

"${R_BIN}" --vanilla -e '
files <- list.files(Sys.getenv("HAC_SCRIPT_DIR"), pattern = "[.]R$", recursive = TRUE,
                    full.names = TRUE)
stopifnot(length(files) >= 10L)
for (path in files) parse(file = path, keep.source = FALSE)
cat("ALL_R_FILES_PARSE\t", length(files), "\n", sep = "")
' 2>&1 | tee "${HAC_OUT_ROOT}/tests/r_parse.log"

"${R_BIN}" --vanilla "${SCRIPT_DIR}/tests/run_tests.R" 2>&1 |
  tee "${HAC_OUT_ROOT}/tests/contract_tests.log"

"${R_BIN}" --vanilla "${SCRIPT_DIR}/00_prepare_inputs.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/01_reproduce_kamzolas_v1.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/02_build_unsupervised_axes.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/03_build_projection_rf_axes.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/04_score_programs_and_validate.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/05_validate_paired_gse193066.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/06_evaluate_gates.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/07_render_candidate_figures.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/08_assemble_report.R"
python3 "${SCRIPT_DIR}/09_validate_candidate.py"

echo "HISTOLOGY_ANCHORED_CONTINUUM_COMPLETE ${HAC_OUT_ROOT}"
