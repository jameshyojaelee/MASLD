#!/usr/bin/env bash
#SBATCH --job-name=hac_ml_freeze
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

: "${HAC_ML_OUT_ROOT:?HAC_ML_OUT_ROOT is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HAC_ML_SCRIPT_DIR="${SCRIPT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${R_BIN}" --vanilla -e '
files <- list.files(Sys.getenv("HAC_ML_SCRIPT_DIR"), pattern = "[.]R$", recursive = TRUE,
                    full.names = TRUE)
stopifnot(length(files) >= 5L)
for (path in files) parse(file = path, keep.source = FALSE)
cat("ALL_MOLECULAR_R_FILES_PARSE\t", length(files), "\n", sep = "")
' 2>&1 | tee "${HAC_ML_OUT_ROOT}/r_parse.log"

"${R_BIN}" --vanilla "${SCRIPT_DIR}/tests/run_tests.R" 2>&1 |
  tee "${HAC_ML_OUT_ROOT}/contract_tests.log"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/01_freeze_inputs.R"

echo "HAC_MOLECULAR_PREFLIGHT_COMPLETE ${HAC_ML_OUT_ROOT}"
