#!/usr/bin/env bash
#SBATCH --job-name=hac_preflight
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/preflight-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

if [[ -z "${HAC_PREFLIGHT_OUT:-}" ]]; then
  echo "HAC_PREFLIGHT_OUT is unset" >&2
  exit 1
fi
if [[ -e "${HAC_PREFLIGHT_OUT}" ]]; then
  echo "Refusing to overwrite preflight namespace: ${HAC_PREFLIGHT_OUT}" >&2
  exit 1
fi

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HAC_SCRIPT_DIR="${SCRIPT_DIR}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir -p "${HAC_PREFLIGHT_OUT}"

"${R_BIN}" --vanilla -e '
files <- list.files(Sys.getenv("HAC_SCRIPT_DIR"), pattern = "[.]R$", recursive = TRUE,
                    full.names = TRUE)
stopifnot(length(files) >= 10L)
for (path in files) parse(file = path, keep.source = FALSE)
cat("ALL_R_FILES_PARSE\t", length(files), "\n", sep = "")
' 2>&1 | tee "${HAC_PREFLIGHT_OUT}/r_parse.log"

"${R_BIN}" --vanilla "${SCRIPT_DIR}/tests/run_tests.R" 2>&1 |
  tee "${HAC_PREFLIGHT_OUT}/contract_tests.log"

echo "HAC_PREFLIGHT_COMPLETE ${HAC_PREFLIGHT_OUT}"
