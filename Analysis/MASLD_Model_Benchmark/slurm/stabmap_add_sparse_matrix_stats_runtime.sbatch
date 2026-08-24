#!/usr/bin/env bash
#SBATCH --job-name=masld-stabmap-sparse-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PARENT="${ROOT}/executions/environments/stabmap-v1.6.0-portable-21066039"
PARENT_ARTIFACTS_SHA256=acbff9a5df6b3d0f770d2ab82b0f7a72c57b1c6beb83a43191e3d8548524fb91
SOURCE_URL=https://bioconductor.org/packages/3.23/bioc/src/contrib/sparseMatrixStats_1.24.0.tar.gz
SOURCE_MD5=e4fa1a5604d206580a5e7fb7ffa753eb
TARGET="${ROOT}/executions/environments/stabmap-v1.6.0-sparse-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed StabMap sparse runtime: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite StabMap sparse runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/source" "${STAGE}/validation"
printf '%s  %s\n' "${PARENT_ARTIFACTS_SHA256}" "${PARENT}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PARENT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
module purge
cp -a "${PARENT}/library" "${STAGE}/library"
curl --fail --location --silent --show-error \
  --output "${STAGE}/source/sparseMatrixStats_1.24.0.tar.gz" "${SOURCE_URL}"
printf '%s  %s\n' "${SOURCE_MD5}" \
  "${STAGE}/source/sparseMatrixStats_1.24.0.tar.gz" | md5sum --check --strict
module load R/4.5.2-mba
R CMD INSTALL --no-multiarch --library="${STAGE}/library" \
  "${STAGE}/source/sparseMatrixStats_1.24.0.tar.gz" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"
R_LIBS="${STAGE}/library" Rscript --vanilla - <<'RS' \
  > "${STAGE}/validation/sparse_method.tsv" \
  2> "${STAGE}/validation/sparse_method.stderr.txt"
suppressPackageStartupMessages({
  library(Matrix)
  library(StabMap)
  library(sparseMatrixStats)
})
if (as.character(packageVersion("StabMap")) != "1.6.0") stop("StabMap version differs")
if (as.character(packageVersion("sparseMatrixStats")) != "1.24.0") {
  stop("sparseMatrixStats version differs")
}
x <- Matrix::Matrix(c(1, 0, 2, 1, 0, 3), nrow = 2, sparse = TRUE)
value <- MatrixGenerics::rowVars(x)
if (length(value) != 2L || any(!is.finite(value))) stop("sparse rowVars dispatch differs")
cat("StabMap\t1.6.0\n")
cat("sparseMatrixStats\t1.24.0\n")
cat("dgCMatrix_rowVars\tTRUE\n")
RS
R_LIBS="${STAGE}/library" Rscript --vanilla - <<'RS' \
  > "${STAGE}/validation/sessionInfo.txt"
suppressPackageStartupMessages(library(sparseMatrixStats))
sessionInfo()
RS
sha256sum "${STAGE}/source/sparseMatrixStats_1.24.0.tar.gz" \
  > "${STAGE}/source.sha256"
if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'StabMap sparse runtime contains a symlink\n' >&2
  exit 1
fi
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
records = {}
for line in (root / "validation/sparse_method.tsv").read_text().splitlines():
    key, value = line.split("\t")
    records[key] = value
if records != {
    "StabMap": "1.6.0",
    "sparseMatrixStats": "1.24.0",
    "dgCMatrix_rowVars": "TRUE",
}:
    raise SystemExit("StabMap sparse runtime contract differs")
freeze_tree(root, {
    "artifact_class": "stabmap_v1_6_0_sparse_portable_runtime",
    "model_id": "stabmap",
    "exact_package": "StabMap 1.6.0",
    "sparse_matrix_stats": "1.24.0",
    "sparse_matrix_stats_source_md5": "e4fa1a5604d206580a5e7fb7ffa753eb",
    "parent_runtime_artifacts_sha256": "acbff9a5df6b3d0f770d2ab82b0f7a72c57b1c6beb83a43191e3d8548524fb91",
    "dgCMatrix_rowVars_dispatch": True,
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'StabMap sparse portable runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
