#!/usr/bin/env bash
#SBATCH --job-name=masld-seurat551-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
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
SOURCE_URL=https://cran.r-project.org/src/contrib/Seurat_5.5.1.tar.gz
SOURCE_MD5=b07be078da2106bfa6cd0572744fdbb6
TARGET="${ROOT}/executions/environments/seurat-5.5.1-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed Seurat 5.5.1 stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed Seurat 5.5.1 stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite Seurat 5.5.1 runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/library" "${STAGE}/source" "${STAGE}/validation"
curl --fail --location --proto '=https' --tlsv1.2 \
  --output "${STAGE}/source/Seurat_5.5.1.tar.gz" "${SOURCE_URL}"
printf '%s  %s\n' "${SOURCE_MD5}" "${STAGE}/source/Seurat_5.5.1.tar.gz" \
  | md5sum --check --strict
python3 - "${STAGE}/source/Seurat_5.5.1.tar.gz" <<'PY'
from pathlib import Path, PurePosixPath
import sys
import tarfile

path = Path(sys.argv[1]).resolve(strict=True)
with tarfile.open(path, "r:gz") as archive:
    members = archive.getmembers()
    if not members:
        raise SystemExit("empty Seurat source archive")
    for member in members:
        value = PurePosixPath(member.name)
        if value.is_absolute() or ".." in value.parts or value.parts[0] != "Seurat":
            raise SystemExit(f"unsafe Seurat source member: {member.name}")
        if member.issym() or member.islnk() or member.isdev() or member.isfifo():
            raise SystemExit(f"unsupported Seurat source member type: {member.name}")
PY

module load R/4.5.2-mba
R CMD INSTALL --no-multiarch --library="${STAGE}/library" \
  "${STAGE}/source/Seurat_5.5.1.tar.gz" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"
R_LIBS="${STAGE}/library" Rscript --vanilla - "${STAGE}" <<'RS' \
  > "${STAGE}/validation/runtime.tsv"
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
if (as.character(packageVersion("Seurat")) != "5.5.1") {
  stop("Seurat version differs")
}
if (!exists("FindMultiModalNeighbors", envir = asNamespace("Seurat"), inherits = FALSE)) {
  stop("FindMultiModalNeighbors is absent")
}
required <- c(
  "Seurat", "SeuratObject", "Matrix", "Rcpp", "RcppEigen", "RcppAnnoy",
  "RcppHNSW", "igraph", "future", "sctransform"
)
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("required Seurat dependency is absent")
}
cat("status\tpass\n")
cat("R\t", R.version.string, "\n", sep = "")
cat("platform\t", R.version$platform, "\n", sep = "")
for (package in required) {
  cat(package, "\t", as.character(packageVersion(package)), "\n", sep = "")
}
cat("FindMultiModalNeighbors\tTRUE\n")
cat("local_library_first\t", identical(normalizePath(.libPaths()[[1]]),
    normalizePath(file.path(root, "library"))), "\n", sep = "")
cat("project_data_read\tFALSE\n")
cat("outcomes_read\tFALSE\n")
cat("model_fit_executed\tFALSE\n")
writeLines(capture.output(sessionInfo()), file.path(root, "validation", "sessionInfo.txt"))
RS
sha256sum "${STAGE}/source/Seurat_5.5.1.tar.gz" \
  > "${STAGE}/validation/source.sha256"
sha256sum "${ROOT}/slurm/seurat_wnn_prepare_v5_5_1_runtime.sbatch" \
  > "${STAGE}/source.sha256"
if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'Seurat 5.5.1 portable runtime contains a symlink\n' >&2
  exit 1
fi

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SOURCE_URL}" "${SOURCE_MD5}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
records = {}
for line in (root / "validation/runtime.tsv").read_text().splitlines():
    key, value = line.split("\t", 1)
    records[key] = value
if (
    records.get("status") != "pass"
    or records.get("Seurat") != "5.5.1"
    or records.get("FindMultiModalNeighbors") != "TRUE"
    or records.get("local_library_first") != "TRUE"
    or records.get("project_data_read") != "FALSE"
    or records.get("outcomes_read") != "FALSE"
    or records.get("model_fit_executed") != "FALSE"
):
    raise SystemExit("Seurat 5.5.1 runtime contract differs")
source_sha256 = (root / "validation/source.sha256").read_text().split()[0]
freeze_tree(
    root,
    {
        "artifact_class": "seurat_5_5_1_portable_runtime",
        "source_url": sys.argv[2],
        "source_md5": sys.argv[3],
        "source_sha256": source_sha256,
        "seurat_version": "5.5.1",
        "find_multimodal_neighbors_available": True,
        "project_data_read": False,
        "outcomes_read": False,
        "model_fit_executed": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Seurat 5.5.1 portable runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
