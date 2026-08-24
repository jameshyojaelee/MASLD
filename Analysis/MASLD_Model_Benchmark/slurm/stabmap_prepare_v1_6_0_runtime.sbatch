#!/usr/bin/env bash
#SBATCH --job-name=masld-stabmap160-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
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
SOURCE_URL=https://bioconductor.org/packages/3.23/bioc/src/contrib/StabMap_1.6.0.tar.gz
SOURCE_MD5=bf6e6db94749669866f0aca6f6e8085d
TARGET="${ROOT}/executions/environments/stabmap-v1.6.0-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed StabMap runtime: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite StabMap runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/library" "${STAGE}/source" "${STAGE}/validation"
curl --fail --location --proto '=https' --tlsv1.2 \
  --output "${STAGE}/source/StabMap_1.6.0.tar.gz" "${SOURCE_URL}"
printf '%s  %s\n' "${SOURCE_MD5}" "${STAGE}/source/StabMap_1.6.0.tar.gz" \
  | md5sum --check --strict
python3 - "${STAGE}/source/StabMap_1.6.0.tar.gz" <<'PY'
from pathlib import Path, PurePosixPath
import sys
import tarfile
path = Path(sys.argv[1]).resolve(strict=True)
with tarfile.open(path, "r:gz") as archive:
    members = archive.getmembers()
    if not members:
        raise SystemExit("empty StabMap source archive")
    for member in members:
        value = PurePosixPath(member.name)
        if value.is_absolute() or ".." in value.parts or value.parts[0] != "StabMap":
            raise SystemExit(f"unsafe StabMap member: {member.name}")
        if member.issym() or member.islnk() or member.isdev() or member.isfifo():
            raise SystemExit(f"unsupported StabMap member type: {member.name}")
PY
module load R/4.5.2-mba
Rscript --vanilla - "${STAGE}/library" <<'RS' \
  > "${STAGE}/validation/dependencies.stdout.txt" \
  2> "${STAGE}/validation/dependencies.stderr.txt"
args <- commandArgs(trailingOnly = TRUE)
library_root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(library_root, .libPaths()))
install.packages(
  c("slam", "BiocNeighbors", "BiocSingular"),
  lib = library_root,
  repos = c(
    BioCsoft = "https://bioconductor.org/packages/3.23/bioc",
    CRAN = "https://cran.r-project.org"
  ),
  dependencies = c("Depends", "Imports", "LinkingTo"),
  type = "source",
  Ncpus = 4L
)
required <- c("slam", "BiocNeighbors", "BiocSingular")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("StabMap dependency installation differs")
}
RS
R CMD INSTALL --no-multiarch --library="${STAGE}/library" \
  "${STAGE}/source/StabMap_1.6.0.tar.gz" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"
R_LIBS="${STAGE}/library" Rscript --vanilla - "${STAGE}" <<'RS' \
  > "${STAGE}/validation/runtime.tsv"
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
if (as.character(packageVersion("StabMap")) != "1.6.0") {
  stop("StabMap version differs")
}
required_functions <- c("stabMap", "mosaicDataTopology", "imputeEmbedding")
if (!all(vapply(required_functions, exists, logical(1),
                envir = asNamespace("StabMap"), inherits = FALSE))) {
  stop("StabMap function census differs")
}
cat("status\tpass\n")
cat("R\t", R.version.string, "\n", sep = "")
cat("StabMap\t", as.character(packageVersion("StabMap")), "\n", sep = "")
cat("stabMap\tTRUE\n")
cat("project_data_read\tFALSE\n")
cat("outcomes_read\tFALSE\n")
cat("model_fit_executed\tFALSE\n")
writeLines(capture.output(sessionInfo()), file.path(root, "validation", "sessionInfo.txt"))
RS
sha256sum "${STAGE}/source/StabMap_1.6.0.tar.gz" \
  > "${STAGE}/validation/source.sha256"
sha256sum "${ROOT}/slurm/stabmap_prepare_v1_6_0_runtime.sbatch" \
  > "${STAGE}/source.sha256"
if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'StabMap portable runtime contains a symlink\n' >&2
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
    or records.get("StabMap") != "1.6.0"
    or records.get("stabMap") != "TRUE"
    or records.get("project_data_read") != "FALSE"
    or records.get("outcomes_read") != "FALSE"
    or records.get("model_fit_executed") != "FALSE"
):
    raise SystemExit("StabMap runtime contract differs")
source_sha256 = (root / "validation/source.sha256").read_text().split()[0]
freeze_tree(root, {
    "artifact_class": "stabmap_v1_6_0_portable_runtime",
    "source_url": sys.argv[2],
    "source_md5": sys.argv[3],
    "source_sha256": source_sha256,
    "stabmap_version": "1.6.0",
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'StabMap v1.6.0 portable runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
