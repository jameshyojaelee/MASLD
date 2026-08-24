#!/usr/bin/env bash
#SBATCH --job-name=masld-bedtools-runtime
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
TARGET="${ROOT}/executions/environments/bedtools-bioconda-2.31.1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
BUILD_ROOT="${TMPDIR:-/tmp}/masld-bedtools-${SLURM_JOB_ID}"
BUILD_ENV="${BUILD_ROOT}/env"
PACKAGE_CACHE="${BUILD_ROOT}/package-cache"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Bedtools runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ -e "${BUILD_ROOT}" ]]; then
  printf 'Refusing to overwrite Bedtools scratch: %s\n' "${BUILD_ROOT}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/validation" "${STAGE}/runtime" "${STAGE}/packages" \
  "${PACKAGE_CACHE}"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict

module load Mamba/23.11.0-0
export CONDA_PKGS_DIRS="${PACKAGE_CACHE}"
mamba --version > "${STAGE}/validation/mamba.version.txt"
mamba create --yes \
  --prefix "${BUILD_ENV}" \
  --override-channels \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority \
  bedtools=2.31.1 \
  > "${STAGE}/validation/mamba-create.stdout.txt" \
  2> "${STAGE}/validation/mamba-create.stderr.txt"
mamba list --prefix "${BUILD_ENV}" --explicit \
  > "${STAGE}/conda-explicit.txt"
mamba list --prefix "${BUILD_ENV}" --json \
  > "${STAGE}/conda-list.json"
cp -aL "${BUILD_ENV}/." "${STAGE}/runtime/"
find "${PACKAGE_CACHE}" -maxdepth 1 -type f \
  \( -name '*.conda' -o -name '*.tar.bz2' \) \
  -exec cp -p '{}' "${STAGE}/packages/" \;
if find "${STAGE}/runtime" -type l -print -quit | grep -q .; then
  printf 'Materialized Bedtools runtime still contains a symlink\n' >&2
  exit 1
fi
"${STAGE}/runtime/bin/bedtools" --version \
  > "${STAGE}/validation/host-bedtools.version.txt"
ldd "${STAGE}/runtime/bin/bedtools" \
  > "${STAGE}/validation/host-bedtools.ldd.txt"

module purge
module load singularity/4.2.2
export SINGULARITYENV_PREPEND_PATH="${STAGE}/runtime/bin"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  "${IMAGE}" bedtools --version \
  > "${STAGE}/validation/container-bedtools.version.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  "${IMAGE}" ldd "${STAGE}/runtime/bin/bedtools" \
  > "${STAGE}/validation/container-bedtools.ldd.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" <<'PY'
from hashlib import sha256
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
metadata_paths = sorted(
    (root / "runtime" / "conda-meta").glob("bedtools-*.json")
)
if len(metadata_paths) != 1:
    raise SystemExit("Bedtools package metadata count differs")
metadata = json.loads(metadata_paths[0].read_text())
if metadata.get("name") != "bedtools" or metadata.get("version") != "2.31.1":
    raise SystemExit("Bedtools package identity differs")
if not metadata.get("license") or "bioconda" not in str(metadata.get("channel")):
    raise SystemExit("Bedtools package license or channel provenance differs")
versions = {
    name: (root / "validation" / name).read_text().strip()
    for name in ("host-bedtools.version.txt", "container-bedtools.version.txt")
}
if set(versions.values()) != {"bedtools v2.31.1"}:
    raise SystemExit(f"Bedtools version probe differs: {versions!r}")
binary = root / "runtime" / "bin" / "bedtools"
digest = sha256()
with binary.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
validation = {
    "schema_version": "masld-bench-bedtools-runtime-v1",
    "status": "pass",
    "package": "bedtools",
    "version": "2.31.1",
    "build": metadata.get("build"),
    "channel": metadata.get("channel"),
    "license": metadata.get("license"),
    "license_family": metadata.get("license_family"),
    "binary_sha256": digest.hexdigest(),
    "host_version": versions["host-bedtools.version.txt"],
    "container_version": versions["container-bedtools.version.txt"],
    "compatible_with_chrombpnet_tf2_8_2_image": True,
}
(root / "validation.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
sha256sum "${ROOT}/slurm/build_bedtools_bioconda_runtime.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
validation = json.loads((root / "validation.json").read_text())
print(
    freeze_tree(
        root,
        {
            "artifact_class": "bedtools_bioconda_runtime",
            "package": validation["package"],
            "version": validation["version"],
            "build": validation["build"],
            "license": validation["license"],
            "binary_sha256": validation["binary_sha256"],
            "compatible_with_chrombpnet_tf2_8_2_image": True,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Bedtools Bioconda runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
