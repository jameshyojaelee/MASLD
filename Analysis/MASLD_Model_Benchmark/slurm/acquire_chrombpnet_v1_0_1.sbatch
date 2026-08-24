#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-acquire
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
REVISION=eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c
WHEEL_URL=https://files.pythonhosted.org/packages/13/57/18efcef9b0e7613e377fdcd2a8091cc8ba8bc76386e255a9249fa1552062/chrombpnet-1.0.1-py3-none-any.whl
SDIST_URL=https://files.pythonhosted.org/packages/6d/fc/89a993fa1537caca28bfc8b38e01c97b4aa3403f02bac20a96826a2acab1/chrombpnet-1.0.1.tar.gz
PYPI_JSON_URL=https://pypi.org/pypi/chrombpnet/1.0.1/json
SOURCE_URL="https://codeload.github.com/kundajelab/ChromBPNet/tar.gz/${REVISION}"
TARGET="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite acquisition output: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/source" "${STAGE}/validation"

curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
  "${WHEEL_URL}" --output "${STAGE}/downloads/chrombpnet-1.0.1-py3-none-any.whl"
curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
  "${SDIST_URL}" --output "${STAGE}/downloads/chrombpnet-1.0.1.tar.gz"
curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
  "${PYPI_JSON_URL}" --output "${STAGE}/downloads/pypi-v1.0.1.json"
curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
  "${SOURCE_URL}" --output "${STAGE}/downloads/ChromBPNet-${REVISION}.tar.gz"

test "$(stat -c '%s' "${STAGE}/downloads/chrombpnet-1.0.1-py3-none-any.whl")" -eq 440753
test "$(stat -c '%s' "${STAGE}/downloads/chrombpnet-1.0.1.tar.gz")" -eq 410041
printf '%s  %s\n' \
  fc2986a04d9cb0c0e4737c823b056ddd97edc77fac9b495fe0db9963531aad59 \
  "${STAGE}/downloads/chrombpnet-1.0.1-py3-none-any.whl" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  640bbf75e34fb2428ddbdb4de6ac646d2353ebe8885f4b02f86c9918294f1059 \
  "${STAGE}/downloads/chrombpnet-1.0.1.tar.gz" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
python - \
  "${STAGE}/downloads/pypi-v1.0.1.json" \
  "${STAGE}/downloads/chrombpnet-1.0.1-py3-none-any.whl" \
  "${STAGE}/downloads/ChromBPNet-${REVISION}.tar.gz" <<'PY'
from pathlib import Path, PurePosixPath
import hashlib
import json
import sys
import tarfile
import zipfile

pypi_path = Path(sys.argv[1]).resolve(strict=True)
wheel_path = Path(sys.argv[2]).resolve(strict=True)
source_archive = Path(sys.argv[3]).resolve(strict=True)
pypi = json.loads(pypi_path.read_text(encoding="utf-8"))
if pypi.get("info", {}).get("version") != "1.0.1":
    raise SystemExit("PyPI release version differs")
expected = {
    "chrombpnet-1.0.1-py3-none-any.whl": (
        440753,
        "fc2986a04d9cb0c0e4737c823b056ddd97edc77fac9b495fe0db9963531aad59",
    ),
    "chrombpnet-1.0.1.tar.gz": (
        410041,
        "640bbf75e34fb2428ddbdb4de6ac646d2353ebe8885f4b02f86c9918294f1059",
    ),
}
observed = {
    item["filename"]: (item["size"], item["digests"]["sha256"])
    for item in pypi.get("urls", [])
}
if observed != expected:
    raise SystemExit("PyPI file roster differs")
requirements = set(pypi.get("info", {}).get("requires_dist") or [])
for requirement in {
    "tensorflow==2.8.0",
    "tensorflow-estimator==2.8.0",
    "tensorflow-probability==0.15.0",
    "numpy==1.23.4",
    "protobuf==3.20",
}:
    if requirement not in requirements:
        raise SystemExit(f"required native dependency changed: {requirement}")

with zipfile.ZipFile(wheel_path) as handle:
    names = handle.namelist()
    for name in names:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"unsafe wheel member: {name!r}")
    if len(names) != len(set(names)):
        raise SystemExit("wheel contains duplicate members")

with tarfile.open(source_archive, mode="r:gz") as handle:
    members = handle.getmembers()
root_names = set()
for member in members:
    path = PurePosixPath(member.name)
    if path.is_absolute() or ".." in path.parts or not (member.isfile() or member.isdir()):
        raise SystemExit(f"unsafe source member: {member.name!r}")
    if path.parts:
        root_names.add(path.parts[0])
if len(root_names) != 1:
    raise SystemExit("source archive does not have one root directory")
Path(sys.argv[3]).with_suffix(Path(sys.argv[3]).suffix + ".root.txt").write_text(
    next(iter(root_names)) + "\n", encoding="utf-8"
)
PY

tar --extract \
  --file "${STAGE}/downloads/ChromBPNet-${REVISION}.tar.gz" \
  --directory "${STAGE}/source" \
  --no-same-owner \
  --no-same-permissions

python - "${STAGE}/source" "${STAGE}/validation/source_records.json" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2])
expected = {
    "LICENSE": "eee7b4d55be619630ce91024410485d25226c387e72315b83712a8dbf89188ef",
    "README.md": "16702b378ea6a158614be001b21fe8a420a5603141fec063edbf9c14d9ca6add",
    "requirements.txt": "a2299b7298690f1bdbccf59172d512ee50edf94a042ba8088553864d05cde358",
    "setup.py": "7dde221fd752bcc421da88d6d151c16c77c4254c92a19d7b6324f13fc93f5e9c",
    "bpnet_model.py": "8ddbab40fa13d4ee7b5e326560e184a2502ddbd009444c5143420bca12d7ee28",
    "chrombpnet_with_bias_model.py": "bedc63a36ee27bb25ab39a00f1d9e3d4d64809754bf46aebe10d17b3a331dd33",
    "train.py": "3a569e319fae561de89d99eab95a15bd369a5ff00922c33bb60342bd6bb9a77f",
    "augment.py": "8b560dd4d9563f2620c2b6c1043f8a8220385c3c22d9e84fc289e1e80cfbdea7",
    "losses.py": "c32e816398274bf7995cf83aab9bebb3b40ddaac7d162cebfa0672c31d259a6d",
}
records = []
for basename, digest in expected.items():
    matches = []
    for path in root.rglob(basename):
        if path.is_file() and not path.is_symlink():
            observed = hashlib.sha256(path.read_bytes()).hexdigest()
            if observed == digest:
                matches.append(path)
    if len(matches) != 1:
        raise SystemExit(
            f"expected exactly one source match for {basename}; observed {len(matches)}"
        )
    path = matches[0]
    records.append(
        {
            "basename": basename,
            "path": path.relative_to(root).as_posix(),
            "sha256": digest,
            "size_bytes": path.stat().st_size,
        }
    )
output.write_text(
    json.dumps(records, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY

sha256sum "${STAGE}"/downloads/* > "${STAGE}/validation/downloads.sha256"
curl --version > "${STAGE}/validation/curl.version.txt"
tar --version > "${STAGE}/validation/tar.version.txt"
cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
model_id	chrombpnet
release	1.0.1
repository	kundajelab/ChromBPNet
repository_revision	${REVISION}
source_url	${SOURCE_URL}
wheel_url	${WHEEL_URL}
wheel_sha256	fc2986a04d9cb0c0e4737c823b056ddd97edc77fac9b495fe0db9963531aad59
sdist_url	${SDIST_URL}
sdist_sha256	640bbf75e34fb2428ddbdb4de6ac646d2353ebe8885f4b02f86c9918294f1059
code_license	MIT
upstream_checkpoint	none_project_training_only
EOF

PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${REVISION}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "chrombpnet_source_acquisition",
        "model_id": "chrombpnet",
        "release": "1.0.1",
        "repository_revision": sys.argv[2],
        "wheel_sha256": sha256_file(
            root / "downloads" / "chrombpnet-1.0.1-py3-none-any.whl"
        ),
        "sdist_sha256": sha256_file(
            root / "downloads" / "chrombpnet-1.0.1.tar.gz"
        ),
        "downloads_performed": True,
        "installed": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'ChromBPNet source acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
