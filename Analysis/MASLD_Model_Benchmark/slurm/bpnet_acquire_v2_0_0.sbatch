#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet-acquire
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=executions/masld-bpnet-acquire-%j.out
#SBATCH --error=executions/masld-bpnet-acquire-%j.err

set -euo pipefail
umask 027

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
REVISION=f4593eedca51741f25b8c5cc0c0d647faf8c8a3c
VERSION=2.0.0
EXPECTED_WHEEL_SHA256=dccff1bfe7c1d726967e95e270b712a9eacd580756d80f6c19ae55ea1be43633
EXPECTED_WHEEL_BYTES=92880
TARGET="${ROOT}/executions/upstream-acquisitions/bpnet-v2.0.0-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite BPNet acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi

mkdir -p "${STAGE}/downloads" "${STAGE}/source" "${STAGE}/validation"
cd "${ROOT}"

curl --fail --location --retry 3 --silent --show-error \
  "https://pypi.org/pypi/bpnet/${VERSION}/json" \
  --output "${STAGE}/downloads/pypi-${VERSION}.json"

python3 - "${STAGE}/downloads/pypi-${VERSION}.json" \
  "${STAGE}/validation/pypi-wheel-url.txt" \
  "${EXPECTED_WHEEL_SHA256}" "${EXPECTED_WHEEL_BYTES}" <<'PY'
import json
import pathlib
import sys

record_path = pathlib.Path(sys.argv[1])
url_path = pathlib.Path(sys.argv[2])
expected_hash = sys.argv[3]
expected_size = int(sys.argv[4])
record = json.loads(record_path.read_text(encoding="utf-8"))
if record["info"]["version"] != "2.0.0":
    raise SystemExit("PyPI release version differs")
candidates = [
    item for item in record["urls"]
    if item["filename"] == "bpnet-2.0.0-py3-none-any.whl"
]
if len(candidates) != 1:
    raise SystemExit("exact BPNet wheel is not unique")
wheel = candidates[0]
if wheel["digests"]["sha256"] != expected_hash:
    raise SystemExit("PyPI BPNet wheel digest differs")
if wheel["size"] != expected_size:
    raise SystemExit("PyPI BPNet wheel size differs")
url_path.write_text(wheel["url"] + "\n", encoding="utf-8")
PY

WHEEL_URL="$(tr -d '\n' < "${STAGE}/validation/pypi-wheel-url.txt")"
curl --fail --location --retry 3 --silent --show-error \
  "${WHEEL_URL}" --output "${STAGE}/downloads/bpnet-2.0.0-py3-none-any.whl"
test "$(stat -c '%s' "${STAGE}/downloads/bpnet-2.0.0-py3-none-any.whl")" \
  -eq "${EXPECTED_WHEEL_BYTES}"
printf '%s  %s\n' "${EXPECTED_WHEEL_SHA256}" \
  "${STAGE}/downloads/bpnet-2.0.0-py3-none-any.whl" | sha256sum -c -

SOURCE_ARCHIVE="${STAGE}/downloads/bpnet-${REVISION}.tar.gz"
curl --fail --location --retry 3 --silent --show-error \
  "https://codeload.github.com/kundajelab/bpnet/tar.gz/${REVISION}" \
  --output "${SOURCE_ARCHIVE}"

python3 - "${SOURCE_ARCHIVE}" "${STAGE}/source" "${REVISION}" <<'PY'
import pathlib
import tarfile
import sys

archive = pathlib.Path(sys.argv[1])
destination = pathlib.Path(sys.argv[2])
revision = sys.argv[3]
expected_root = f"bpnet-{revision}"
with tarfile.open(archive, "r:gz") as handle:
    members = handle.getmembers()
    if not members:
        raise SystemExit("BPNet source archive is empty")
    for member in members:
        path = pathlib.PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"unsafe source member: {member.name}")
        if not path.parts or path.parts[0] != expected_root:
            raise SystemExit(f"unexpected archive root: {member.name}")
        if member.issym() or member.islnk() or member.isdev():
            raise SystemExit(f"link/device source member rejected: {member.name}")
    handle.extractall(destination, members=members)
PY

SOURCE="${STAGE}/source/bpnet-${REVISION}"
python3 - "${SOURCE}" <<'PY'
import hashlib
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
expected = {
    "LICENSE": (1068, "07cc7fd1890c0f369dfc2d35bfd5fc5bc8bde189401c35ae0f4e4ef63f088a96"),
    "README.md": (16499, "f1a54bf17f0dac0e17ae6986d35a0130bf7e9fd04d4f946f73c3314ce30b0e7d"),
    "setup.py": (2011, "c9f5c42f741e688805e8dcc487214cb28aef40bd92bb34b2816f22f344aad8ab"),
    "bpnet/model/arch.py": (28883, "cbc05bd2e84c5e68b30785ab47ccce207ea4218e414ecce5bf995f523cd16b86"),
    "bpnet/model/bpnetdefaults.py": (2071, "3cff612b084420b4c45ab70073b429638acd4958a3b2672a561a5194c15b7a88"),
    "bpnet/model/losses.py": (1378, "fe5f0481ebb6cf33833f7cc029e31719806fad0345ac19309faedbf7cab2a023"),
    "bpnet/train/training.py": (29668, "bd7aa3eb4975cc91cd1a9a4c59e97c5962c131a2c24afb29090669eec1962092"),
    "bpnet/generators/generators.py": (49861, "34ceebd89052cba57c836001cebcfda251fa6ffdc7345d839142023c3b96c3b2"),
}
for relative, (size, digest) in expected.items():
    path = root / relative
    if not path.is_file() or path.is_symlink():
        raise SystemExit(f"missing or unsafe source record: {relative}")
    if path.stat().st_size != size:
        raise SystemExit(f"source size differs: {relative}")
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != digest:
        raise SystemExit(f"source digest differs: {relative}")
setup_text = (root / "setup.py").read_text(encoding="utf-8")
if "version='2.0.0'" not in setup_text:
    raise SystemExit("source package version differs")
PY

python3 - "${STAGE}/downloads/bpnet-2.0.0-py3-none-any.whl" \
  "${STAGE}/validation/wheel-members.json" <<'PY'
import json
import pathlib
import sys
import zipfile

wheel = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
with zipfile.ZipFile(wheel) as handle:
    names = sorted(handle.namelist())
    for name in names:
        path = pathlib.PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"unsafe wheel member: {name}")
required = {
    "bpnet/model/arch.py",
    "bpnet/model/bpnetdefaults.py",
    "bpnet/model/custommodel.py",
    "bpnet/model/losses.py",
}
if not required.issubset(names):
    raise SystemExit("BPNet wheel member roster differs")
output.write_text(json.dumps({"members": names}, indent=2, sort_keys=True) + "\n")
PY

sha256sum "${SOURCE_ARCHIVE}" > "${STAGE}/validation/source-archive.sha256"
sha256sum "${STAGE}/downloads/bpnet-2.0.0-py3-none-any.whl" \
  > "${STAGE}/validation/wheel.sha256"
sha256sum "${ROOT}/config/artifacts/models/bpnet/checkpoints.json" \
  "${ROOT}/config/artifacts/models/bpnet/development_crosswalk.json" \
  "${ROOT}/config/artifacts/models/bpnet/exposure_audit.json" \
  "${ROOT}/slurm/bpnet_acquire_v2_0_0.sbatch" \
  > "${STAGE}/validation/input-locks.sha256"

cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
model_id	bpnet
repository	kundajelab/bpnet
revision	${REVISION}
pypi_version	${VERSION}
wheel_sha256	${EXPECTED_WHEEL_SHA256}
upstream_checkpoint	none_project_trained
EOF

python3 - "${STAGE}" "${REVISION}" "${EXPECTED_WHEEL_SHA256}" <<'PY'
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
revision = sys.argv[2]
wheel_hash = sys.argv[3]
items = []
for path in sorted(root.rglob("*")):
    if path.is_file() and path.name not in {"ARTIFACTS.json", "COMPLETE"}:
        items.append({
            "path": path.relative_to(root).as_posix(),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
        })
payload = {
    "schema_version": "masld-bench-artifacts-v1",
    "artifacts": items,
    "metadata": {
        "artifact_class": "bpnet_source_acquisition",
        "model_id": "bpnet",
        "repository": "kundajelab/bpnet",
        "revision": revision,
        "pypi_version": "2.0.0",
        "wheel_sha256": wheel_hash,
        "upstream_checkpoint": None,
    },
}
(root / "ARTIFACTS.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

ARTIFACTS_SHA256="$(sha256sum "${STAGE}/ARTIFACTS.json" | cut -d' ' -f1)"
cat > "${STAGE}/COMPLETE" <<EOF
status	complete
artifact_class	bpnet_source_acquisition
artifacts_sha256	${ARTIFACTS_SHA256}
EOF

mv "${STAGE}" "${TARGET}"
printf 'BPNet source acquisition: %s\n' "${TARGET}"
