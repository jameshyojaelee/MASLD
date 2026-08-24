#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-acquire
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
REVISION=aed3a6f713091fd988196b297e9c06e092ff1d22
SOURCE_URL="https://codeload.github.com/calico/scBasset/tar.gz/${REVISION}"
TARGET="${ROOT}/executions/upstream-acquisitions/scbasset-${REVISION}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scBasset acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/source" "${STAGE}/validation"

curl --fail --location --silent --show-error --retry 5 --retry-delay 5 \
  "${SOURCE_URL}" --output "${STAGE}/downloads/scBasset-${REVISION}.tar.gz"

module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}/downloads/scBasset-${REVISION}.tar.gz" "${STAGE}/source" \
  "${STAGE}/validation/source_records.json" <<'PY'
from __future__ import annotations

from pathlib import Path, PurePosixPath
import hashlib
import json
import sys
import tarfile

archive = Path(sys.argv[1]).resolve(strict=True)
destination = Path(sys.argv[2]).resolve(strict=True)
records_path = Path(sys.argv[3])

with tarfile.open(archive, mode="r:gz") as handle:
    members = handle.getmembers()
    roots: set[str] = set()
    for member in members:
        path = PurePosixPath(member.name)
        if (
            path.is_absolute()
            or ".." in path.parts
            or not (member.isfile() or member.isdir())
        ):
            raise SystemExit(f"unsafe source member: {member.name!r}")
        if path.parts:
            roots.add(path.parts[0])
    if len(roots) != 1:
        raise SystemExit(f"source archive root roster differs: {sorted(roots)!r}")
    root_name = next(iter(roots))
    if not root_name.endswith("aed3a6f713091fd988196b297e9c06e092ff1d22"):
        raise SystemExit(f"source archive revision differs: {root_name!r}")
    handle.extractall(destination, filter="data")

source = destination / root_name
expected = {
    "LICENSE": "4534e268bc50d35b1272461ba82979cac2cf90c0ac3cc8bf97b91fcba1b40747",
    "README.md": "a5a53c6501c3bfa50927621aa0237d8575be0dc829846ed0a7ef1c461f521b87",
    "setup.py": "b7b54078e11ad8bc594512cc33ee38db80c190b06bfde66e3420cb2248ab19d1",
    "requirements-with-version.txt": "48b2e8adebb2320311111a114a25239783884ef746266b09498fe2d5908c1aeb",
    "scbasset/utils.py": "7b653646fc83867cf38c51ae14fdcdf523f0dc99685d1b9446ea5940de71543c",
    "scbasset/basenji_utils.py": "165fdc3c2dd85afd300bf71e27199ae907b005e325c33fa0cc1e9b96be9a92b2",
    "bin/scbasset_preprocess.py": "9d3c47a62182395b1344497140651569691c63f22b6cc7acaf6ed6a7efd347f7",
    "bin/scbasset_train.py": "5773d15e9fdc2a9f679ede6c7a14d6354f067ea92691fd3b7f99a37ca389a155",
}
records = []
for relative, expected_sha in expected.items():
    path = source / relative
    observed_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed_sha != expected_sha:
        raise SystemExit(
            f"registered source hash differs for {relative}: {observed_sha}"
        )
    records.append(
        {
            "path": relative,
            "sha256": observed_sha,
            "size_bytes": path.stat().st_size,
        }
    )
records_path.write_text(
    json.dumps(records, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
(records_path.parent / "source-root.txt").write_text(
    root_name + "\n", encoding="utf-8"
)
PY

sha256sum "${STAGE}/downloads/scBasset-${REVISION}.tar.gz" \
  > "${STAGE}/validation/download.sha256"
curl --version > "${STAGE}/validation/curl.version.txt"
tar --version > "${STAGE}/validation/tar.version.txt"
cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
model_id	scbasset
repository	calico/scBasset
repository_revision	${REVISION}
source_url	${SOURCE_URL}
code_license	Apache-2.0
upstream_checkpoint	none_project_training_only
public_tutorial_weights_downloaded	false
EOF

PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${REVISION}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
revision = sys.argv[2]
archive = root / "downloads" / f"scBasset-{revision}.tar.gz"
freeze_tree(
    root,
    {
        "artifact_class": "scbasset_source_acquisition",
        "model_id": "scbasset",
        "repository_revision": revision,
        "source_archive_sha256": sha256_file(archive),
        "downloads_performed": True,
        "installed": False,
        "public_tutorial_weights_downloaded": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset source acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
