#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-bias-acq
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
TARGET="${ROOT}/executions/chrombpnet-bias-models-zenodo7443683-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
RECORD_URL=https://zenodo.org/api/records/7443683
FILE_URL='https://zenodo.org/records/7443683/files/bias_models.zip?download=1'
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ChromBPNet bias acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir -m 750 "${STAGE}"
module load curl/8.7.1-GCCcore-13.3.0
curl --fail --location --retry 4 --retry-delay 5 \
  "${RECORD_URL}" --output "${STAGE}/zenodo-record-7443683.json"
curl --fail --location --retry 4 --retry-delay 5 \
  "${FILE_URL}" --output "${STAGE}/bias_models.zip"
printf '%s  %s\n' 8205306b23dbf54b2081cba4ad1fabbc \
  "${STAGE}/bias_models.zip" | md5sum --check --strict
curl --version > "${STAGE}/curl.version.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" <<'PY'
from pathlib import Path, PurePosixPath
import json
import sys
import zipfile

root = Path(sys.argv[1]).resolve(strict=True)
metadata = json.loads((root / "zenodo-record-7443683.json").read_text())
if (
    metadata.get("id") != 7443683
    or metadata.get("doi") != "10.5281/zenodo.7443683"
    or metadata.get("metadata", {}).get("version") != "1"
):
    raise SystemExit("Zenodo record metadata differs")
archive = root / "bias_models.zip"
model_root = root / "models"
model_root.mkdir(mode=0o750)
with zipfile.ZipFile(archive) as handle:
    members = handle.infolist()
    if not members:
        raise SystemExit("bias-model archive is empty")
    for member in members:
        path = PurePosixPath(member.filename)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit("bias-model archive contains an unsafe path")
    handle.extractall(model_root)
h5_files = sorted(model_root.rglob("*.h5"))
if not h5_files:
    raise SystemExit("bias-model archive contains no HDF5 weights")
license_value = metadata.get("metadata", {}).get("license")
validation = {
    "record_id": 7443683,
    "doi": metadata["doi"],
    "version": metadata["metadata"]["version"],
    "archive_members": len(members),
    "h5_files": [str(path.relative_to(root)) for path in h5_files],
    "zenodo_license_field": license_value,
    "weight_license_status": "unknown_not_redistributable",
    "allowed_role": "restricted_smoke_comparator",
}
(root / "validation.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY
python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum "${ROOT}/slurm/acquire_chrombpnet_bias_models.sbatch" \
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
            "artifact_class": "chrombpnet_pretrained_bias_models",
            "record_id": validation["record_id"],
            "doi": validation["doi"],
            "h5_files": len(validation["h5_files"]),
            "weight_license_status": validation["weight_license_status"],
            "allowed_role": validation["allowed_role"],
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet bias models: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
