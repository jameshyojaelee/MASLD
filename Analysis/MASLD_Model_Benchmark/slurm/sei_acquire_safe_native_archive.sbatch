#!/usr/bin/env bash
#SBATCH --job-name=masld-sei-archive
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/alphagenome-sei-admission-sources-21064512"
CONFIG="${ROOT}/config/alphagenome_sei_admission.json"
TARGET="${ROOT}/executions/sei-safe-native-archive-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
ARCHIVE="${STAGE}/checkpoint/sei_model.tar.gz"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Sei archive artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/checkpoint" "${STAGE}/inventory" "${STAGE}/validation"
printf '%s  %s\n' 84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785 \
  "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

source = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(source)
receipt = json.loads((source / "review" / "source_admission_receipt.json").read_text())
sei = receipt["model_dispositions"]["sei"]
if (
    receipt["status"] != "pass"
    or receipt["checkpoint_bytes_downloaded"]
    or not sei["safe_archive_acquisition_allowed"]
    or sei["checkpoint_deserialization_allowed"]
    or sei["gpu_forward_allowed"]
):
    raise SystemExit("Sei source admission does not permit safe archive acquisition")
print("Sei source/license/exposure gate passed")
PY

curl \
  --fail \
  --location \
  --retry 8 \
  --retry-all-errors \
  --connect-timeout 30 \
  --continue-at - \
  --output "${ARCHIVE}" \
  "https://zenodo.org/api/records/4906997/files/sei_model.tar.gz/content"

python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${ARCHIVE}" \
  --output "${STAGE}/inventory/safe_archive_inventory.json" \
  --expected-size 3281639255 \
  --expected-md5 4297aafb711aec4ecccb645b8928ea26 \
  --maximum-members 10000 \
  --maximum-uncompressed-bytes 17179869184 \
  > "${STAGE}/validation/inventory.stdout.txt" \
  2> "${STAGE}/validation/inventory.stderr.txt"

python - "${STAGE}/inventory/safe_archive_inventory.json" <<'PY'
from pathlib import Path
import json
import sys

inventory = json.loads(Path(sys.argv[1]).read_text())
members = {row["path"]: row for row in inventory["members"]}
expected = {
    "model/sei.pth",
    "model/projvec_targets.npy",
    "model/histone_inds.npy",
    "model/seqclass.names",
    "model/target.names",
}
if (
    inventory["status"] != "pass"
    or inventory["archive_extracted"]
    or not inventory["safe_inventory_passed"]
    or inventory["archive_size_bytes"] != 3281639255
    or inventory["archive_md5"] != "4297aafb711aec4ecccb645b8928ea26"
    or not expected.issubset(members)
    or any(members[name]["type"] != "regular_file" for name in expected)
):
    raise SystemExit("Sei safe archive or native member contract differs")
summary = {
    "schema_version": "masld-bench-sei-safe-archive-admission-v1",
    "status": "pass",
    "archive_sha256": inventory["archive_sha256"],
    "archive_md5": inventory["archive_md5"],
    "archive_size_bytes": inventory["archive_size_bytes"],
    "member_count": inventory["member_count"],
    "regular_file_count": inventory["regular_file_count"],
    "uncompressed_regular_bytes": inventory["uncompressed_regular_bytes"],
    "expected_native_members": {
        name: {
            "size_bytes": members[name]["size_bytes"],
            "sha256": members[name]["sha256"],
        }
        for name in sorted(expected)
    },
    "archive_extracted": False,
    "checkpoint_deserialized": False,
    "pickle_loaded": False,
    "model_forward_executed": False,
    "observed_outcomes_loaded": False,
    "next_allowed_action": "isolated_no_network_tensor_only_checkpoint_inspection",
}
(Path(sys.argv[1]).parent / "admission_summary.json").write_text(
    json.dumps(summary, sort_keys=True, indent=2) + "\n"
)
print(json.dumps(summary, sort_keys=True))
PY

curl --version > "${STAGE}/validation/curl.version.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  "${ROOT}/slurm/sei_acquire_safe_native_archive.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "inventory" / "admission_summary.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "sei_safe_native_archive",
        "model_id": "sei",
        "checkpoint_sha256": summary["archive_sha256"],
        "checkpoint_deserialized": False,
        "archive_extracted": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Sei safe native archive artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
