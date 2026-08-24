#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-port-admit
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/enformer-crested-restricted-port-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
CRESTED_REV=934eea18c97992cd098d547cc09edb20b79c28fd
ARCHIVE_URL=https://resources.aertslab.org/CREsted/models/enformer_human.tar.gz
ARCHIVE_SIZE=915174055
ARCHIVE_SHA256=628f67f540304d4d0e143176dc824ed72b3413f78f2fa2efe5d4f0ab51ea1bcc

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Enformer port admission: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/inventories" \
  "${STAGE}/metadata" "${STAGE}/validation"

curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)
curl "${curl_flags[@]}" --output "${STAGE}/downloads/enformer_human.tar.gz" \
  "${ARCHIVE_URL}"
curl "${curl_flags[@]}" --output "${STAGE}/downloads/crested-source.tar.gz" \
  "https://codeload.github.com/aertslab/CREsted/tar.gz/${CRESTED_REV}"

if [[ $(stat -c '%s' "${STAGE}/downloads/enformer_human.tar.gz") != "${ARCHIVE_SIZE}" ]]; then
  printf 'CREsted Enformer archive size differs\n' >&2
  exit 1
fi
printf '%s  %s\n' "${ARCHIVE_SHA256}" \
  "${STAGE}/downloads/enformer_human.tar.gz" | sha256sum --check --strict
gzip -t "${STAGE}/downloads/enformer_human.tar.gz"
gzip -t "${STAGE}/downloads/crested-source.tar.gz"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/enformer_human.tar.gz" \
  --expected-size "${ARCHIVE_SIZE}" \
  --maximum-members 1000 \
  --maximum-uncompressed-bytes 4294967296 \
  --output "${STAGE}/inventories/enformer_human.inventory.json" \
  > "${STAGE}/validation/enformer_human.inventory.stdout.txt"
python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/crested-source.tar.gz" \
  --maximum-members 10000 \
  --maximum-uncompressed-bytes 1073741824 \
  --output "${STAGE}/inventories/crested_source.inventory.json" \
  > "${STAGE}/validation/crested_source.inventory.stdout.txt"

python - "${ROOT}" "${STAGE}" "${CRESTED_REV}" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys

project = Path(sys.argv[1]).resolve(strict=True)
root = Path(sys.argv[2]).resolve(strict=True)
revision = sys.argv[3]
model = json.loads((root / "inventories/enformer_human.inventory.json").read_text())
source = json.loads((root / "inventories/crested_source.inventory.json").read_text())
source_files = {
    row["path"]: row
    for row in source["members"]
    if row["type"] == "regular_file"
}
prefix = f"CREsted-{revision}/"
expected = {
    "LICENSE": "c29d9c8a8b9e76a4a2b0a42f508b71a964f8d1244e74b29d8bb6f552c9ffa496",
    "src/crested/tl/zoo/_enformer.py": "c7a2e74a91ca3f418b986df6b7fb301fb321ed3953248d08ee695bc2e224d474",
    "src/crested/_datasets.py": "e645fb812478838c8963217dcaa75ee6f705fe2e03b7d14f33da909e76adb63e",
}
for relative, digest in expected.items():
    if source_files.get(prefix + relative, {}).get("sha256") != digest:
        raise SystemExit(f"pinned CREsted source differs: {relative}")

model_files = [row for row in model["members"] if row["type"] == "regular_file"]
suffixes = sorted({Path(str(row["path"])).suffix.lower() for row in model_files})
authority_paths = [
    project / "config/artifacts/models/enformer/checkpoints.json",
    project / "config/artifacts/models/enformer/exposure_audit.json",
    project / "config/artifacts/models/enformer/development_crosswalk.json",
    project / "config/models/regulatory_sequence.toml",
]
authority = {
    path.relative_to(project).as_posix(): sha256(path.read_bytes()).hexdigest()
    for path in authority_paths
}
receipt = {
    "schema_version": "masld-bench-enformer-crested-port-admission-v1",
    "status": "pass",
    "model_id": "enformer_crested_restricted_port",
    "registered_scientific_identity": "restricted_conversion_not_native_sonnet",
    "archive_url": "https://resources.aertslab.org/CREsted/models/enformer_human.tar.gz",
    "archive_sha256": model["archive_sha256"],
    "archive_size_bytes": model["archive_size_bytes"],
    "archive_extracted": False,
    "safe_member_inventory_passed": model["safe_inventory_passed"],
    "member_count": model["member_count"],
    "regular_file_count": model["regular_file_count"],
    "uncompressed_regular_bytes": model["uncompressed_regular_bytes"],
    "regular_file_suffixes": suffixes,
    "crested_revision": revision,
    "crested_source_sha256": source["archive_sha256"],
    "authority_sha256": authority,
    "checkpoint_deserialized": False,
    "model_forward_executed": False,
    "project_data_read": False,
    "outcomes_read": False,
    "sealed_features_or_labels_read": False,
    "native_sonnet_identity": False,
    "native_parity_established": False,
    "terms": "internal_academic_noncommercial_nontransferable_no_modification",
    "open_champion_eligible": False,
    "allowed_next_action": "isolated_read_only_strict_runtime_probe_of_this_restricted_port",
    "forbidden_claims": [
        "official_native_Enformer_parity",
        "donor_or_MASLD_context_conditioning",
        "standalone_signed_eQTL_or_ieQTL_prediction",
        "open_champion_eligibility",
    ],
    "terminal_disposition": "restricted_port_bytes_admitted_safe_runtime_and_native_parity_pending",
}
(root / "enformer_crested_port_admission_receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n"
)
PY

python --version > "${STAGE}/validation/python.version.txt"
pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  "${ROOT}/slurm/enformer_crested_acquire_restricted_port.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "enformer_crested_port_admission_receipt.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "restricted_sequence_model_checkpoint_admission",
        "model_id": receipt["model_id"],
        "archive_sha256": receipt["archive_sha256"],
        "native_sonnet_identity": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Enformer CREsted restricted-port admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
