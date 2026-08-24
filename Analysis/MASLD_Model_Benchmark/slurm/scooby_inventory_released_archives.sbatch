#!/usr/bin/env bash
#SBATCH --job-name=masld-scooby-archive-inventory
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
REVIEW="${ROOT}/executions/scooby-cross-review-21064420"
REVIEW_ARTIFACTS_SHA256=27668bff0943a8baac6489f073e000c2ba181e729a6f3360ea1dc0a3cfd6aa7c
TARGET="${ROOT}/executions/scooby-released-archive-inventory-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scooby archive inventory: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ $(sha256sum "${REVIEW}/ARTIFACTS.json" | cut -d ' ' -f 1) != "${REVIEW_ARTIFACTS_SHA256}" ]]; then
  printf 'Pinned scooby cross-review artifact differs\n' >&2
  exit 1
fi

install -d -m 0750 "${STAGE}/inventories" "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/scooby_remote_zip_inventory.py" \
  --url 'https://zenodo.org/records/15517072/files/training_data.zip?download=1' \
  --expected-size 49746574762 \
  --expected-md5 bac5ca14612afad2456de2d711fabff0 \
  --output "${STAGE}/inventories/training_data.inventory.json" \
  > "${STAGE}/validation/training_data.stdout.json"
python "${ROOT}/scripts/scooby_remote_zip_inventory.py" \
  --url 'https://zenodo.org/records/13891693/files/submission_data.zip?download=1' \
  --expected-size 24639962286 \
  --expected-md5 be466ec760807465ec5172d7f88c332c \
  --output "${STAGE}/inventories/submission_data.inventory.json" \
  > "${STAGE}/validation/submission_data.stdout.json"

python - "${STAGE}" "${REVIEW}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
review_root = Path(sys.argv[2]).resolve(strict=True)
review = json.loads((review_root / "scooby_cross_review_receipt.json").read_text())
training = json.loads((root / "inventories/training_data.inventory.json").read_text())
submission = json.loads((root / "inventories/submission_data.inventory.json").read_text())
if (
    review["status"] != "pass"
    or review["checkpoint_forward_allowed"]
    or review["rna_conditioned_atac_models"]
    or review["production_allowed"] != {
        "scooby_onek1k": "metadata_and_bounded_header_only",
        "scooby_epicardioids": "metadata_and_bounded_header_only",
        "scooby_neurips": "metadata_and_pointer_only",
    }
):
    raise SystemExit("scooby cross-review disposition differs")
for inventory in (training, submission):
    if (
        inventory["status"] != "pass"
        or inventory["archive_payload_downloaded"]
        or inventory["archive_extracted"]
        or inventory["declared_archive_md5_verified_from_full_payload"]
    ):
        raise SystemExit("remote ZIP inventory disposition differs")
receipt = {
    "schema_version": "masld-bench-scooby-released-archive-inventory-v1",
    "status": "pass",
    "cross_review_artifact": str(review_root),
    "archives": {
        "training_data": {
            "entry_count": training["entry_count"],
            "context_candidate_count": training["context_candidate_count"],
            "http_bytes_downloaded": training["http_bytes_downloaded"],
        },
        "submission_data": {
            "entry_count": submission["entry_count"],
            "context_candidate_count": submission["context_candidate_count"],
            "http_bytes_downloaded": submission["http_bytes_downloaded"],
        },
    },
    "archive_payload_downloaded": False,
    "archive_extracted": False,
    "checkpoint_bytes_downloaded": False,
    "project_data_read": False,
    "outcomes_read": False,
    "terminal_disposition": "central_directories_inventoried_context_artifact_adjudication_pending",
}
(root / "scooby_released_archive_inventory_receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY

sha256sum \
  "${ROOT}/scripts/scooby_remote_zip_inventory.py" \
  "${ROOT}/slurm/scooby_inventory_released_archives.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (root / "scooby_released_archive_inventory_receipt.json").read_text()
)
freeze_tree(
    root,
    {
        "artifact_class": "scooby_released_archive_inventory",
        "archives": receipt["archives"],
        "archive_payload_downloaded": False,
        "checkpoint_bytes_downloaded": False,
        "project_data_read": False,
        "outcomes_read": False,
        "terminal_disposition": receipt["terminal_disposition"],
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scooby released archive inventory: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
