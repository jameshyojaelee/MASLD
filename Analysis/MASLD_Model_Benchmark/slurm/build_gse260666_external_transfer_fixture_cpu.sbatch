#!/usr/bin/env bash
#SBATCH --job-name=model-data-068
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ACTIVATION=${ROOT}/executions/gse260666-activation-21066072
MATRIX=${ROOT}/executions/bulk-expansion-matrix-audit-21065853
TRAINING_MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
TRAINING_JOIN=${ROOT}/executions/gse267145-authoritative-join-21064930
TASK=${ROOT}/config/evaluation/gse260666_external_histology_transfer_task.toml
PRIVATE_SALT_ROOT=${ROOT}/executions/private/fingerprint-salts
PRIVATE_SALT=${PRIVATE_SALT_ROOT}/gse260666-gse267145-${SLURM_JOB_ID}.bin
TARGET=${ROOT}/executions/model-data-068-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" || -e "${PRIVATE_SALT}" ]]; then
  printf 'Refusing to overwrite GSE260666 external-transfer artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0700 "${PRIVATE_SALT_ROOT}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/build_gse260666_external_transfer_fixture_cpu.sbatch

"${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse260666_external_transfer_fixture \
  tests.unit.test_evaluate_gse260666_external_histology_transfer \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" - "${PRIVATE_SALT}" <<'PY'
from pathlib import Path
import os
import secrets
import sys

path = Path(sys.argv[1])
descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
with os.fdopen(descriptor, "wb") as handle:
    handle.write(secrets.token_bytes(32))
    handle.flush()
    os.fsync(handle.fileno())
PY
"${PYTHON}" - "${ACTIVATION}" "${MATRIX}" "${TRAINING_MOLECULAR}" "${TRAINING_JOIN}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("all upstream frozen-tree verifications passed")
PY
"${PYTHON}" scripts/build_gse260666_external_transfer_fixture.py \
  --activation "${ACTIVATION}" \
  --activation-artifacts-sha256 734d45aee49221ae9b72c2b79de6b4959ae1d3d01ac85420c89135c05469e8f7 \
  --matrix "${MATRIX}" \
  --matrix-artifacts-sha256 3f934ab9c6ee730c6bd0887b10d963f8210f274c36dc0878cb814c85026dc038 \
  --training-molecular "${TRAINING_MOLECULAR}" \
  --training-molecular-artifacts-sha256 2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456 \
  --training-join "${TRAINING_JOIN}" \
  --training-join-artifacts-sha256 e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  --task-spec "${TASK}" \
  --task-spec-sha256 dcd4a9e83c83d625d03c19b73de21185380aa66118fb51d108ce74069b813e40 \
  --project-root "${PROJECT}" \
  --benchmark-root "${ROOT}" \
  --fingerprint-salt-file "${PRIVATE_SALT}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/build_gse260666_external_transfer_fixture.py \
  scripts/evaluate_gse260666_external_histology_transfer.py \
  tests/unit/test_build_gse260666_external_transfer_fixture.py \
  tests/unit/test_evaluate_gse260666_external_histology_transfer.py \
  slurm/build_gse260666_external_transfer_fixture_cpu.sbatch \
  "${TASK}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${MATRIX}/ARTIFACTS.json" \
  "${TRAINING_MOLECULAR}/ARTIFACTS.json" \
  "${TRAINING_JOIN}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
fixture = root / "fixture"
model = fixture / "model_input"
evaluator = fixture / "evaluator_only"
audit = fixture / "contamination_audit"
counts = np.load(model / "rna_counts.npy", mmap_mode="r", allow_pickle=False)
mask = np.load(model / "rna_observed_mask.npy", allow_pickle=False)
with (model / "participant_axis.tsv").open(encoding="utf-8", newline="") as handle:
    participants = list(csv.DictReader(handle, delimiter="\t"))
with (evaluator / "labels.tsv").open(encoding="utf-8", newline="") as handle:
    labels = list(csv.DictReader(handle, delimiter="\t"))
model_receipt = json.loads((model / "receipt.json").read_text(encoding="utf-8"))
audit_receipt = json.loads((audit / "audit.json").read_text(encoding="utf-8"))
if (
    counts.shape[0] != 16
    or counts.dtype != np.uint64
    or mask.dtype != np.bool_
    or not mask.all()
    or len(participants) != 16
    or len(labels) != 16
    or [row["row_id"] for row in participants] != [row["row_id"] for row in labels]
    or model_receipt.get("labels_included") is not False
    or audit_receipt.get("status") != "passed"
    or audit_receipt.get("near_duplicate_detected") is not False
    or audit_receipt.get("accession_bioproject_biosample_sra_or_participant_alias_overlap_count") != 0
    or audit_receipt.get("permissioned_random_salt_outside_release") is not True
    or audit_receipt.get("raw_expression_fingerprints_retained") is not False
):
    raise SystemExit("GSE260666 external-transfer fixture verification differs")
freeze_tree(
    model,
    {
        "artifact_class": "gse260666_external_transfer_model_input",
        "participants": 16,
        "biological_unit": "participant",
        "labels_included": False,
        "external_fit_or_calibration_performed": False,
        "status": "passed",
    },
)
verify_frozen_tree(model)
freeze_tree(
    evaluator,
    {
        "artifact_class": "gse260666_external_transfer_evaluator_only_outcomes",
        "participants": 16,
        "biological_unit": "participant",
        "model_environment_access": False,
        "nafl_nash_pooled": False,
        "status": "passed",
    },
)
verify_frozen_tree(evaluator)
freeze_tree(
    audit,
    {
        "artifact_class": "gse260666_gse267145_cross_cohort_contamination_audit",
        "participants_external": 16,
        "participants_training": 99,
        "raw_expression_fingerprints_retained": False,
        "status": "passed",
    },
)
verify_frozen_tree(audit)
freeze_tree(
    fixture,
    {
        "artifact_class": "gse260666_external_histology_transfer_fixture",
        "participants": 16,
        "evaluation_role": "external_development",
        "project_sealed": False,
        "selection_lock_present": False,
        "outcomes_scored": False,
        "status": "passed",
    },
)
verify_frozen_tree(fixture)
freeze_tree(
    root,
    {
        "artifact_class": "gse260666_external_histology_transfer_fixture_campaign",
        "participants": 16,
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/fixture/ARTIFACTS.json" \
  "${TARGET}/fixture/model_input/ARTIFACTS.json" \
  "${TARGET}/fixture/evaluator_only/ARTIFACTS.json" \
  "${TARGET}/fixture/contamination_audit/ARTIFACTS.json"
