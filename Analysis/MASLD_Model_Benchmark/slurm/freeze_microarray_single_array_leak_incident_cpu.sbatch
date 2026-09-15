#!/usr/bin/env bash
#SBATCH --job-name=model-check-090
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Freeze the incident record for the fRMA multi-array summarization leak.
#
# The lesson generalizes past GPL570: any summarizer that fits a robust or
# otherwise shared model across the columns it is handed borrows across held
# arrays, however frozen its reference vectors are. oligo RMA and SCAN on
# GPL16686/GSE83452 are in that class and are named here as untested.
#
# Every number in the record is re-derived from the bound evidence rather than
# restated, so the conclusion can be checked instead of trusted.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
PYTHON_MODULE=python/3.11.5-GCCcore-13.2.0

FAILING_RECEIPT=${ROOT}/executions/model-check-085-21107058-firewall.failed/firewall/frma_single_array_firewall.json
FAILING_RECEIPT_SHA256=65e9ab1ff2021b9d3ba38849c4e73cecb9469942e7e0fffdcd8d4b0aac5d4903
MECHANISM=${ROOT}/executions/model-check-087-21107439-diagnostic/frma_batch_dependence_diagnostic.json
MECHANISM_SHA256=27663ce8395b3b25506f7a6098839dd4539fb6631377a7c6b84d083a25d38122
AGREEMENT=${ROOT}/executions/model-check-088-21107642-agreement/frma_summarize_rule_agreement.json
AGREEMENT_SHA256=a2844be4f80ca4135fbf7b26c778ea26a5a949809c96261adc480cab70aefa6d

FIREWALL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-089-21107882-firewall/firewall
FIREWALL_SHA256=c292dadc8454f2d13c682205a3dbf2e738082e591fa10ff0d6a9c4806f3a0e29
PASSING_VERDICT=${FIREWALL}/firewall_verdict.json

PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e

TARGET=${ROOT}/executions/model-check-090-${SLURM_JOB_ID}-incident
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed incident stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite incident artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation"

printf '%s  %s\n' \
  "${FAILING_RECEIPT_SHA256}" "${FAILING_RECEIPT}" \
  "${MECHANISM_SHA256}"       "${MECHANISM}" \
  "${AGREEMENT_SHA256}"       "${AGREEMENT}" \
  "${FIREWALL_SHA256}"        "${FIREWALL}/ARTIFACTS.json" \
  "${PREPROCESSING_SHA256}"   "${PREPROCESSING}" \
  | sha256sum --check --strict

cd "${ROOT}"
module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"

python3 -m unittest -v \
  tests.unit.test_build_microarray_single_array_leak_incident \
  tests.unit.test_gse49541_frma_firewall_receipt \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

python3 scripts/build_microarray_single_array_leak_incident.py \
  --failing-receipt "${FAILING_RECEIPT}" \
  --mechanism-diagnostic "${MECHANISM}" \
  --agreement-diagnostic "${AGREEMENT}" \
  --passing-verdict "${PASSING_VERDICT}" \
  --output "${STAGE}/incident" \
  > "${STAGE}/validation/incident_stdout.json" \
  2> "${STAGE}/validation/incident_stderr.txt"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

sha256sum \
  scripts/build_microarray_single_array_leak_incident.py \
  tests/unit/test_build_microarray_single_array_leak_incident.py \
  slurm/freeze_microarray_single_array_leak_incident_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  "${FAILING_RECEIPT}" \
  "${MECHANISM}" \
  "${AGREEMENT}" \
  "${FIREWALL}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'Incident stage contains a symlink\n' >&2
  exit 1
fi

python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
incident_dir = stage / "incident"
incident = json.loads(
    (incident_dir / "single_array_leak_incident.json").read_text(encoding="utf-8")
)
if (
    incident["status"] != "resolved_by_holding_the_frma_batch_at_one_array"
    or incident["series"] != "GSE49541"
    or incident["platform_id"] != "GPL570"
    or incident["resolution"]["summarize_rule_retained"] != "robust_weighted_average"
    or incident["mechanism"]["batch_dependent_summarize_rules"] != ["robust_weighted_average"]
    or incident["failure_evidence"]["batch_membership_matters"] is not True
    or incident["generalization"]["frozen_reference_vectors_are_not_sufficient"] is not True
    or incident["generalization"]["single_array_capability_probe_is_not_sufficient"] is not True
    or len(incident["bound_evidence"]) != 4
    or incident["labels_read"] is not False
    or incident["model_training_activated"] is not False
    or incident["sealed_outcomes_read"] is not False
):
    raise SystemExit("incident record differs from the bound evidence")

metadata = {
    "artifact_class": "microarray_single_array_summarization_leak_incident",
    "incident_id": incident["incident_id"],
    "series": "GSE49541",
    "platform_id": "GPL570",
    "method": "frma_with_exact_hgu133plus2frmavecs",
    "batch_dependent_summarize_rules": incident["mechanism"]["batch_dependent_summarize_rules"],
    "batch_invariant_summarize_rules": incident["mechanism"]["batch_invariant_summarize_rules"],
    "maximum_absolute_difference_log2": incident["failure_evidence"][
        "maximum_absolute_difference_log2"
    ],
    "maximum_differing_elements": incident["failure_evidence"]["maximum_differing_elements"],
    "alternative_rule_cost_relative_to_leak": incident["alternative_rule_cost_relative_to_leak"],
    "resolution": "hold_the_frma_batch_at_one_array",
    "untested_platforms_in_the_same_hazard_class": incident["generalization"][
        "must_be_tested_before_use"
    ],
    "labels_read": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
}
freeze_tree(incident_dir, metadata)
verify_frozen_tree(incident_dir)
freeze_tree(stage, metadata | {"artifact_class": "microarray_single_array_leak_incident_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Microarray single-array leak incident: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/incident/ARTIFACTS.json"
