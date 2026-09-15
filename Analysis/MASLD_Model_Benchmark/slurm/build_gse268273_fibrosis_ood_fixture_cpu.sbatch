#!/usr/bin/env bash
#SBATCH --job-name=model-data-074
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
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
SOURCE=${ROOT}/executions/bulk-expansion-source-audit-21065762
MATRIX=${ROOT}/executions/bulk-expansion-matrix-audit-21065853
ACTIVATION=${ROOT}/executions/gse268273-activation-21066074
TRAINING_MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
TRAINING_JOIN=${ROOT}/executions/gse267145-authoritative-join-21064930
GSE260666_ACTIVATION=${ROOT}/executions/gse260666-activation-21066072
MICROARRAY_PARTICIPANTS=${ROOT}/executions/model-data-067-21080193-participants
MICROARRAY_OVERLAP=${ROOT}/executions/model-data-072-21080651-overlap
TASK=${ROOT}/config/evaluation/gse268273_fibrosis_ood_transfer_task.toml
GATE=${ROOT}/config/evaluation/gse268273_fibrosis_ood_promotion_gate.json
EXPANSION=${ROOT}/config/evaluation/cross_cohort_expansion.toml
PRIVATE_SALT_ROOT=${ROOT}/executions/private/fingerprint-salts
PRIVATE_SALT=${PRIVATE_SALT_ROOT}/gse268273-ood-${SLURM_JOB_ID}.bin
TARGET=${ROOT}/executions/model-data-074-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" || -e "${PRIVATE_SALT}" ]]; then
  printf 'Refusing to overwrite GSE268273 OOD fixture artifacts.\n' >&2
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
bash -n slurm/build_gse268273_fibrosis_ood_fixture_cpu.sbatch

"${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse268273_fibrosis_ood_fixture \
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
"${PYTHON}" - \
  "${SOURCE}" \
  "${MATRIX}" \
  "${ACTIVATION}" \
  "${TRAINING_MOLECULAR}" \
  "${TRAINING_JOIN}" \
  "${GSE260666_ACTIVATION}" \
  "${MICROARRAY_PARTICIPANTS}" \
  "${MICROARRAY_OVERLAP}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("all upstream frozen-tree verifications passed")
PY
"${PYTHON}" scripts/build_gse268273_fibrosis_ood_fixture.py \
  --source "${SOURCE}" \
  --matrix "${MATRIX}" \
  --activation "${ACTIVATION}" \
  --training-molecular "${TRAINING_MOLECULAR}" \
  --training-join "${TRAINING_JOIN}" \
  --gse260666-activation "${GSE260666_ACTIVATION}" \
  --microarray-participants "${MICROARRAY_PARTICIPANTS}" \
  --microarray-overlap "${MICROARRAY_OVERLAP}" \
  --task-spec "${TASK}" \
  --promotion-gate "${GATE}" \
  --expansion-config "${EXPANSION}" \
  --project-root "${PROJECT}" \
  --benchmark-root "${ROOT}" \
  --fingerprint-salt-file "${PRIVATE_SALT}" \
  --output "${STAGE}/fixture" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/build_gse268273_fibrosis_ood_fixture.py \
  tests/unit/test_build_gse268273_fibrosis_ood_fixture.py \
  slurm/build_gse268273_fibrosis_ood_fixture_cpu.sbatch \
  "${TASK}" \
  "${GATE}" \
  "${EXPANSION}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${MATRIX}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${TRAINING_MOLECULAR}/ARTIFACTS.json" \
  "${TRAINING_JOIN}/ARTIFACTS.json" \
  "${GSE260666_ACTIVATION}/ARTIFACTS.json" \
  "${MICROARRAY_PARTICIPANTS}/ARTIFACTS.json" \
  "${MICROARRAY_OVERLAP}/ARTIFACTS.json" \
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
source = fixture / "source_evidence"
preprocessing = fixture / "preprocessing_input"
model = fixture / "model_input"
evaluator = fixture / "evaluator_only"
audit = fixture / "contamination_audit"
summary = json.loads((fixture / "receipt.json").read_text(encoding="utf-8"))
source_receipt = json.loads((source / "receipt.json").read_text(encoding="utf-8"))
preprocessing_receipt = json.loads((preprocessing / "receipt.json").read_text(encoding="utf-8"))
model_receipt = json.loads((model / "receipt.json").read_text(encoding="utf-8"))
evaluator_receipt = json.loads((evaluator / "receipt.json").read_text(encoding="utf-8"))
audit_receipt = json.loads((audit / "audit.json").read_text(encoding="utf-8"))
mask = np.load(model / "rna_observed_mask.npy", allow_pickle=False)
with (model / "participant_axis.tsv").open(encoding="utf-8", newline="") as handle:
    participants = list(csv.DictReader(handle, delimiter="\t"))
with (model / "candidate_gene_axis.tsv").open(encoding="utf-8", newline="") as handle:
    genes = list(csv.DictReader(handle, delimiter="\t"))
with (preprocessing / "raw_run_manifest.tsv").open(encoding="utf-8", newline="") as handle:
    runs = list(csv.DictReader(handle, delimiter="\t"))
with (evaluator / "outcomes.tsv").open(encoding="utf-8", newline="") as handle:
    outcomes = list(csv.DictReader(handle, delimiter="\t"))
if (
    summary.get("participants") != 109
    or summary.get("geo_records_are_unique_participants") is not True
    or summary.get("technical_runs_are_replicates") is not False
    or summary.get("model_input_contains_labels") is not False
    or summary.get("model_input_contains_expression_values") is not False
    or summary.get("deposited_voom_scored_input_allowed") is not False
    or summary.get("processed_differential_expression_opened") is not False
    or summary.get("external_outcomes_scored") is not False
    or source_receipt.get("unique_participants") != 109
    or source_receipt.get("ena", {}).get("runs") != 824
    or source_receipt.get("multiple_runs_are_technical_partitions") is not True
    or preprocessing_receipt.get("raw_fastq_materialized") is not False
    or preprocessing_receipt.get("technical_runs_collapsed_to_participant") is not False
    or model_receipt.get("labels_included") is not False
    or model_receipt.get("rna_values_present") is not False
    or model_receipt.get("candidate_unique_v49_genes") != 14078
    or model_receipt.get("candidate_aggregated_v49_genes") != 11
    or model_receipt.get("duplicate_mapping_resolution") != "deterministic_sum_raw_counts_before_normalization"
    or evaluator_receipt.get("mapped_to_gse267145_stage3") is not False
    or evaluator_receipt.get("outcomes_scored") is not False
    or audit_receipt.get("status") != "passed"
    or audit_receipt.get("cohort_family_overlap_detected") is not False
    or audit_receipt.get("near_duplicate_source_ids") != []
    or audit_receipt.get("within_gse268273", {}).get("near_duplicate_detected") is not False
    or audit_receipt.get("raw_expression_fingerprints_retained") is not False
    or len(participants) != 109
    or len(genes) != 14078
    or sum(int(row["source_feature_count"]) > 1 for row in genes) != 11
    or len(runs) != 824
    or len(outcomes) != 109
    or mask.shape != (109,)
    or mask.dtype != np.bool_
    or mask.any()
):
    raise SystemExit("GSE268273 OOD fixture verification differs")
freeze_tree(
    source,
    {
        "artifact_class": "gse268273_official_source_topology_evidence",
        "participants": 109,
        "technical_runs": 824,
        "status": "passed",
    },
)
verify_frozen_tree(source)
freeze_tree(
    preprocessing,
    {
        "artifact_class": "gse268273_raw_preprocessing_input",
        "participants": 109,
        "technical_runs": 824,
        "labels_included": False,
        "raw_materialized": False,
        "status": "blocked_raw_materialization",
    },
)
verify_frozen_tree(preprocessing)
freeze_tree(
    model,
    {
        "artifact_class": "gse268273_fibrosis_ood_model_input",
        "participants": 109,
        "labels_included": False,
        "rna_values_present": False,
        "scored_prediction_allowed": False,
        "status": "blocked_raw_reprocessing",
    },
)
verify_frozen_tree(model)
freeze_tree(
    evaluator,
    {
        "artifact_class": "gse268273_fibrosis_ood_evaluator_only_outcomes",
        "participants": 109,
        "model_environment_access": False,
        "outcomes_scored": False,
        "status": "passed",
    },
)
verify_frozen_tree(evaluator)
freeze_tree(
    audit,
    {
        "artifact_class": "gse268273_cross_cohort_contamination_audit",
        "participants": 109,
        "raw_expression_fingerprints_retained": False,
        "status": "passed",
    },
)
verify_frozen_tree(audit)
freeze_tree(
    fixture,
    {
        "artifact_class": "gse268273_fibrosis_etiology_ood_fixture",
        "participants": 109,
        "evaluation_role": "project_exposed_external_development",
        "project_sealed": False,
        "model_fitted": False,
        "metrics_calculated": False,
        "scored_use_blocked": True,
        "status": "passed",
    },
)
verify_frozen_tree(fixture)
freeze_tree(
    root,
    {
        "artifact_class": "gse268273_fibrosis_etiology_ood_fixture_campaign",
        "participants": 109,
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
  "${TARGET}/fixture/source_evidence/ARTIFACTS.json" \
  "${TARGET}/fixture/preprocessing_input/ARTIFACTS.json" \
  "${TARGET}/fixture/model_input/ARTIFACTS.json" \
  "${TARGET}/fixture/evaluator_only/ARTIFACTS.json" \
  "${TARGET}/fixture/contamination_audit/ARTIFACTS.json"
