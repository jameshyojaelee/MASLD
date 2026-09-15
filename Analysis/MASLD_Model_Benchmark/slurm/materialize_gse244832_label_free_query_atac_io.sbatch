#!/usr/bin/env bash
#SBATCH --job-name=model-data-405
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
CONFIG="${ROOT}/config/campaigns/gse244832_label_free_query_atac_materialization_20260825.json"
: "${VALIDATION_JOB_ID:?VALIDATION_JOB_ID must identify the afterok validation job}"
VALIDATION="${ROOT}/executions/model-check-335-${VALIDATION_JOB_ID}"
VALID_TARGET="${ROOT}/executions/gse244832-label-free-atac-query-valid-context-20260825"
TEST_TARGET="${ROOT}/executions/gse244832-label-free-atac-query-test-context-20260825"
VALID_STAGE="${VALID_TARGET}.staging-${SLURM_JOB_ID}"
TEST_STAGE="${TEST_TARGET}.staging-${SLURM_JOB_ID}"
VALID_FAILED="${VALID_TARGET}.failed-${SLURM_JOB_ID}"
TEST_FAILED="${TEST_TARGET}.failed-${SLURM_JOB_ID}"
preserve_failed_stages() {
  if [[ -d "${VALID_STAGE}" && ! -e "${VALID_FAILED}" ]]; then mv "${VALID_STAGE}" "${VALID_FAILED}"; fi
  if [[ -d "${TEST_STAGE}" && ! -e "${TEST_FAILED}" ]]; then mv "${TEST_STAGE}" "${TEST_FAILED}"; fi
}
trap preserve_failed_stages EXIT
for path in "${VALID_TARGET}" "${TEST_TARGET}" "${VALID_STAGE}" "${TEST_STAGE}" "${VALID_FAILED}" "${TEST_FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite query-ATAC artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
cd "${ROOT}"

for pair in \
  "dfbf1a2db5b9d0eb4768081e61145360616e932720b036eed7fcfd09bebb4788 config/campaigns/gse244832_label_free_query_atac_materialization_20260825.json" \
  "1c6a99e7eb69b8c358da804d1b00730502b00c7a7b29fa9df2f790d62458749f scripts/materialize_gse244832_label_free_query_atac.py" \
  "9db6a001dde8170d7eea214b7dbbddf7c3b20b317ae65ca1c2b56242de8e6c81 tests/unit/test_materialize_gse244832_label_free_query_atac.py" \
  "ff1a506d5cca81382911355978097ecd49fcb1b2b48826e8e760a5f2ad1ad35e executions/gse244832-reference-guarded-atac-exchange-axis-20260825/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${VALIDATION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
metadata = verify_frozen_tree(root).get("metadata", {})
if (
    metadata.get("artifact_class") != "gse244832_label_free_query_atac_validation"
    or metadata.get("sequence_execution_authorized") is not False
    or metadata.get("evaluator_outcome_artifact_read") is not False
    or metadata.get("status") != "passed"
):
    raise SystemExit("query validation authority differs")
PY

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" \
  scripts/materialize_gse244832_label_free_query_atac.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --valid-output "${VALID_STAGE}" \
  --test-output "${TEST_STAGE}" \
  --workers 6

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${ROOT}" "${VALID_STAGE}" "${TEST_STAGE}" <<'PY'
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
valid_root = Path(sys.argv[2]).resolve(strict=True)
test_root = Path(sys.argv[3]).resolve(strict=True)
axis_root = root / "executions/gse244832-reference-guarded-atac-exchange-axis-20260825"
axis_roles: dict[str, set[str]] = {"valid": set(), "test": set()}
with (axis_root / "windows.tsv").open("r", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        axis_roles[row["role"]].add(row["contig"])
if axis_roles["valid"] & axis_roles["test"]:
    raise SystemExit("axis role contigs overlap")

for artifact, role, target_role in (
    (valid_root, "valid", "test"),
    (test_root, "test", "valid"),
):
    metadata = verify_frozen_tree(artifact).get("metadata", {})
    contract = json.loads((artifact / "query_contract.json").read_text(encoding="utf-8"))
    authorities = json.loads((artifact / "input_authorities.json").read_text(encoding="utf-8"))
    with (artifact / "query_axis.tsv").open("r", encoding="utf-8", newline="") as handle:
        axis_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (artifact / "query_windows.tsv").open("r", encoding="utf-8", newline="") as handle:
        window_rows = list(csv.DictReader(handle, delimiter="\t"))
    profiles = np.load(artifact / "query_unit_fragment_profile.uint32.npy", mmap_mode="r", allow_pickle=False)
    counts = np.load(artifact / "query_unit_fragment_total.uint32.npy", mmap_mode="r", allow_pickle=False)
    if profiles.shape != (18, 4, 16000, 20) or profiles.dtype != np.uint32:
        raise SystemExit(f"profile contract differs: {artifact}")
    if counts.shape != (18, 4, 16000) or counts.dtype != np.uint32:
        raise SystemExit(f"count contract differs: {artifact}")
    if not np.array_equal(profiles.sum(axis=3, dtype=np.uint64), counts):
        raise SystemExit(f"profile/count conservation differs: {artifact}")
    observed_contigs = {row["contig"] for row in window_rows}
    if observed_contigs != axis_roles[role] or observed_contigs & axis_roles[target_role]:
        raise SystemExit(f"context/scored-target contigs differ: {artifact}")
    if len(window_rows) != 16000 or len(axis_rows) != 72:
        raise SystemExit(f"query axes differ: {artifact}")
    if (
        metadata.get("artifact_class") != "gse244832_label_free_query_atac"
        or metadata.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
        or metadata.get("condition_values_read") is not False
        or metadata.get("rna_assay_opened") is not False
        or metadata.get("evaluator_outcome_artifact_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("sequence_execution_authorized") is not False
        or contract.get("scored_target_contigs_exported_in_query_values") is not False
        or contract.get("cross_assay_join_inferred") is not False
        or "processing_qc.tsv" in contract.get("model_input_files", [])
        or authorities.get("condition_or_phenotype_authority_included") is not False
        or authorities.get("rna_or_cross_assay_authority_included") is not False
        or authorities.get("evaluator_outcome_authority_included") is not False
        or authorities.get("sequence_execution_authorized") is not False
        or sorted(authorities.get("fragment_sha256_by_donor", {}))
        != [f"D{index:02d}" for index in range(1, 19)]
    ):
        raise SystemExit(f"query firewall or authority differs: {artifact}")
PY

trap - EXIT
mv "${VALID_STAGE}" "${VALID_TARGET}"
mv "${TEST_STAGE}" "${TEST_TARGET}"
printf 'Valid-context query artifact: %s\n' "${VALID_TARGET}"
sha256sum "${VALID_TARGET}/ARTIFACTS.json"
printf 'Test-context query artifact: %s\n' "${TEST_TARGET}"
sha256sum "${TEST_TARGET}/ARTIFACTS.json"
