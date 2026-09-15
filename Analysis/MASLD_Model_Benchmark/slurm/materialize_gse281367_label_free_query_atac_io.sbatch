#!/usr/bin/env bash
#SBATCH --job-name=model-data-402
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
CONFIG="${ROOT}/config/campaigns/gse281367_label_free_query_atac_materialization_20260825.json"
VALIDATION="${ROOT}/executions/gse281367-label-free-query-atac-validation-21100549"
VALID_TARGET="${ROOT}/executions/gse281367-label-free-atac-query-valid-context-20260825"
TEST_TARGET="${ROOT}/executions/gse281367-label-free-atac-query-test-context-20260825"
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
  "01c2c38d7583fa16c06444b8870e6bccd181f4e2d1ad57f05749da820cdbe51b config/campaigns/gse281367_label_free_query_atac_materialization_20260825.json" \
  "164847977304d0d545e6448ed8a7371bce4a353ca0d64ff52273f0557e985776 scripts/materialize_gse281367_label_free_query_atac.py" \
  "88e957fc3ff3ab4dcb70eed25c94d4d8601369e352655587ff53c9da5e4daefa tests/unit/test_materialize_gse281367_label_free_query_atac.py" \
  "02e7d6c2a317f8b539f9abf28c7e02ec6d9c4a4af46c40341d58746bf2c8e7ec executions/gse281367-label-free-query-atac-validation-21100549/ARTIFACTS.json" \
  "c7a5eaf49784aa76b13a1cd9edaa94220775c7378ab9272fb0c2bdcd1d3a9873 executions/gse281367-label-free-atac-exchange-axis-20260825/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" \
  scripts/materialize_gse281367_label_free_query_atac.py \
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
axis_root = root / "executions/gse281367-label-free-atac-exchange-axis-20260825"

axis_roles: dict[str, set[str]] = {"valid": set(), "test": set()}
with (axis_root / "windows.tsv").open("r", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        axis_roles[row["role"]].add(row["contig"])
if axis_roles["valid"] & axis_roles["test"]:
    raise SystemExit("axis role contigs overlap")

expected_campaign_sha = "01c2c38d7583fa16c06444b8870e6bccd181f4e2d1ad57f05749da820cdbe51b"
expected_builder_sha = "164847977304d0d545e6448ed8a7371bce4a353ca0d64ff52273f0557e985776"
for artifact, role, target_role in (
    (valid_root, "valid", "test"),
    (test_root, "test", "valid"),
):
    manifest = verify_frozen_tree(artifact)
    metadata = manifest.get("metadata", {})
    contract = json.loads((artifact / "query_contract.json").read_text(encoding="utf-8"))
    authorities = json.loads((artifact / "input_authorities.json").read_text(encoding="utf-8"))
    with (artifact / "query_axis.tsv").open("r", encoding="utf-8", newline="") as handle:
        axis_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (artifact / "query_windows.tsv").open("r", encoding="utf-8", newline="") as handle:
        window_rows = list(csv.DictReader(handle, delimiter="\t"))
    profiles = np.load(artifact / "query_unit_fragment_profile.uint32.npy", mmap_mode="r", allow_pickle=False)
    counts = np.load(artifact / "query_unit_fragment_total.uint32.npy", mmap_mode="r", allow_pickle=False)
    if profiles.shape != (12, 4, 16000, 20) or profiles.dtype != np.uint32:
        raise SystemExit(f"profile contract differs: {artifact}")
    if counts.shape != (12, 4, 16000) or counts.dtype != np.uint32:
        raise SystemExit(f"count contract differs: {artifact}")
    if not np.array_equal(profiles.sum(axis=3, dtype=np.uint64), counts):
        raise SystemExit(f"profile/count conservation differs: {artifact}")
    observed_contigs = {row["contig"] for row in window_rows}
    if observed_contigs != axis_roles[role] or observed_contigs & axis_roles[target_role]:
        raise SystemExit(f"context/scored-target contigs differ: {artifact}")
    if len(window_rows) != 16000 or len(axis_rows) != 48:
        raise SystemExit(f"query axes differ: {artifact}")
    if (
        metadata.get("artifact_class") != "gse281367_label_free_query_atac"
        or metadata.get("context_role") != role
        or metadata.get("scored_target_role") != target_role
        or metadata.get("condition_values_read") is not False
        or metadata.get("evaluator_outcome_artifact_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("champion_claim_allowed") is not False
        or contract.get("scored_target_contigs_exported_in_query_values") is not False
        or contract.get("condition_values_read") is not False
        or contract.get("phenotype_values_read") is not False
        or contract.get("evaluator_outcome_artifact_read") is not False
        or contract.get("metrics_calculated") is not False
        or contract.get("champion_claim_allowed") is not False
        or "processing_qc.tsv" in contract.get("model_input_files", [])
        or authorities.get("campaign_sha256") != expected_campaign_sha
        or authorities.get("builder_sha256") != expected_builder_sha
        or authorities.get("condition_or_phenotype_authority_included") is not False
        or authorities.get("evaluator_outcome_authority_included") is not False
        or sorted(authorities.get("fragment_sha256_by_donor", {})) != [f"Z{index:02d}" for index in range(1, 13)]
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
