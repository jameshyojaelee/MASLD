#!/usr/bin/env bash
#SBATCH --job-name=model-check-336
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
CONFIG="${ROOT}/config/campaigns/peakvi_gse296875_to_gse281367_fit_predict_20260825.json"
CONFIG_SHA256=0b10777b59d84b81672f89c1977c2159aae74205bc61acdbe08a990f2a1f80c3
TARGET="${ROOT}/executions/peakvi-cross-cohort-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite PeakVI cross-cohort admission\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${CONFIG_SHA256}" "${CONFIG}" | sha256sum --check --strict
cd "${ROOT}"
PYTHONPATH="src:." "${PYTHON}" -m unittest -v \
  tests.unit.test_peakvi_cross_cohort_fit_predict \
  tests.unit.test_commit_gse281367_atac_transport_predictions \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"
sha256sum \
  "${CONFIG}" \
  scripts/peakvi_cross_cohort_fit_predict.py \
  scripts/commit_gse281367_atac_transport_predictions.py \
  tests/unit/test_peakvi_cross_cohort_fit_predict.py \
  tests/unit/test_commit_gse281367_atac_transport_predictions.py \
  slurm/run_peakvi_cross_cohort_bundle.sbatch \
  slurm/admit_peakvi_cross_cohort_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src:." "${PYTHON}" - \
  "${ROOT}" "${CONFIG}" "${STAGE}" <<'PY'
from pathlib import Path
from hashlib import sha256
import json
import sys

from scipy import sparse

from masld_bench.artifacts import freeze_tree
from scripts.peakvi_cross_cohort_fit_predict import (
    ROTATIONS,
    cell_rows,
    feature_rows,
    validate_config,
)

root = Path(sys.argv[1]).resolve(strict=True)
config_path = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
config, paths = validate_config(root, config_path)
source_cells = cell_rows(paths["source"] / "cell_axis.tsv", 39)
checks = []
for rotation_id, (context_role, target_role) in ROTATIONS.items():
    target_root = paths[f"target_{context_role}"]
    target_cells = cell_rows(target_root / "cell_axis.tsv", 12)
    source_context_features = feature_rows(
        paths["source"] / f"{context_role}.windows.tsv"
    )
    target_context_features = feature_rows(target_root / "context.windows.tsv")
    source_target_features = feature_rows(
        paths["source"] / f"{target_role}.windows.tsv"
    )
    source_context = sparse.load_npz(
        paths["source"] / f"{context_role}.counts.npz"
    ).tocsr()
    source_targets = sparse.load_npz(
        paths["source"] / f"{target_role}.counts.npz"
    ).tocsr()
    target_context = sparse.load_npz(target_root / "context.counts.npz").tocsr()
    if (
        source_context_features != target_context_features
        or {row["contig"] for row in source_context_features}
        & {row["contig"] for row in source_target_features}
        or source_context.shape != (len(source_cells), 16_000)
        or source_targets.shape != source_context.shape
        or target_context.shape != (len(target_cells), 16_000)
        or any(
            (matrix.getnnz(axis=1) == 0).any()
            for matrix in (source_context, source_targets, target_context)
        )
    ):
        raise SystemExit(f"PeakVI admitted matrix contract differs: {rotation_id}")
    checks.append(
        {
            "rotation_id": rotation_id,
            "source_cells": len(source_cells),
            "target_cells": len(target_cells),
            "context_features": len(source_context_features),
            "target_features": len(source_target_features),
            "source_context_nnz": int(source_context.nnz),
            "source_target_nnz": int(source_targets.nnz),
            "target_context_nnz": int(target_context.nnz),
            "context_target_contigs_disjoint": True,
        }
    )
receipt = {
    "schema_version": "masld-bench-peakvi-cross-cohort-admission-receipt-v1",
    "status": "passed",
    "campaign_sha256": sha256(config_path.read_bytes()).hexdigest(),
    "rotations": checks,
    "seeds": config["seeds"],
    "native_output_family": "joint_representation",
    "projected_output_family": "masked_accessibility_count",
    "strict_state_restore_required": True,
    "target_query_adaptation_allowed": False,
    "condition_labels_read": False,
    "phenotype_values_read": False,
    "development_outcomes_read": False,
    "sealed_data_read": False,
    "metrics_calculated": False,
    "scorer_invoked": False,
    "champion_claim_allowed": False,
    "next_gate": "dependency_queue_registration_then_source_fit_blind_prediction_and_commit",
}
(stage / "receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(
    stage,
    {
        "artifact_class": "peakvi_cross_cohort_fit_predict_admission",
        "campaign_sha256": receipt["campaign_sha256"],
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
