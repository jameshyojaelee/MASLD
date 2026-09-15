#!/usr/bin/env bash
#SBATCH --job-name=model-check-338
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
CONFIG="${ROOT}/config/campaigns/peakvi_gse296875_to_gse281367_fit_predict_revision2_20260825.json"
CONFIG_SHA256=a8c949c33469970d8f07fe98fd1464c2622d1ebee3f6928fa6213d62c2d3c795
TARGET="${ROOT}/executions/peakvi-cross-cohort-revision2-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite PeakVI revision-2 admission\n' >&2
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
bash -n slurm/run_peakvi_cross_cohort_revision2_bundle.sbatch
sha256sum \
  "${CONFIG}" \
  scripts/peakvi_cross_cohort_fit_predict.py \
  scripts/commit_gse281367_atac_transport_predictions.py \
  tests/unit/test_peakvi_cross_cohort_fit_predict.py \
  tests/unit/test_commit_gse281367_atac_transport_predictions.py \
  slurm/run_peakvi_cross_cohort_bundle.sbatch \
  slurm/run_peakvi_cross_cohort_revision2_bundle.sbatch \
  slurm/admit_peakvi_cross_cohort_revision2_cpu.sbatch \
  config/artifacts/incidents/peakvi_target_cell_hash_collision_21101008.json \
  > "${STAGE}/source.sha256"
PYTHONPATH="src:." "${PYTHON}" - \
  "${ROOT}" "${CONFIG}" "${STAGE}" <<'PY'
from hashlib import sha256
from pathlib import Path
import json
import sys

import numpy as np
from scipy import sparse

from masld_bench.artifacts import freeze_tree
from scripts.peakvi_cross_cohort_fit_predict import (
    ROTATIONS,
    cell_rows,
    digest,
    feature_rows,
    validate_config,
)

root = Path(sys.argv[1]).resolve(strict=True)
config_path = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
config, paths = validate_config(root, config_path)
incident_record = config["incident_authority"]
incident_path = (root / incident_record["path"]).resolve(strict=True)
incident = json.loads(incident_path.read_text(encoding="utf-8"))
if (
    digest(incident_path) != incident_record["sha256"]
    or incident.get("status") != "terminal_failed_closed"
    or incident.get("correction", {}).get("old_target_artifacts_eligible_for_modeling")
    is not False
    or incident.get("firewall", {}).get("model_fit_started") is not False
    or incident.get("firewall", {}).get("outcomes_read") is not False
):
    raise SystemExit("PeakVI row-identity incident authority differs")
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
        or any(
            not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0)
            for matrix in (source_context, source_targets, target_context)
        )
    ):
        raise SystemExit(f"PeakVI revision-2 matrix contract differs: {rotation_id}")
    checks.append(
        {
            "rotation_id": rotation_id,
            "source_cells": len(source_cells),
            "target_cells": len(target_cells),
            "source_cell_hashes_unique": len(
                {row["cell_hash"] for row in source_cells}
            )
            == len(source_cells),
            "target_cell_hashes_unique": len(
                {row["cell_hash"] for row in target_cells}
            )
            == len(target_cells),
            "context_features": len(source_context_features),
            "target_features": len(source_target_features),
            "source_context_nnz": int(source_context.nnz),
            "source_target_nnz": int(source_targets.nnz),
            "target_context_nnz": int(target_context.nnz),
            "source_zero_context_features": int(
                np.count_nonzero(source_context.getnnz(axis=0) == 0)
            ),
            "source_zero_target_features": int(
                np.count_nonzero(source_targets.getnnz(axis=0) == 0)
            ),
            "context_target_contigs_disjoint": True,
        }
    )
receipt = {
    "schema_version": "masld-bench-peakvi-cross-cohort-revision2-admission-v1",
    "status": "passed",
    "campaign_sha256": sha256(config_path.read_bytes()).hexdigest(),
    "incident_sha256": incident_record["sha256"],
    "rotations": checks,
    "unit_tests": 8,
    "seeds": config["seeds"],
    "native_output_family": "joint_representation",
    "projected_output_family": "masked_accessibility_count",
    "strict_state_restore_required": True,
    "target_query_adaptation_allowed": False,
    "query_atac_context_values_read_for_contract": True,
    "source_training_values_read_for_contract": True,
    "target_role_values_read": False,
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
        "revision": 2,
        "campaign_sha256": receipt["campaign_sha256"],
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "cell_hashes_globally_unique": True,
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
