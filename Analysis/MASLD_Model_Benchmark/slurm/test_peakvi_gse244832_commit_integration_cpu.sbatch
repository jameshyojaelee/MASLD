#!/usr/bin/env bash
#SBATCH --job-name=model-check-342
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
OUTPUT=${ROOT}/executions/model-check-342-${SLURM_JOB_ID}
STAGE=${OUTPUT}.staging

if [[ -e ${OUTPUT} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite GSE244832 commit integration fixture\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
sha256sum \
  scripts/peakvi_gse244832_inference_only.py \
  scripts/commit_peakvi_gse244832_atac_transport_predictions.py \
  scripts/commit_gse244832_atac_transport_predictions.py \
  slurm/test_peakvi_gse244832_commit_integration_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src:. "${PYTHON}" - "${ROOT}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import numpy as np
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.commit_peakvi_gse244832_atac_transport_predictions import commit_peakvi
from scripts.peakvi_gse244832_inference_only import digest, make_prediction_manifest

root, stage = map(Path, sys.argv[1:])
registration = root / "executions/model-check-337-21101061"
axis = root / "executions/gse244832-reference-guarded-atac-exchange-axis-20260825"
registration_sha = "194a0a9d96ce036ec1e67409dfd59623c079b706ca3547accf67d9f9278c0a6f"
axis_sha = "ff1a506d5cca81382911355978097ecd49fcb1b2b48826e8e760a5f2ad1ad35e"
contract = json.loads((registration / "execution_contract.json").read_text(encoding="utf-8"))
states = np.load(
    registration / contract["query_lineage_missing_state_path"], allow_pickle=False
)
config = {
    "seeds": [1103, 2207, 3301],
    "exchange_axis": {"artifacts_sha256": axis_sha},
    "execution_registration": {"artifacts_sha256": registration_sha},
    "query_authorities": {
        "valid": {"artifacts_sha256": "801411d062f9d9f199c0468b6c8c6256263af85ea2805d75e3d166299d190742"},
        "test": {"artifacts_sha256": "970992d3dc84f2d2c6c2c5357149579fa09b16051ddafb7eed0c1f2949a9ee8a"},
    },
}
queries = {
    "valid": root / "executions/gse244832-label-free-atac-query-valid-context-20260825",
    "test": root / "executions/gse244832-label-free-atac-query-test-context-20260825",
}
rotations = {
    "valid_context_test_target": "valid",
    "test_context_valid_target": "test",
}
source = stage / "source"
source_seed_records = {}
local_seed_records = {}
for rotation_id in rotations:
    source_rotation = source / "predictions" / rotation_id
    source_seed_records[rotation_id] = []
    local_seed_records[rotation_id] = []
    for seed in config["seeds"]:
        seed_root = source_rotation / f"seed-{seed}"
        seed_root.mkdir(parents=True)
        peakvi = seed_root / "peakvi.state_dict.pt"
        decoder = seed_root / "decoder.state_dict.pt"
        peakvi.write_bytes(f"peakvi-{rotation_id}-{seed}".encode())
        decoder.write_bytes(f"decoder-{rotation_id}-{seed}".encode())
        peakvi_sha = digest(peakvi)
        decoder_sha = digest(decoder)
        source_seed_records[rotation_id].append(
            {
                "seed": seed,
                "peakvi_state_file_sha256": peakvi_sha,
                "decoder_state_file_sha256": decoder_sha,
                "peakvi_parameter_sha256": f"peakvi-parameter-{rotation_id}-{seed}",
                "decoder_parameter_sha256": f"decoder-parameter-{rotation_id}-{seed}",
            }
        )
        local_seed_records[rotation_id].append(
            {
                "seed": seed,
                "source_peakvi_state_path": str(peakvi.resolve()),
                "source_peakvi_state_sha256": peakvi_sha,
                "source_decoder_state_path": str(decoder.resolve()),
                "source_decoder_state_sha256": decoder_sha,
                "source_peakvi_parameter_sha256": f"peakvi-parameter-{rotation_id}-{seed}",
                "source_decoder_parameter_sha256": f"decoder-parameter-{rotation_id}-{seed}",
                "strict_state_restore_passed": True,
                "target_encoding_repeat_bit_identical": True,
                "fit_or_refit_performed": False,
            }
        )
    write_json_exclusive(
        source_rotation / "fit_receipt.json",
        {
            "status": "passed",
            "rotation_id": rotation_id,
            "target_query_adaptation_performed": False,
            "metrics_calculated": False,
            "seeds": source_seed_records[rotation_id],
        },
    )
source_sha = freeze_tree(
    source,
    {
        "artifact_class": "peakvi_cross_cohort_fit_predict_and_blind_commit",
        "model_id": "peakvi_source_fit_fixed_window_decoder",
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "development_outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
commits = []
for rotation_id, context_role in rotations.items():
    prediction = stage / "predictions" / rotation_id
    prediction.mkdir(parents=True)
    values = np.full((18, 4, 16000), np.nan, dtype=np.float32)
    values[states == 0] = np.float32(0.5)
    missing = np.broadcast_to(states[:, :, None], values.shape).copy()
    predictions_path = prediction / "predictions.float32.npy"
    missing_path = prediction / "missing_state.uint8.npy"
    np.save(predictions_path, values, allow_pickle=False)
    np.save(missing_path, missing, allow_pickle=False)
    model_path = prediction / "model_execution_manifest.json"
    write_json_exclusive(
        model_path,
        {
            "schema_version": "masld-bench-peakvi-gse244832-model-execution-v1",
            "model_id": "peakvi",
            "rotation_id": rotation_id,
            "source_dataset_id": "gse296875",
            "execution_registration_artifacts_sha256": registration_sha,
            "target_dataset_id": "gse244832",
            "fit_scope": "gse296875_source_training_only",
            "source_fit_job_id": "21101210",
            "source_fit_result_path": str(source.resolve()),
            "source_fit_result_artifacts_sha256": source_sha,
            "sequence_features_consumed": False,
            "rna_assay_consumed": False,
            "target_dataset_fit_or_adaptation": False,
            "per_donor_context_operation": "not_applicable",
            "fit_or_refit_performed": False,
            "registered_output_family": "masked_accessibility_count",
            "prediction_scale": "mean_predicted_single_nucleus_binary_accessibility_probability",
            "prediction_is_observed_fragment_count": False,
            "seeds": local_seed_records[rotation_id],
            "metrics_calculated": False,
            "fixture_only": True,
        },
    )
    write_json_exclusive(
        prediction / "prediction_bundle.json",
        make_prediction_manifest(
            config=config,
            config_sha256="f" * 64,
            rotation_id=rotation_id,
            model_manifest=model_path,
            query_authority=queries[context_role],
            predictions=predictions_path,
            missing_state=missing_path,
            output=prediction,
        ),
    )
    prediction_sha = freeze_tree(
        prediction,
        {
            "artifact_class": "peakvi_gse244832_commit_integration_fixture_prediction",
            "rotation_id": rotation_id,
            "fixture_only": True,
            "outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    commit_root = stage / "commits" / rotation_id
    receipt = commit_peakvi(
        prediction_root=prediction,
        prediction_artifacts_sha256=prediction_sha,
        axis_root=axis,
        axis_artifacts_sha256=axis_sha,
        registration_root=registration,
        registration_artifacts_sha256=registration_sha,
        source_fit_root=source,
        source_fit_artifacts_sha256=source_sha,
        output=commit_root,
    )
    verify_frozen_tree(commit_root)
    commits.append(
        {
            "rotation_id": rotation_id,
            "prediction_artifacts_sha256": prediction_sha,
            "commit_artifacts_sha256": digest(commit_root / "ARTIFACTS.json"),
            "status": receipt["status"],
            "source_fit_artifacts_sha256": source_sha,
        }
    )
write_json_exclusive(
    stage / "receipt.json",
    {
        "schema_version": "masld-bench-peakvi-gse244832-commit-integration-fixture-v1",
        "status": "passed",
        "fixture_only": True,
        "commits": commits,
        "eligible_units": 48,
        "structurally_missing_units": 4,
        "below_qc_units": 20,
        "fit_or_refit_performed": False,
        "source_fit_artifacts_sha256": source_sha,
        "rna_assay_read": False,
        "cross_assay_join_consumed": False,
        "outcomes_read": False,
        "metrics_calculated": False,
    },
)
freeze_tree(
    stage,
    {
        "artifact_class": "peakvi_gse244832_commit_integration_fixture",
        "fixture_only": True,
        "prediction_commits": 2,
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${OUTPUT}"
printf 'PeakVI GSE244832 commit integration fixture: %s\n' "${OUTPUT}"
sha256sum "${OUTPUT}/ARTIFACTS.json"
