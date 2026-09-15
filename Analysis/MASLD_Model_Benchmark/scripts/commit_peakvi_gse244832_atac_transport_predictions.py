#!/usr/bin/env python3
"""Commit GSE244832 PeakVI predictions only after frozen source-state proof."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.commit_gse244832_atac_transport_predictions import (
    ATACPredictionCommitError,
    _digest,
    _load_json,
    resolve_model_manifest,
    validate_bundle,
)


SEEDS = (1103, 2207, 3301)
SOURCE_ARTIFACT_CLASS = "peakvi_cross_cohort_fit_predict_and_blind_commit"
SOURCE_MODEL_ID = "peakvi_source_fit_fixed_window_decoder"


def validate_peakvi_source_reuse(
    prediction_root: Path,
    source_fit_root: Path,
    source_fit_artifacts_sha256: str,
) -> dict[str, Any]:
    prediction_root = prediction_root.resolve(strict=True)
    source_fit_root = source_fit_root.resolve(strict=True)
    source_artifacts = verify_frozen_tree(source_fit_root)
    source_metadata = source_artifacts.get("metadata", {})
    if (
        source_fit_root.is_symlink()
        or _digest(source_fit_root / "ARTIFACTS.json")
        != source_fit_artifacts_sha256
        or source_metadata.get("artifact_class") != SOURCE_ARTIFACT_CLASS
        or source_metadata.get("model_id") != SOURCE_MODEL_ID
        or source_metadata.get("source_dataset_id") != "gse296875"
        or source_metadata.get("target_dataset_id") != "gse281367"
        or source_metadata.get("development_outcomes_read") is not False
        or source_metadata.get("metrics_calculated") is not False
    ):
        raise ATACPredictionCommitError("PeakVI source-fit authority differs")

    prediction_manifest = _load_json(prediction_root / "prediction_bundle.json")
    model_record = prediction_manifest.get("model_manifest")
    if not isinstance(model_record, dict):
        raise ATACPredictionCommitError("PeakVI model execution manifest is absent")
    model = _load_json(resolve_model_manifest(prediction_root, model_record))
    rotation_id = str(prediction_manifest.get("rotation_id", ""))
    seed_records = model.get("seeds")
    if (
        model.get("schema_version")
        != "masld-bench-peakvi-gse244832-model-execution-v1"
        or model.get("model_id") != "peakvi"
        or model.get("rotation_id") != rotation_id
        or model.get("source_dataset_id") != "gse296875"
        or model.get("target_dataset_id") != "gse244832"
        or model.get("source_fit_job_id") != "21101210"
        or Path(str(model.get("source_fit_result_path", ""))).resolve(strict=True)
        != source_fit_root
        or model.get("source_fit_result_artifacts_sha256")
        != source_fit_artifacts_sha256
        or model.get("fit_or_refit_performed") is not False
        or model.get("target_dataset_fit_or_adaptation") is not False
        or model.get("registered_output_family") != "masked_accessibility_count"
        or model.get("prediction_scale")
        != "mean_predicted_single_nucleus_binary_accessibility_probability"
        or model.get("prediction_is_observed_fragment_count") is not False
        or not isinstance(seed_records, list)
        or [record.get("seed") for record in seed_records] != list(SEEDS)
    ):
        raise ATACPredictionCommitError("PeakVI no-refit execution contract differs")

    source_rotation = source_fit_root / "predictions" / rotation_id
    source_receipt = _load_json(source_rotation / "fit_receipt.json")
    source_seeds = source_receipt.get("seeds")
    if (
        source_receipt.get("status") != "passed"
        or source_receipt.get("rotation_id") != rotation_id
        or source_receipt.get("target_query_adaptation_performed") is not False
        or source_receipt.get("metrics_calculated") is not False
        or not isinstance(source_seeds, list)
        or [record.get("seed") for record in source_seeds] != list(SEEDS)
    ):
        raise ATACPredictionCommitError("PeakVI source-fit receipt differs")

    state_hashes: list[dict[str, Any]] = []
    for local, source in zip(seed_records, source_seeds, strict=True):
        seed = int(local["seed"])
        state_root = source_rotation / f"seed-{seed}"
        peakvi_state = (state_root / "peakvi.state_dict.pt").resolve(strict=True)
        decoder_state = (state_root / "decoder.state_dict.pt").resolve(strict=True)
        for state_path in (peakvi_state, decoder_state):
            try:
                state_path.relative_to(source_fit_root)
            except ValueError as error:
                raise ATACPredictionCommitError(
                    "PeakVI source state escaped frozen authority"
                ) from error
            if state_path.is_symlink() or not state_path.is_file():
                raise ATACPredictionCommitError("PeakVI source state differs")
        peakvi_sha = _digest(peakvi_state)
        decoder_sha = _digest(decoder_state)
        if (
            str(peakvi_state) != local.get("source_peakvi_state_path")
            or str(decoder_state) != local.get("source_decoder_state_path")
            or peakvi_sha != local.get("source_peakvi_state_sha256")
            or decoder_sha != local.get("source_decoder_state_sha256")
            or peakvi_sha != source.get("peakvi_state_file_sha256")
            or decoder_sha != source.get("decoder_state_file_sha256")
            or local.get("source_peakvi_parameter_sha256")
            != source.get("peakvi_parameter_sha256")
            or local.get("source_decoder_parameter_sha256")
            != source.get("decoder_parameter_sha256")
            or local.get("strict_state_restore_passed") is not True
            or local.get("target_encoding_repeat_bit_identical") is not True
            or local.get("fit_or_refit_performed") is not False
        ):
            raise ATACPredictionCommitError("PeakVI source-state reuse proof differs")
        state_hashes.append(
            {
                "seed": seed,
                "peakvi_state_sha256": peakvi_sha,
                "decoder_state_sha256": decoder_sha,
            }
        )
    return {
        "source_fit_job_id": "21101210",
        "source_fit_artifacts_sha256": source_fit_artifacts_sha256,
        "rotation_id": rotation_id,
        "state_hashes": state_hashes,
        "strict_state_restore_passed": True,
        "fit_or_refit_performed": False,
    }


def commit_peakvi(
    prediction_root: Path,
    prediction_artifacts_sha256: str,
    axis_root: Path,
    axis_artifacts_sha256: str,
    registration_root: Path,
    registration_artifacts_sha256: str,
    source_fit_root: Path,
    source_fit_artifacts_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise ATACPredictionCommitError("refusing to overwrite PeakVI prediction commit")
    source_reuse = validate_peakvi_source_reuse(
        prediction_root,
        source_fit_root,
        source_fit_artifacts_sha256,
    )
    receipt = validate_bundle(
        prediction_root,
        prediction_artifacts_sha256,
        axis_root,
        axis_artifacts_sha256,
        registration_root,
        registration_artifacts_sha256,
    )
    if receipt.get("model_id") != "peakvi":
        raise ATACPredictionCommitError("PeakVI-specific commit received another model")
    receipt["peakvi_source_reuse"] = source_reuse
    output.mkdir(parents=True, mode=0o750)
    write_json_exclusive(output / "commit.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse244832_atac_transport_prediction_commit",
            "bundle_id": receipt["bundle_id"],
            "model_id": "peakvi",
            "prediction_artifacts_sha256": prediction_artifacts_sha256,
            "axis_artifacts_sha256": axis_artifacts_sha256,
            "execution_registration_artifacts_sha256": registration_artifacts_sha256,
            "source_fit_artifacts_sha256": source_fit_artifacts_sha256,
            "fit_or_refit_performed": False,
            "development_outcomes_read": False,
            "metrics_calculated": False,
            "sequence_execution_authorized": False,
            "status": "passed",
        },
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-root", required=True, type=Path)
    parser.add_argument("--prediction-artifacts-sha256", required=True)
    parser.add_argument("--axis-root", required=True, type=Path)
    parser.add_argument("--axis-artifacts-sha256", required=True)
    parser.add_argument("--registration-root", required=True, type=Path)
    parser.add_argument("--registration-artifacts-sha256", required=True)
    parser.add_argument("--source-fit-root", required=True, type=Path)
    parser.add_argument("--source-fit-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = commit_peakvi(**vars(arguments))
    import json

    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
