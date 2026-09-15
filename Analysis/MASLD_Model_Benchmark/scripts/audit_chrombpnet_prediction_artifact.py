#!/usr/bin/env python3
"""Admit frozen ChromBPNet predictions without reading benchmark outcomes."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
from typing import Any

import h5py
import numpy as np

try:
    import tomllib  # noqa: F401
except ModuleNotFoundError:
    import tomli

    sys.modules["tomllib"] = tomli

from masld_bench.artifacts import verify_frozen_tree


COUNT_FIELDS = (
    "window_id",
    "selection_hash",
    "contig",
    "output_start",
    "output_end",
    "genomic_fold",
    "role",
    "ccre_class",
    "forward_logcount",
    "reverse_logcount",
    "strand_averaged_mass",
)
MODEL_CONTRACTS = {
    "full_model": {
        "role": "assay_qc_score",
        "regional_mass_transform": "log1p_absolute",
        "regional_mass_semantics": "predicted_count_after_log1p_inverse",
    },
    "nobias_model": {
        "role": "primary_biological_score",
        "regional_mass_transform": "log_component",
        "regional_mass_semantics": (
            "positive_sequence_component_mass_from_exp_log_component"
        ),
    },
}


class AdmissionError(ValueError):
    """Raised when a frozen prediction file violates inclusion."""


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AdmissionError(f"{path.name} is not a JSON object")
    return value


def _decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values.tolist()
    ]


def _artifact_records(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    records = manifest.get("artifacts")
    if not isinstance(records, list):
        raise AdmissionError("artifact manifest records differ")
    by_path = {str(record.get("path")): record for record in records}
    if len(by_path) != len(records):
        raise AdmissionError("artifact manifest contains duplicate paths")
    return by_path


def _audit_prediction(
    root: Path,
    model_id: str,
    manifest_records: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], list[tuple[str, str, str]]]:
    contract = MODEL_CONTRACTS[model_id]
    prediction_root = root / "predictions" / model_id
    summary = _json(prediction_root / "summary.json")
    if (
        summary.get("schema_version")
        != "masld-bench-chrombpnet-ccre-predictions-v1"
        or summary.get("status") != "pass"
        or summary.get("windows") != 32_000
        or summary.get("role_counts") != {"test": 16_000, "valid": 16_000}
        or summary.get("regional_mass_transform")
        != contract["regional_mass_transform"]
        or summary.get("regional_mass_semantics")
        != contract["regional_mass_semantics"]
        or summary.get("observed_atac_input_exposed") is not False
        or summary.get("benchmark_metrics_calculated") is not False
    ):
        raise AdmissionError(f"{model_id} prediction summary differs")

    counts_path = prediction_root / "regional_counts.tsv"
    with counts_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != COUNT_FIELDS:
            raise AdmissionError(f"{model_id} count fields differ")
        rows = [dict(row) for row in reader]
    if len(rows) != 32_000:
        raise AdmissionError(f"{model_id} count census differs")
    inventory: list[tuple[str, str, str]] = []
    role_counts = {"valid": 0, "test": 0}
    for row in rows:
        role = row["role"]
        if role not in role_counts:
            raise AdmissionError(f"{model_id} contains an ineligible role")
        role_counts[role] += 1
        inventory.append((row["window_id"], row["selection_hash"], role))
        for field in (
            "forward_logcount",
            "reverse_logcount",
            "strand_averaged_mass",
        ):
            value = float(row[field])
            if not math.isfinite(value):
                raise AdmissionError(f"{model_id} contains a non-finite count field")
        if float(row["strand_averaged_mass"]) < 0:
            raise AdmissionError(f"{model_id} contains negative regional mass")
    if role_counts != {"valid": 16_000, "test": 16_000}:
        raise AdmissionError(f"{model_id} role census differs")
    if len({row[0] for row in inventory}) != len(inventory):
        raise AdmissionError(f"{model_id} window IDs are not unique")

    profile_path = prediction_root / "profile_probabilities.h5"
    minimum_probability = float("inf")
    maximum_probability = float("-inf")
    maximum_row_sum_error = 0.0
    with h5py.File(profile_path, "r") as handle:
        if set(handle.keys()) != {
            "profile_probability",
            "selection_hash",
            "window_id",
        }:
            raise AdmissionError(f"{model_id} HDF5 datasets differ")
        profiles = handle["profile_probability"]
        if profiles.shape != (32_000, 1_000) or profiles.dtype != np.dtype("float32"):
            raise AdmissionError(f"{model_id} profile tensor differs")
        h5_windows = _decode(handle["window_id"][:])
        h5_selection = _decode(handle["selection_hash"][:])
        if h5_windows != [row[0] for row in inventory] or h5_selection != [
            row[1] for row in inventory
        ]:
            raise AdmissionError(f"{model_id} HDF5 and TSV row IDs differ")
        for start in range(0, profiles.shape[0], 512):
            block = profiles[start : start + 512]
            if not np.isfinite(block).all():
                raise AdmissionError(f"{model_id} profile contains non-finite values")
            minimum_probability = min(minimum_probability, float(block.min()))
            maximum_probability = max(maximum_probability, float(block.max()))
            error = np.abs(block.sum(axis=1, dtype=np.float64) - 1.0)
            maximum_row_sum_error = max(maximum_row_sum_error, float(error.max()))
        if minimum_probability < -1.0e-7 or maximum_probability > 1.0 + 1.0e-7:
            raise AdmissionError(f"{model_id} profile is outside probability bounds")
        if maximum_row_sum_error > 1.0e-6:
            raise AdmissionError(f"{model_id} profile does not sum to one")

    relative_counts = f"predictions/{model_id}/regional_counts.tsv"
    relative_profiles = f"predictions/{model_id}/profile_probabilities.h5"
    for relative in (relative_counts, relative_profiles):
        if relative not in manifest_records:
            raise AdmissionError(f"manifest omits {relative}")
    report = {
        "model_view": model_id,
        "scientific_role": contract["role"],
        "prediction_units": 32_000,
        "role_counts": role_counts,
        "row_id": "window_id",
        "join_guard": "selection_hash",
        "profile_schema": {
            "path": relative_profiles,
            "sha256": manifest_records[relative_profiles]["sha256"],
            "dataset": "profile_probability",
            "shape": [32_000, 1_000],
            "dtype": "float32",
            "minimum": minimum_probability,
            "maximum": maximum_probability,
            "maximum_row_sum_error": maximum_row_sum_error,
            "semantics": "reverse_complement_strand_averaged_probability",
        },
        "count_schema": {
            "path": relative_counts,
            "sha256": manifest_records[relative_counts]["sha256"],
            "fields": list(COUNT_FIELDS),
            "regional_mass_transform": contract["regional_mass_transform"],
            "regional_mass_semantics": contract["regional_mass_semantics"],
        },
        "observed_atac_input_exposed": False,
        "benchmark_metrics_calculated": False,
    }
    return report, inventory


def audit(artifact: Path, expected_artifacts_sha256: str, output: Path) -> dict[str, Any]:
    root = artifact.resolve(strict=True)
    if output.exists() or output.is_symlink():
        raise AdmissionError("output already exists")
    if any("outcome" in part.lower() for part in root.parts):
        raise AdmissionError("outcome artifact entered admission")
    artifacts_path = root / "ARTIFACTS.json"
    if _sha256(artifacts_path) != expected_artifacts_sha256:
        raise AdmissionError("source ARTIFACTS hash differs")
    verified = verify_frozen_tree(root)
    metadata = verified["metadata"]
    if (
        metadata.get("artifact_class")
        != "chrombpnet_full_depth_training_and_predictions"
        or metadata.get("model_id") != "chrombpnet"
        or metadata.get("champion_eligible") is not False
        or metadata.get("test_outcomes_used") is not False
        or metadata.get("benchmark_metrics_calculated") is not False
    ):
        raise AdmissionError("source artifact metadata differs")
    validation = _json(root / "validation" / "summary.json")
    if (
        validation.get("status") != "pass"
        or validation.get("primary_biological_score") != "nobias_model"
        or validation.get("assay_qc_score") != "full_model"
        or validation.get("test_outcomes_used") is not False
        or validation.get("observed_held_donor_atac_used_for_inference") is not False
        or validation.get("benchmark_metrics_calculated") is not False
        or validation.get("champion_eligible") is not False
    ):
        raise AdmissionError("training validation firewall differs")

    manifest = _json(artifacts_path)
    records = _artifact_records(manifest)
    model_paths = {
        "full_model": "model/chrombpnet.h5",
        "nobias_model": "model/chrombpnet_nobias.h5",
    }
    for model_id, relative in model_paths.items():
        if records.get(relative, {}).get("sha256") != validation["model_sha256"].get(
            model_id
        ):
            raise AdmissionError(f"{model_id} checkpoint hash differs")

    reports: dict[str, Any] = {}
    inventories: dict[str, list[tuple[str, str, str]]] = {}
    for model_id in MODEL_CONTRACTS:
        reports[model_id], inventories[model_id] = _audit_prediction(
            root, model_id, records
        )
    if inventories["full_model"] != inventories["nobias_model"]:
        raise AdmissionError("full and no-bias prediction row inventories differ")

    result = {
        "schema_version": "masld-bench-chrombpnet-prediction-admission-v1",
        "status": "pass",
        "source_artifact": str(root),
        "source_artifacts_sha256": expected_artifacts_sha256,
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "split_id": "donor0_genomic0",
        "seed": 20260824,
        "adaptation_lane": "native_full_head",
        "primary_biological_score": "nobias_model",
        "assay_qc_score": "full_model",
        "prediction_views": reports,
        "admission_disposition": (
            "admitted_for_independent_development_evaluation_only"
        ),
        "champion_eligible": False,
        "champion_ineligibility_reasons": [
            "single_development_split",
            "single_seed",
            "development_outcomes_not_scored",
            "sealed_external_evaluation_not_performed",
            "source_bias_model_is_a_restricted_smoke_comparator",
        ],
        "observed_held_donor_atac_used_for_inference": False,
        "evaluator_outcomes_exposed": False,
        "benchmark_metrics_calculated": False,
    }
    output.mkdir(mode=0o750)
    (output / "admission.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--expected-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = audit(
        arguments.artifact,
        arguments.expected_artifacts_sha256,
        arguments.output,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
