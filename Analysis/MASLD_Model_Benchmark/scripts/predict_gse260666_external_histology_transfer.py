#!/usr/bin/env python3
"""Predict GSE260666 from frozen source models without label access or query fit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.contracts import ArtifactRef, MissingState, PredictionBundle
from masld_bench.hashing import canonical_sha256


CLASSES = ("NOR", "NAFL", "NASH")
MODEL_IDS = (
    "training_stage_distribution",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
)
TASK_ID = "gse260666_external_histology_transfer"
DATASET_IDS = ("gse260666_bulk_rna",)
SPLIT_ID = "gse267145_full_development_refit_to_gse260666_external_v1"
UNIT_NAMESPACE = "gse260666_outcome_blind_row_id"
PREDICTION_FIELDS = (
    "prediction_row_id",
    "row_id",
    "endpoint_id",
    "prediction_state",
    "model_id",
    "probability_NOR",
    "probability_NAFL",
    "probability_NASH",
    "predicted_stage3",
)


class ExternalPredictionError(RuntimeError):
    """Raised when label-blind external prediction would violate its selection record."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExternalPredictionError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def softmax(logits: np.ndarray) -> np.ndarray:
    values = np.asarray(logits, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)):
        raise ExternalPredictionError("prediction logits are invalid")
    shifted = values - values.max(axis=1, keepdims=True)
    exponential = np.exp(shifted)
    probabilities = exponential / exponential.sum(axis=1, keepdims=True)
    if not np.allclose(probabilities.sum(axis=1), 1.0, rtol=0.0, atol=1e-12):
        raise ExternalPredictionError("prediction probabilities do not sum to one")
    return probabilities


def log1p_cpm(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ExternalPredictionError("external RNA matrix is not finite and nonnegative")
    totals = values.sum(axis=1)
    if np.any(totals <= 0):
        raise ExternalPredictionError("external participant RNA library is empty")
    return np.log1p(values * (1_000_000.0 / totals[:, None]))


def infer_seed(
    selected_log_cpm: np.ndarray,
    state: Mapping[str, np.ndarray],
    *,
    centroid_temperature: float = 1.0,
) -> dict[str, np.ndarray]:
    required = {
        "center",
        "scale",
        "pca_mean",
        "pca_components",
        "svm_coef",
        "svm_intercept",
        "svm_calibrator_coef",
        "svm_calibrator_intercept",
        "elastic_coef",
        "elastic_intercept",
        "centroid_centroids",
    }
    if set(state) != required or centroid_temperature <= 0:
        raise ExternalPredictionError("source model-state schema differs")
    values = np.asarray(selected_log_cpm, dtype=np.float64)
    center = np.asarray(state["center"], dtype=np.float64)
    scale = np.asarray(state["scale"], dtype=np.float64)
    components = np.asarray(state["pca_components"], dtype=np.float64)
    pca_mean = np.asarray(state["pca_mean"], dtype=np.float64)
    if (
        values.ndim != 2
        or center.shape != (values.shape[1],)
        or scale.shape != center.shape
        or components.shape[1] != values.shape[1]
        or pca_mean.shape != center.shape
        or np.any(scale <= 0)
    ):
        raise ExternalPredictionError("source preprocessing state differs")
    scaled = (values - center) / scale
    projection = (scaled - pca_mean) @ components.T
    elastic_logits = (
        projection @ np.asarray(state["elastic_coef"], dtype=np.float64).T
        + np.asarray(state["elastic_intercept"], dtype=np.float64)
    )
    svm_raw = (
        projection @ np.asarray(state["svm_coef"], dtype=np.float64).T
        + np.asarray(state["svm_intercept"], dtype=np.float64)
    )
    svm_logits = (
        svm_raw @ np.asarray(state["svm_calibrator_coef"], dtype=np.float64).T
        + np.asarray(state["svm_calibrator_intercept"], dtype=np.float64)
    )
    centroids = np.asarray(state["centroid_centroids"], dtype=np.float64)
    distances = np.sum(
        (projection[:, None, :] - centroids[None, :, :]) ** 2, axis=2
    )
    return {
        "rna_hvg_pca_elastic_net": softmax(elastic_logits),
        "rna_hvg_pca_linear_svm": softmax(svm_logits),
        "rna_hvg_pca_nearest_centroid": softmax(
            -distances / centroid_temperature
        ),
    }


def _prediction_row_id(run_id: str, model_id: str, row_id: str) -> str:
    return hashlib.sha256(
        f"{run_id}\0{model_id}\0stage3\0{row_id}".encode("utf-8")
    ).hexdigest()


def _freeze_bundle(
    *,
    output: Path,
    model_id: str,
    row_ids: Sequence[str],
    probabilities: np.ndarray,
    run_id: str,
    source_fit_artifacts_sha256: str,
    external_model_input_artifacts_sha256: str,
    model_seeds: Sequence[int],
) -> dict[str, Any]:
    target = output / model_id
    target.mkdir()
    rows: list[dict[str, Any]] = []
    row_id_rows: list[dict[str, str]] = []
    for row_id, probability in zip(row_ids, probabilities, strict=True):
        prediction_row_id = _prediction_row_id(run_id, model_id, row_id)
        rows.append(
            {
                "prediction_row_id": prediction_row_id,
                "row_id": row_id,
                "endpoint_id": "stage3",
                "prediction_state": "observed",
                "model_id": model_id,
                "probability_NOR": format(float(probability[0]), ".17g"),
                "probability_NAFL": format(float(probability[1]), ".17g"),
                "probability_NASH": format(float(probability[2]), ".17g"),
                "predicted_stage3": CLASSES[int(np.argmax(probability))],
            }
        )
        row_id_rows.append(
            {"prediction_row_id": prediction_row_id, "row_id": row_id}
        )
    table = target / "predictions.tsv"
    row_inventory = target / "row_ids.tsv"
    write_tsv(table, PREDICTION_FIELDS, rows)
    write_tsv(
        row_inventory,
        ("prediction_row_id", "row_id"),
        row_id_rows,
    )
    table_ref = ArtifactRef.from_path(
        table,
        relative_to=target,
        media_type="text/tab-separated-values",
        role=f"standardized_prediction_table:{TASK_ID}",
    )
    row_ref = ArtifactRef.from_path(
        row_inventory,
        relative_to=target,
        media_type="text/tab-separated-values",
        role=f"prediction_row_ids:{TASK_ID}",
    )
    join_key = canonical_sha256(
        {
            "task_id": TASK_ID,
            "dataset_ids": list(DATASET_IDS),
            "split_id": SPLIT_ID,
            "row_id_field": "prediction_row_id",
            "unit_id_field": "row_id",
            "unit_id_namespace": UNIT_NAMESPACE,
            "biological_unit": "participant",
        }
    )
    identity = {
        "run_id": run_id,
        "task_id": TASK_ID,
        "model_id": model_id,
        "dataset_ids": list(DATASET_IDS),
        "split_id": SPLIT_ID,
        "table_sha256": table_ref.sha256,
        "row_ids_sha256": row_ref.sha256,
        "n_predictions": len(row_ids),
        "source_join_key_sha256": join_key,
    }
    bundle = PredictionBundle(
        schema_version="masld-bench-prediction-bundle-v1",
        bundle_id=canonical_sha256(identity),
        run_id=run_id,
        task_id=TASK_ID,
        model_id=model_id,
        dataset_ids=DATASET_IDS,
        split_id=SPLIT_ID,
        artifacts=(table_ref, row_ref),
        standardized_table=table_ref,
        row_ids=row_ref,
        n_predictions=len(row_ids),
        row_id_field="prediction_row_id",
        unit_id_field="row_id",
        unit_id_namespace=UNIT_NAMESPACE,
        biological_unit="participant",
        table_schema_sha256=canonical_sha256(
            {"format": "tsv", "fields": list(PREDICTION_FIELDS)}
        ),
        source_join_key_sha256=join_key,
        format_version="gse260666-external-stage3-prediction-tsv-v1",
        missing_state=MissingState.OBSERVED,
        metadata={
            "endpoint_id": "stage3",
            "prediction_role": "five_seed_source_full_fit_external_development",
            "prediction_aggregation": "mean_probabilities_across_five_source_model_seeds",
            "model_seed_roster": list(model_seeds),
            "source_fit_artifacts_sha256": source_fit_artifacts_sha256,
            "external_model_input_artifacts_sha256": external_model_input_artifacts_sha256,
            "external_labels_read": False,
            "query_fit_or_calibration_performed": False,
            "prediction_frozen_before_evaluator_label_join": True,
            "project_sealed": False,
            "champion_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "clinical_claim_allowed": False,
        },
    )
    document = target / "prediction_bundle.json"
    write_json_exclusive(document, bundle.to_dict(), mode=0o440)
    loaded = PredictionBundle.load_json(document)
    loaded.validate_artifacts(target)
    artifacts_sha256 = freeze_tree(
        target,
        {
            "artifact_class": "gse260666_external_histology_prediction_bundle",
            "bundle_id": bundle.bundle_id,
            "model_id": model_id,
            "participants": len(row_ids),
            "biological_unit": "participant",
            "external_labels_read": False,
            "query_fit_or_calibration_performed": False,
            "metrics_calculated": False,
            "project_sealed": False,
            "status": "passed_unscored",
        },
    )
    return {
        "model_id": model_id,
        "bundle_path": f"{model_id}/prediction_bundle.json",
        "bundle_id": bundle.bundle_id,
        "bundle_document_sha256": sha256_file(document),
        "bundle_artifacts_sha256": artifacts_sha256,
        "prediction_table_sha256": table_ref.sha256,
    }


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise ExternalPredictionError(f"{label} SHA-256 differs")


def predict_external(
    *, benchmark_root: Path, contract_path: Path, contract_sha256: str, output: Path
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise ExternalPredictionError(f"refusing to overwrite predictions: {output}")
    _check_hash(contract_path, contract_sha256, "external prediction contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_prediction_only_pending"
        or contract.get("external_labels_read") is not False
        or contract.get("query_fit_or_calibration_allowed") is not False
    ):
        raise ExternalPredictionError("external prediction contract differs")
    source_spec = contract["source_fit"]
    source_fit = benchmark_root / source_spec["path"]
    _check_hash(
        source_fit / "ARTIFACTS.json",
        source_spec["artifacts_sha256"],
        "source full fit",
    )
    verify_frozen_tree(source_fit)
    fit_receipt = json.loads(
        (source_fit / "fit_receipt.json").read_text(encoding="utf-8")
    )
    if (
        fit_receipt.get("status") != source_spec["required_status"]
        or fit_receipt.get("external_labels_read") is not False
        or fit_receipt.get("external_expression_values_read") is not False
    ):
        raise ExternalPredictionError("source full-fit receipt differs")
    external_spec = contract["external_model_input"]
    model_input = benchmark_root / external_spec["path"]
    _check_hash(
        model_input / "ARTIFACTS.json",
        external_spec["artifacts_sha256"],
        "external model input",
    )
    verify_frozen_tree(model_input)
    model_receipt = json.loads(
        (model_input / "receipt.json").read_text(encoding="utf-8")
    )
    if (
        model_receipt.get("labels_included") is not False
        or model_receipt.get("query_preprocessing_fitted") is not False
        or model_receipt.get("participants") != external_spec["participants"]
    ):
        raise ExternalPredictionError("external model-input receipt differs")
    participant_fields, participants = read_tsv(model_input / "participant_axis.tsv")
    feature_fields, external_features = read_tsv(model_input / "rna_feature_axis.tsv")
    selected_fields, selected_features = read_tsv(source_fit / "selected_feature_axis.tsv")
    if (
        "row_id" not in participant_fields
        or "stable_gene_id" not in feature_fields
        or "stable_gene_id" not in selected_fields
        or len(participants) != external_spec["participants"]
    ):
        raise ExternalPredictionError("external or selected feature axis differs")
    row_ids = [row["row_id"] for row in participants]
    if len(row_ids) != len(set(row_ids)):
        raise ExternalPredictionError("external row ID is duplicated")
    external_ids = [row["stable_gene_id"] for row in external_features]
    selected_ids = [row["stable_gene_id"] for row in selected_features]
    external_lookup = {value: index for index, value in enumerate(external_ids)}
    if (
        len(external_lookup) != len(external_ids)
        or len(selected_ids) != source_spec["selected_features"]
        or any(value not in external_lookup for value in selected_ids)
    ):
        raise ExternalPredictionError("selected source gene is not observable externally")
    counts = np.load(model_input / "rna_counts.npy", mmap_mode="r", allow_pickle=False)
    if counts.shape != (len(row_ids), len(external_ids)):
        raise ExternalPredictionError("external RNA count matrix shape differs")
    transformed = log1p_cpm(counts)
    selected_values = transformed[:, [external_lookup[value] for value in selected_ids]]
    seeds = [int(value) for value in source_spec["model_seeds"]]
    by_model: dict[str, list[np.ndarray]] = {
        model_id: [] for model_id in MODEL_IDS if model_id != "training_stage_distribution"
    }
    for seed in seeds:
        seed_root = source_fit / f"seed_models/seed_{seed}"
        verify_frozen_tree(seed_root)
        with np.load(seed_root / "model_state.npz", allow_pickle=False) as loaded:
            state = {key: np.asarray(loaded[key]) for key in loaded.files}
        predictions = infer_seed(
            selected_values,
            state,
            centroid_temperature=float(contract["centroid_temperature"]),
        )
        for model_id, values in predictions.items():
            by_model[model_id].append(values)
    probabilities: dict[str, np.ndarray] = {
        model_id: np.mean(np.stack(values, axis=0), axis=0)
        for model_id, values in by_model.items()
    }
    source_priors = fit_receipt["training_stage_distribution"]
    prior = np.asarray([float(source_priors[label]) for label in CLASSES])
    probabilities["training_stage_distribution"] = np.repeat(
        prior[None, :], len(row_ids), axis=0
    )
    if set(probabilities) != set(MODEL_IDS):
        raise ExternalPredictionError("external prediction model roster differs")
    for values in probabilities.values():
        if (
            values.shape != (len(row_ids), 3)
            or not np.all(np.isfinite(values))
            or np.any(values < 0)
            or not np.allclose(values.sum(axis=1), 1.0, rtol=0.0, atol=1e-8)
        ):
            raise ExternalPredictionError("external ensemble probabilities are invalid")
    output.mkdir(parents=True)
    index_rows = [
        _freeze_bundle(
            output=output,
            model_id=model_id,
            row_ids=row_ids,
            probabilities=probabilities[model_id],
            run_id=source_spec["artifacts_sha256"],
            source_fit_artifacts_sha256=source_spec["artifacts_sha256"],
            external_model_input_artifacts_sha256=external_spec["artifacts_sha256"],
            model_seeds=seeds,
        )
        for model_id in MODEL_IDS
    ]
    write_tsv(
        output / "prediction_bundle_index.tsv",
        tuple(index_rows[0]),
        index_rows,
    )
    receipt = {
        "schema_version": "masld-bench-gse260666-external-prediction-only-v1",
        "status": "passed_frozen_prediction_only_no_label_access",
        "participants": len(row_ids),
        "biological_unit": "participant",
        "models": list(MODEL_IDS),
        "prediction_bundles": len(index_rows),
        "model_seeds": seeds,
        "prediction_aggregation": "mean_probabilities_across_five_source_model_seeds",
        "seed_level_results_emitted": False,
        "external_labels_read": False,
        "query_fit_or_calibration_performed": False,
        "query_participants_jointly_normalized": False,
        "external_metrics_calculated": False,
        "prediction_frozen_before_evaluator_label_join": True,
        "project_sealed": False,
        "champion_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
    }
    with (output / "prediction_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = predict_external(
        benchmark_root=arguments.benchmark_root,
        contract_path=arguments.contract,
        contract_sha256=arguments.contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
