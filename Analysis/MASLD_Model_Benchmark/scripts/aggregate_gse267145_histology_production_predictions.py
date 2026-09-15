#!/usr/bin/env python3
"""Assemble frozen GSE267145 fold-by-seed predictions without reading outcomes."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import ArtifactError, verify_frozen_tree, write_json_exclusive


OUTER_FOLDS = tuple(range(5))
SEEDS = (1701, 1709, 1721, 1723, 1733)
STAGE3 = ("NOR", "NAFL", "NASH")
FIBROSIS_GROUP3 = ("F0", "F1", "F2_3")
MODEL_IDS = (
    "training_stage_distribution",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_knn",
    "h3_variance_pca_elastic_net",
    "h3_variance_pca_linear_svm",
    "h3_variance_pca_nearest_centroid",
    "h3_variance_pca_knn",
    "block_pca_elastic_net",
    "calibrated_late_fusion",
)
PREDICTION_FIELDS = (
    "participant_id",
    "outer_fold",
    "probability_NOR",
    "probability_NAFL",
    "probability_NASH",
    "predicted_stage3",
    "predicted_nash_crn_component_sum",
    "predicted_fibrosis_cumulative_expected",
    "predicted_fibrosis_regression",
    "probability_fibrosis_F0",
    "probability_fibrosis_F1",
    "probability_fibrosis_F2_3",
    "predicted_fibrosis_group3",
)
NUMERIC_FIELDS = (
    "probability_NOR",
    "probability_NAFL",
    "probability_NASH",
    "predicted_nash_crn_component_sum",
    "predicted_fibrosis_cumulative_expected",
    "predicted_fibrosis_regression",
    "probability_fibrosis_F0",
    "probability_fibrosis_F1",
    "probability_fibrosis_F2_3",
)


class AggregationError(RuntimeError):
    """Raised when a production unit or prediction bundle does not meet its requirements."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AggregationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
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


def _load_axis(folds: Path) -> tuple[list[str], dict[str, int]]:
    fields, rows = read_tsv(folds / "participant_outer_folds.tsv")
    if fields != ("participant_id", "outer_fold") or len(rows) != 99:
        raise AggregationError("frozen participant/fold axis differs")
    participant_ids = [row["participant_id"] for row in rows]
    if len(participant_ids) != len(set(participant_ids)):
        raise AggregationError("participant/fold axis contains duplicate participants")
    assignments = {row["participant_id"]: int(row["outer_fold"]) for row in rows}
    if set(assignments.values()) != set(OUTER_FOLDS):
        raise AggregationError("participant/fold axis does not contain five outer folds")
    return participant_ids, assignments


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AggregationError(f"invalid JSON document: {path}: {error}") from error
    if not isinstance(value, dict):
        raise AggregationError(f"JSON document is not an object: {path}")
    return value


def _verify_tree(path: Path) -> dict[str, Any]:
    try:
        return verify_frozen_tree(path)
    except ArtifactError as error:
        raise AggregationError(f"invalid frozen production tree: {path}: {error}") from error


def _validate_row(row: Mapping[str, str], *, path: Path) -> None:
    try:
        values = [float(row[field]) for field in NUMERIC_FIELDS]
    except (KeyError, ValueError) as error:
        raise AggregationError(f"prediction has invalid numeric value: {path}") from error
    if not all(math.isfinite(value) for value in values):
        raise AggregationError(f"prediction has non-finite numeric value: {path}")
    stage = [float(row[f"probability_{label}"]) for label in STAGE3]
    fibrosis = [float(row[f"probability_fibrosis_{label}"]) for label in FIBROSIS_GROUP3]
    if (
        any(value < 0.0 or value > 1.0 for value in stage + fibrosis)
        or not math.isclose(math.fsum(stage), 1.0, abs_tol=1e-6)
        or not math.isclose(math.fsum(fibrosis), 1.0, abs_tol=1e-6)
    ):
        raise AggregationError(f"prediction probabilities are invalid: {path}")
    if row["predicted_stage3"] != STAGE3[max(range(3), key=stage.__getitem__)]:
        raise AggregationError(f"stage label does not match probabilities: {path}")
    if row["predicted_fibrosis_group3"] != FIBROSIS_GROUP3[
        max(range(3), key=fibrosis.__getitem__)
    ]:
        raise AggregationError(f"fibrosis label does not match probabilities: {path}")


def _load_unit(
    unit: Path,
    *,
    outer_fold: int,
    seed: int,
    assignments: Mapping[str, int],
) -> tuple[dict[str, dict[str, dict[str, str]]], dict[str, str]]:
    unit_manifest = _verify_tree(unit)
    metadata = unit_manifest["metadata"]
    if (
        metadata.get("artifact_class")
        != "gse267145_histology_production_outer_seed_unit"
        or metadata.get("outer_fold") != outer_fold
        or metadata.get("seed") != seed
        or metadata.get("status") != "passed_unscored"
        or metadata.get("outcomes_read") is not False
        or metadata.get("metrics_calculated") is not False
    ):
        raise AggregationError(f"production unit metadata differs: {unit}")
    receipt = _read_json_object(unit / "unit_receipt.json")
    expected_receipt = {
        "outer_fold": outer_fold,
        "seed": seed,
        "model_family_fits": len(MODEL_IDS),
        "prediction_files": len(MODEL_IDS),
        "outer_test_outcomes_read": False,
        "outer_test_metrics_calculated": False,
        "production_unit_complete": True,
    }
    for field, expected in expected_receipt.items():
        if receipt.get(field) != expected:
            raise AggregationError(f"production unit receipt field differs: {unit}: {field}")
    predictions = unit / "fit/predictions"
    prediction_manifest = _verify_tree(predictions)
    prediction_metadata = prediction_manifest["metadata"]
    if (
        prediction_metadata.get("outer_fold") != outer_fold
        or prediction_metadata.get("seed") != seed
        or prediction_metadata.get("prediction_files") != len(MODEL_IDS)
        or prediction_metadata.get("outcomes_read") is not False
        or prediction_metadata.get("metrics_calculated") is not False
    ):
        raise AggregationError(f"unit prediction metadata differs: {predictions}")
    observed_files = sorted(path.name for path in predictions.glob("*.tsv"))
    expected_files = sorted(f"{model_id}.tsv" for model_id in MODEL_IDS)
    if observed_files != expected_files:
        raise AggregationError(f"unit prediction-file roster differs: {predictions}")

    models: dict[str, dict[str, dict[str, str]]] = {}
    expected_participants = {
        participant for participant, fold in assignments.items() if fold == outer_fold
    }
    for model_id in MODEL_IDS:
        path = predictions / f"{model_id}.tsv"
        fields, rows = read_tsv(path)
        if fields != PREDICTION_FIELDS:
            raise AggregationError(f"prediction schema differs: {path}")
        indexed: dict[str, dict[str, str]] = {}
        for row in rows:
            participant = row["participant_id"]
            if participant in indexed:
                raise AggregationError(f"duplicate participant prediction: {path}")
            if assignments.get(participant) != outer_fold or int(row["outer_fold"]) != outer_fold:
                raise AggregationError(f"prediction participant is outside held fold: {path}")
            _validate_row(row, path=path)
            indexed[participant] = row
        if set(indexed) != expected_participants:
            raise AggregationError(f"held-fold participant roster differs: {path}")
        models[model_id] = indexed
    hashes = {
        "unit_artifacts_sha256": sha256_file(unit / "ARTIFACTS.json"),
        "predictions_artifacts_sha256": sha256_file(predictions / "ARTIFACTS.json"),
    }
    return models, hashes


def _ensemble_rows(seed_rows: Sequence[Mapping[str, str]]) -> dict[str, str]:
    if len(seed_rows) != len(SEEDS):
        raise AggregationError("seed ensemble has incomplete membership")
    first = seed_rows[0]
    if any(
        row["participant_id"] != first["participant_id"]
        or row["outer_fold"] != first["outer_fold"]
        for row in seed_rows[1:]
    ):
        raise AggregationError("seed ensemble participant axes differ")
    averaged = {
        field: math.fsum(float(row[field]) for row in seed_rows) / len(seed_rows)
        for field in NUMERIC_FIELDS
    }
    stage = [averaged[f"probability_{label}"] for label in STAGE3]
    fibrosis = [averaged[f"probability_fibrosis_{label}"] for label in FIBROSIS_GROUP3]
    return {
        "participant_id": first["participant_id"],
        "outer_fold": first["outer_fold"],
        **{field: format(averaged[field], ".17g") for field in NUMERIC_FIELDS[:3]},
        "predicted_stage3": STAGE3[max(range(3), key=stage.__getitem__)],
        **{field: format(averaged[field], ".17g") for field in NUMERIC_FIELDS[3:6]},
        **{field: format(averaged[field], ".17g") for field in NUMERIC_FIELDS[6:]},
        "predicted_fibrosis_group3": FIBROSIS_GROUP3[
            max(range(3), key=fibrosis.__getitem__)
        ],
    }


def aggregate(*, units: Path, folds: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise AggregationError(f"refusing to overwrite aggregation output: {output}")
    participant_ids, assignments = _load_axis(folds)
    output.mkdir(parents=True)
    prediction_root = output / "predictions"
    prediction_root.mkdir()
    assembled: dict[tuple[str, int], dict[str, dict[str, str]]] = {
        (model, seed): {} for model in MODEL_IDS for seed in SEEDS
    }
    unit_manifest_rows: list[dict[str, Any]] = []
    for outer_fold in OUTER_FOLDS:
        for seed in SEEDS:
            unit = units / f"outer_{outer_fold}" / f"seed_{seed}"
            models, hashes = _load_unit(
                unit,
                outer_fold=outer_fold,
                seed=seed,
                assignments=assignments,
            )
            for model_id, rows in models.items():
                overlap = set(assembled[(model_id, seed)]).intersection(rows)
                if overlap:
                    raise AggregationError("participant appears in multiple outer folds")
                assembled[(model_id, seed)].update(rows)
            unit_manifest_rows.append(
                {
                    "outer_fold": outer_fold,
                    "seed": seed,
                    "unit_path": unit.as_posix(),
                    **hashes,
                    "status": "passed_unscored",
                }
            )

    for model_id in MODEL_IDS:
        for seed in SEEDS:
            rows = assembled[(model_id, seed)]
            if set(rows) != set(participant_ids):
                raise AggregationError("assembled participant roster is incomplete")
            write_tsv(
                prediction_root / f"{model_id}--seed-{seed}.tsv",
                PREDICTION_FIELDS,
                (rows[participant] for participant in participant_ids),
            )
        write_tsv(
            prediction_root / f"{model_id}.tsv",
            PREDICTION_FIELDS,
            (
                _ensemble_rows(
                    [assembled[(model_id, seed)][participant] for seed in SEEDS]
                )
                for participant in participant_ids
            ),
        )

    write_tsv(
        output / "unit_manifest.tsv",
        (
            "outer_fold",
            "seed",
            "unit_path",
            "unit_artifacts_sha256",
            "predictions_artifacts_sha256",
            "status",
        ),
        unit_manifest_rows,
    )
    receipt = {
        "schema_version": "masld-bench-gse267145-production-aggregation-v1",
        "status": "passed_predictions_unscored",
        "participants": len(participant_ids),
        "outer_folds": len(OUTER_FOLDS),
        "model_seeds": list(SEEDS),
        "model_ids": list(MODEL_IDS),
        "logical_outer_seed_units": len(OUTER_FOLDS) * len(SEEDS),
        "model_family_unit_fits": len(OUTER_FOLDS) * len(SEEDS) * len(MODEL_IDS),
        "seed_prediction_files": len(MODEL_IDS) * len(SEEDS),
        "ensemble_prediction_files": len(MODEL_IDS),
        "prediction_files": len(MODEL_IDS) * (len(SEEDS) + 1),
        "outcomes_read": False,
        "metrics_calculated": False,
        "scorer_called": False,
        "seeds_are_biological_replicates": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--units", required=True, type=Path)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--folds-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    _verify_tree(arguments.folds)
    if sha256_file(arguments.folds / "ARTIFACTS.json") != arguments.folds_artifacts_sha256:
        raise AggregationError("fold artifact manifest SHA-256 differs")
    receipt = aggregate(units=arguments.units, folds=arguments.folds, output=arguments.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
