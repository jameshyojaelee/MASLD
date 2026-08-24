#!/usr/bin/env python
"""Independently evaluate one donor-held-out cell-state smoke bundle."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
from typing import Any, Iterable, Mapping, Sequence


class CellEvaluationError(RuntimeError):
    """Raised when a prediction or outcome join violates the frozen contract."""


TASK_ID = "cell_state_mapping"
DATASET_ID = "resource_atlas_current"
EVALUATOR_ID = "cell_donor_balanced_macro_f1_v1"
_HEX = frozenset("0123456789abcdef")


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _canonical_hash(value: Any) -> str:
    return sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_sha256(value: Any, label: str) -> str:
    text = str(value)
    if len(text) != 64 or any(character not in _HEX for character in text):
        raise CellEvaluationError(f"{label} must be a lowercase SHA-256")
    return text


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(_canonical_json(value))
        handle.write("\n")


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CellEvaluationError(f"TSV has no header: {path}")
        fields = tuple(reader.fieldnames)
        rows = [dict(row) for row in reader]
    return fields, rows


def _fold_index(unit_id: str, *, seed: int, outer_folds: int) -> int:
    digest = sha256(f"{seed}\0{unit_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % outer_folds


def _join_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode("utf-8")).hexdigest()


def donor_class_balanced_weights(
    donors: Sequence[str], labels: Sequence[str]
) -> list[float]:
    if not donors or len(donors) != len(labels):
        raise CellEvaluationError("donors and labels must be non-empty and aligned")
    counts: dict[tuple[str, str], int] = {}
    donors_by_class: dict[str, set[str]] = {}
    for donor, label in zip(donors, labels):
        counts[(donor, label)] = counts.get((donor, label), 0) + 1
        donors_by_class.setdefault(label, set()).add(donor)
    n_classes = len(donors_by_class)
    return [
        1.0
        / (
            n_classes
            * len(donors_by_class[label])
            * counts[(donor, label)]
        )
        for donor, label in zip(donors, labels)
    ]


def score_rows(
    *,
    observed: Sequence[str],
    predicted: Sequence[str],
    donors: Sequence[str],
    probabilities: Sequence[Mapping[str, float]],
    class_roster: Sequence[str],
) -> dict[str, Any]:
    if not (
        observed
        and len(observed) == len(predicted) == len(donors) == len(probabilities)
    ):
        raise CellEvaluationError("evaluation rows are empty or misaligned")
    roster = tuple(class_roster)
    if tuple(sorted(set(roster))) != roster or set(observed) != set(roster):
        raise CellEvaluationError("observations do not exactly cover the frozen roster")
    if not set(predicted).issubset(roster):
        raise CellEvaluationError("prediction falls outside the frozen roster")
    weights = donor_class_balanced_weights(donors, observed)
    per_class: dict[str, float] = {}
    for class_id in roster:
        true_positive = sum(
            weight
            for truth, guess, weight in zip(observed, predicted, weights)
            if truth == class_id and guess == class_id
        )
        false_positive = sum(
            weight
            for truth, guess, weight in zip(observed, predicted, weights)
            if truth != class_id and guess == class_id
        )
        false_negative = sum(
            weight
            for truth, guess, weight in zip(observed, predicted, weights)
            if truth == class_id and guess != class_id
        )
        denominator = 2.0 * true_positive + false_positive + false_negative
        per_class[class_id] = 0.0 if denominator == 0 else 2.0 * true_positive / denominator
    macro_f1 = sum(per_class.values()) / len(roster)
    total_weight = sum(weights)
    balanced_accuracy = sum(
        weight
        for truth, guess, weight in zip(observed, predicted, weights)
        if truth == guess
    ) / total_weight
    brier = 0.0
    for truth, row, weight in zip(observed, probabilities, weights):
        if set(row) != set(roster):
            raise CellEvaluationError("probability row differs from the frozen roster")
        values = [float(row[class_id]) for class_id in roster]
        if (
            any(not math.isfinite(value) or value < 0 for value in values)
            or abs(sum(values) - 1.0) > 1e-6
        ):
            raise CellEvaluationError("probability row is invalid")
        brier += weight * sum(
            (row[class_id] - float(class_id == truth)) ** 2 for class_id in roster
        )
    return {
        "donor_class_balanced_macro_f1": macro_f1,
        "per_class_f1": per_class,
        "donor_balanced_accuracy": balanced_accuracy,
        "multiclass_brier_score": brier / total_weight,
    }


def _validate_environment(path: Path, expected_sha256: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or _sha256_file(path) != expected_sha256:
        raise CellEvaluationError("environment lock changed")
    try:
        lock = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CellEvaluationError("environment lock is invalid") from error
    if (
        not isinstance(lock, dict)
        or lock.get("schema_version") != "masld-bench-python-environment-lock-v1"
        or os.path.realpath(sys_executable()) != lock.get("python_executable")
        or os.path.realpath(os.sys.prefix) != lock.get("environment_prefix")
        or platform.python_version() != lock.get("python_version")
        or os.environ.get("PYTHONNOUSERSITE") != "1"
    ):
        raise CellEvaluationError("active evaluator environment differs from the lock")
    critical = lock.get("critical_versions")
    if not isinstance(critical, dict) or {
        name: importlib.metadata.version(name) for name in critical
    } != critical:
        raise CellEvaluationError("evaluator scientific packages differ from the lock")
    return lock


def sys_executable() -> str:
    return os.sys.executable


def _verified_bundle_artifact(
    root: Path, record: Mapping[str, Any], expected_role: str
) -> Path:
    relative = Path(str(record.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts or record.get("role") != expected_role:
        raise CellEvaluationError(f"bundle artifact role/path differs: {expected_role}")
    path = root / relative
    if path.is_symlink() or not path.is_file():
        raise CellEvaluationError(f"bundle artifact is missing: {relative}")
    if path.stat().st_size != record.get("size_bytes") or _sha256_file(path) != record.get(
        "sha256"
    ):
        raise CellEvaluationError(f"bundle artifact changed: {relative}")
    return path


def _load_prediction_bundle(
    path: Path, class_roster: Sequence[str], expected_model_id: str
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if path.is_symlink() or not path.is_file():
        raise CellEvaluationError("prediction bundle is missing")
    try:
        bundle = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CellEvaluationError("prediction bundle is invalid JSON") from error
    if not isinstance(bundle, dict):
        raise CellEvaluationError("prediction bundle must be an object")
    if (
        not expected_model_id
        or any(character in expected_model_id for character in "\t\r\n/\\")
        or bundle.get("model_id") != expected_model_id
    ):
        raise CellEvaluationError("prediction bundle model identity differs")
    exact = {
        "schema_version": "masld-bench-prediction-bundle-v1",
        "task_id": TASK_ID,
        "dataset_ids": [DATASET_ID],
        "split_id": "donor_outer",
        "row_id_field": "row_hash",
        "unit_id_field": "unit_hash",
        "biological_unit": "donor",
        "format_version": "tsv-v1",
        "missing_state": "observed",
    }
    for field, expected in exact.items():
        if bundle.get(field) != expected:
            raise CellEvaluationError(f"prediction bundle {field} differs")
    root = path.parent
    prediction_path = _verified_bundle_artifact(
        root,
        bundle["standardized_table"],
        f"standardized_prediction_table:{TASK_ID}",
    )
    row_path = _verified_bundle_artifact(
        root, bundle["row_ids"], f"prediction_row_ids:{TASK_ID}"
    )
    artifacts = bundle.get("artifacts")
    if not isinstance(artifacts, list):
        raise CellEvaluationError("prediction bundle artifact inventory is invalid")
    probability_records = [
        item
        for item in artifacts
        if isinstance(item, Mapping)
        and item.get("role") == f"class_probabilities:{TASK_ID}"
    ]
    if len(probability_records) != 1:
        raise CellEvaluationError("prediction bundle has no unique probability artifact")
    probability_path = _verified_bundle_artifact(
        root, probability_records[0], f"class_probabilities:{TASK_ID}"
    )
    prediction_fields, prediction_rows = _read_tsv(prediction_path)
    if prediction_fields != ("row_hash", "unit_hash", "predicted_class"):
        raise CellEvaluationError("prediction table schema differs or exposes labels")
    row_fields, row_rows = _read_tsv(row_path)
    if row_fields != ("row_hash", "unit_hash"):
        raise CellEvaluationError("prediction row-ID schema differs")
    probability_fields, probability_rows = _read_tsv(probability_path)
    expected_probability_fields = (
        "row_hash",
        "unit_hash",
        *(f"probability::{class_id}" for class_id in class_roster),
    )
    if probability_fields != expected_probability_fields:
        raise CellEvaluationError("probability table schema differs")
    if bundle.get("table_schema_sha256") != _canonical_hash(
        {"format": "tsv", "fields": list(prediction_fields)}
    ):
        raise CellEvaluationError("prediction table schema hash differs")
    prediction_index = {row["row_hash"]: row for row in prediction_rows}
    probability_index = {row["row_hash"]: row for row in probability_rows}
    row_index = {row["row_hash"]: row for row in row_rows}
    if (
        len(prediction_index) != len(prediction_rows)
        or len(probability_index) != len(probability_rows)
        or len(row_index) != len(row_rows)
        or not (set(prediction_index) == set(probability_index) == set(row_index))
    ):
        raise CellEvaluationError("prediction row inventories differ")
    if bundle.get("n_predictions") != len(prediction_rows):
        raise CellEvaluationError("prediction count differs")
    rows: list[dict[str, Any]] = []
    for row_hash in sorted(prediction_index):
        prediction = prediction_index[row_hash]
        probability = probability_index[row_hash]
        row_id = row_index[row_hash]
        unit_hash = prediction["unit_hash"]
        if probability["unit_hash"] != unit_hash or row_id["unit_hash"] != unit_hash:
            raise CellEvaluationError("prediction unit hashes differ across artifacts")
        values = {
            class_id: float(probability[f"probability::{class_id}"])
            for class_id in class_roster
        }
        winner = min(class_roster, key=lambda item: (-values[item], item))
        if prediction["predicted_class"] != winner:
            raise CellEvaluationError("hard class differs from probability argmax")
        rows.append(
            {
                "row_hash": row_hash,
                "unit_hash": unit_hash,
                "predicted_class": winner,
                "probabilities": values,
            }
        )
    return bundle, rows


def evaluate(
    *,
    bundle_path: Path,
    h5ad_path: Path,
    h5ad_sha256: str,
    environment_lock: Path,
    environment_sha256: str,
    class_roster: Sequence[str],
    expected_model_id: str,
    fold: int,
    outer_folds: int,
    split_seed: int,
    join_namespace: str,
    output_root: Path,
) -> Path:
    import anndata

    roster = tuple(class_roster)
    if tuple(sorted(set(roster))) != roster:
        raise CellEvaluationError("class roster must be sorted and unique")
    _validate_environment(environment_lock, _require_sha256(environment_sha256, "environment SHA"))
    if h5ad_path.is_symlink() or not h5ad_path.is_file() or _sha256_file(h5ad_path) != _require_sha256(
        h5ad_sha256, "H5AD SHA"
    ):
        raise CellEvaluationError("outcome-source H5AD changed")
    bundle, predictions = _load_prediction_bundle(
        bundle_path, roster, expected_model_id
    )
    expected_source_join = _canonical_hash(
        {
            "task_id": TASK_ID,
            "dataset_ids": [DATASET_ID],
            "split_id": "donor_outer",
            "row_id_field": "row_hash",
            "unit_id_field": "unit_hash",
            "unit_id_namespace": join_namespace,
            "biological_unit": "donor",
        }
    )
    if (
        bundle.get("unit_id_namespace") != join_namespace
        or bundle.get("source_join_key_sha256") != expected_source_join
    ):
        raise CellEvaluationError("prediction source-join binding differs")
    adata = anndata.read_h5ad(h5ad_path)
    if not {"donor_id", "broad_label"}.issubset(adata.obs.columns):
        raise CellEvaluationError("outcome-source H5AD lacks donor or class labels")
    outcomes: dict[str, dict[str, str]] = {}
    train_units: set[str] = set()
    held_units: set[str] = set()
    for row_id, donor, observed in zip(
        map(str, adata.obs_names),
        map(str, adata.obs["donor_id"]),
        map(str, adata.obs["broad_label"]),
    ):
        assigned = _fold_index(donor, seed=split_seed, outer_folds=outer_folds)
        unit_hash = _join_hash(join_namespace, "unit", donor)
        if assigned == fold:
            row_hash = _join_hash(join_namespace, "row", row_id)
            if row_hash in outcomes:
                raise CellEvaluationError("outcome row hash repeats")
            outcomes[row_hash] = {
                "unit_hash": unit_hash,
                "observed_class": observed,
            }
            held_units.add(unit_hash)
        else:
            train_units.add(unit_hash)
    if train_units & held_units:
        raise CellEvaluationError("one donor crosses train and held folds")
    prediction_index = {row["row_hash"]: row for row in predictions}
    if set(prediction_index) != set(outcomes):
        raise CellEvaluationError("predictions do not exactly cover the held fold")
    observed: list[str] = []
    predicted: list[str] = []
    donors: list[str] = []
    probabilities: list[Mapping[str, float]] = []
    evaluated_rows: list[dict[str, Any]] = []
    for row_hash in sorted(outcomes):
        outcome = outcomes[row_hash]
        prediction = prediction_index[row_hash]
        if prediction["unit_hash"] != outcome["unit_hash"]:
            raise CellEvaluationError("prediction donor hash differs from the outcome join")
        observed.append(outcome["observed_class"])
        predicted.append(prediction["predicted_class"])
        donors.append(outcome["unit_hash"])
        probabilities.append(prediction["probabilities"])
        evaluated_rows.append(
            {
                "row_hash": row_hash,
                "unit_hash": outcome["unit_hash"],
                "observed_class": outcome["observed_class"],
                "predicted_class": prediction["predicted_class"],
            }
        )
    metrics = score_rows(
        observed=observed,
        predicted=predicted,
        donors=donors,
        probabilities=probabilities,
        class_roster=roster,
    )
    identity = {
        "bundle_sha256": _sha256_file(bundle_path),
        "model_id": expected_model_id,
        "h5ad_sha256": h5ad_sha256,
        "environment_sha256": environment_sha256,
        "evaluator_source_sha256": _sha256_file(Path(__file__).resolve(strict=True)),
        "class_roster": list(roster),
        "fold": fold,
        "outer_folds": outer_folds,
        "split_seed": split_seed,
        "join_namespace": join_namespace,
    }
    result_id = _canonical_hash(identity)
    output_root.mkdir(parents=True, exist_ok=True)
    target = output_root / f"cell-state-smoke--{result_id[:16]}"
    target.mkdir(exist_ok=False)
    rows_path = target / "evaluated_rows.tsv"
    _write_tsv(
        rows_path,
        ("row_hash", "unit_hash", "observed_class", "predicted_class"),
        evaluated_rows,
    )
    result_path = target / "evaluation.json"
    _write_json(
        result_path,
        {
            "schema_version": "masld-bench-development-cell-evaluation-v1",
            "result_id": result_id,
            "task_id": TASK_ID,
            "model_id": expected_model_id,
            "evaluator_id": EVALUATOR_ID,
            "metrics": metrics,
            "independent_unit": "donor",
            "independent_unit_count": len(set(donors)),
            "row_count": len(observed),
            "smoke_only": True,
            "confirmatory": False,
            "champion_eligible": False,
            "identity": identity,
        },
    )
    artifacts = [
        {
            "path": path.name,
            "sha256": _sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in (result_path, rows_path)
    ]
    manifest_path = target / "ARTIFACTS.json"
    _write_json(
        manifest_path,
        {
            "schema_version": "masld-bench-artifacts-v1",
            "metadata": {
                "artifact_class": "development_cell_evaluation",
                "result_id": result_id,
            },
            "artifacts": artifacts,
        },
    )
    _write_json(
        target / "COMPLETE",
        {
            "schema_version": "masld-bench-complete-v1",
            "manifest_sha256": _sha256_file(manifest_path),
            "artifact_count": len(artifacts),
        },
    )
    for path in target.iterdir():
        path.chmod(0o440)
    target.chmod(0o550)
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--h5ad-sha256", required=True)
    parser.add_argument("--environment-lock", required=True, type=Path)
    parser.add_argument("--environment-sha256", required=True)
    parser.add_argument("--class", dest="classes", action="append", required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--fold", required=True, type=int)
    parser.add_argument("--outer-folds", required=True, type=int)
    parser.add_argument("--split-seed", required=True, type=int)
    parser.add_argument("--join-namespace", required=True)
    parser.add_argument("--output-root", required=True, type=Path)
    arguments = parser.parse_args(argv)
    target = evaluate(
        bundle_path=arguments.bundle.resolve(strict=True),
        h5ad_path=arguments.h5ad.resolve(strict=True),
        h5ad_sha256=arguments.h5ad_sha256,
        environment_lock=arguments.environment_lock.resolve(strict=True),
        environment_sha256=arguments.environment_sha256,
        class_roster=arguments.classes,
        expected_model_id=arguments.model_id,
        fold=arguments.fold,
        outer_folds=arguments.outer_folds,
        split_seed=arguments.split_seed,
        join_namespace=arguments.join_namespace,
        output_root=arguments.output_root,
    )
    print(_canonical_json({"evaluation": target.as_posix()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
