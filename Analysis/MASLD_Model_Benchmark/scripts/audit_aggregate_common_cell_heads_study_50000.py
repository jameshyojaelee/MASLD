#!/usr/bin/env python3
"""Audit and aggregate frozen 50,000-cell common-head prediction shards."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
HEADS = ("linear", "two_layer_mlp")
SCREEN_SEEDS = (1103, 1201, 1301)
OUTER_FOLDS = (0, 1, 2, 3, 4)
EXPECTED_ROWS = 50_000
EXPECTED_DONORS = 102
EXPECTED_STUDIES = 7
FIXED_EPOCHS = 30
MODEL_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9._-]{2,79}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
EXPOSURE_STATUSES = frozenset(
    {
        "clean_declared",
        "target_label_unexposed",
        "encoder_seen",
        "continual_seen",
        "reference_only",
        "downstream_demo",
        "unknown",
    }
)
PREDICTION_FIELDS = (
    "row_id",
    "donor_id",
    "dataset",
    "outer_fold",
    "predicted_class",
    *(f"probability::{label}" for label in ROSTER),
)


class CommonHeadAuditError(ValueError):
    """Raised when a frozen shard bundle does not meet the common-lane requirements."""


def require_model_id(value: str) -> str:
    if MODEL_ID_PATTERN.fullmatch(value) is None:
        raise CommonHeadAuditError("expected model_id is invalid")
    return value


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CommonHeadAuditError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


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


def validate_donor_safe_split(rows: Sequence[Mapping[str, str]]) -> None:
    if len(rows) != EXPECTED_ROWS:
        raise CommonHeadAuditError("split row count differs")
    row_ids = [row["row_id"] for row in rows]
    donor_folds: dict[str, set[int]] = {}
    study_folds: dict[str, set[int]] = {}
    for row in rows:
        try:
            fold = int(row["outer_fold"])
        except (KeyError, ValueError) as error:
            raise CommonHeadAuditError("split outer fold is invalid") from error
        donor_folds.setdefault(row["donor_id"], set()).add(fold)
        study_folds.setdefault(row["dataset"], set()).add(fold)
    if (
        len(set(row_ids)) != EXPECTED_ROWS
        or len(donor_folds) != EXPECTED_DONORS
        or len(study_folds) != EXPECTED_STUDIES
        or set().union(*donor_folds.values()) != set(OUTER_FOLDS)
        or any(len(value) != 1 for value in donor_folds.values())
        or any(len(value) != 1 for value in study_folds.values())
    ):
        raise CommonHeadAuditError("split is not donor- and study-safe")


def validate_probabilities(rows: Sequence[Mapping[str, str]]) -> None:
    for row in rows:
        try:
            values = np.asarray(
                [float(row[f"probability::{label}"]) for label in ROSTER],
                dtype=np.float64,
            )
        except (KeyError, TypeError, ValueError) as error:
            raise CommonHeadAuditError("prediction probabilities are invalid") from error
        if (
            not np.all(np.isfinite(values))
            or np.any(values < 0.0)
            or not np.isclose(values.sum(), 1.0, rtol=0.0, atol=1.0e-12)
            or row.get("predicted_class") != ROSTER[int(np.argmax(values))]
        ):
            raise CommonHeadAuditError("prediction probability contract differs")


def expected_shard_keys() -> set[tuple[str, int, int]]:
    return {
        (head, seed, fold)
        for head in HEADS
        for seed in SCREEN_SEEDS
        for fold in OUTER_FOLDS
    }


def run(
    *,
    bundle: Path,
    split: Path,
    output: Path,
    expected_bundle_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    expected_source_artifacts_sha256: str,
    expected_model_id: str,
) -> None:
    if output.exists():
        raise CommonHeadAuditError("refusing to overwrite audited predictions")
    expected_model_id = require_model_id(expected_model_id)
    if any(
        SHA256_PATTERN.fullmatch(value) is None
        for value in (
            expected_bundle_artifacts_sha256,
            expected_split_artifacts_sha256,
            expected_source_artifacts_sha256,
        )
    ):
        raise CommonHeadAuditError("expected ARTIFACTS SHA-256 is invalid")
    for root, expected in (
        (bundle, expected_bundle_artifacts_sha256),
        (split, expected_split_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise CommonHeadAuditError("input ARTIFACTS SHA-256 differs")
        verify_frozen_tree(root)

    bundle_metadata = verify_frozen_tree(bundle)["metadata"]
    split_metadata = verify_frozen_tree(split)["metadata"]
    if bundle_metadata != {
        "artifact_class": "cell_foundation_common_head_shard_bundle",
        "head_ids": list(HEADS),
        "logical_shards": len(expected_shard_keys()),
        "metrics_calculated": False,
        "outer_folds": list(OUTER_FOLDS),
        "screen_seeds": list(SCREEN_SEEDS),
        "sealed_outcomes_read": False,
        "status": "passed",
    }:
        raise CommonHeadAuditError("bundle metadata differs")
    if (
        split_metadata.get("artifact_class") != "cell_state_study_outer_split"
        or split_metadata.get("dataset_view_id") != DATASET_VIEW_ID
        or split_metadata.get("split_id") != SPLIT_ID
        or split_metadata.get("rows") != EXPECTED_ROWS
        or split_metadata.get("donors") != EXPECTED_DONORS
        or split_metadata.get("studies") != EXPECTED_STUDIES
        or split_metadata.get("target_labels_used_for_assignment") is not False
        or split_metadata.get("sealed_outcomes_read") is not False
        or split_metadata.get("histology_read") is not False
    ):
        raise CommonHeadAuditError("split metadata differs")

    split_fields, split_rows = read_tsv(split / "row_outer_folds.tsv")
    if split_fields != ["row_id", "donor_id", "dataset", "outer_fold"]:
        raise CommonHeadAuditError("split table schema differs")
    validate_donor_safe_split(split_rows)
    row_ids = [row["row_id"] for row in split_rows]
    row_position = {row_id: index for index, row_id in enumerate(row_ids)}
    expected_by_row = {
        row["row_id"]: (row["donor_id"], row["dataset"], int(row["outer_fold"]))
        for row in split_rows
    }
    fold_rows = {
        fold: sum(int(row["outer_fold"]) == fold for row in split_rows)
        for fold in OUTER_FOLDS
    }
    fold_donors = {
        fold: len(
            {row["donor_id"] for row in split_rows if int(row["outer_fold"]) == fold}
        )
        for fold in OUTER_FOLDS
    }
    fold_studies = {
        fold: sorted(
            {row["dataset"] for row in split_rows if int(row["outer_fold"]) == fold}
        )
        for fold in OUTER_FOLDS
    }

    observed: dict[tuple[str, int, int], list[dict[str, str]]] = {}
    audit_rows: list[dict[str, Any]] = []
    identity: dict[str, Any] | None = None
    shard_roots = sorted((bundle / "shards").glob("*"))
    if any(not path.is_dir() for path in shard_roots):
        raise CommonHeadAuditError("shard root contains a non-directory entry")
    for shard_root in shard_roots:
        verify_frozen_tree(shard_root)
        receipt = json.loads(
            (shard_root / "prediction_receipt.json").read_text(encoding="utf-8")
        )
        try:
            key = (
                str(receipt["head_id"]),
                int(receipt["screen_seed"]),
                int(receipt["outer_fold"]),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise CommonHeadAuditError("shard identity is invalid") from error
        if key not in expected_shard_keys() or key in observed:
            raise CommonHeadAuditError("shard census contains an unexpected duplicate")
        if shard_root.name != f"{key[0]}__seed{key[1]}__fold{key[2]}":
            raise CommonHeadAuditError("shard directory and receipt identities differ")
        if (
            receipt.get("status") != "pass_development_prediction_shard"
            or receipt.get("schema_version")
            != "masld-bench-common-cell-head-shard-v1"
            or receipt.get("model_id") != expected_model_id
            or receipt.get("representation_id") != "common"
            or receipt.get("dataset_view_id") != DATASET_VIEW_ID
            or receipt.get("split_id") != SPLIT_ID
            or receipt.get("fixed_epochs") != FIXED_EPOCHS
            or receipt.get("checkpoint_selection_uses_validation_labels") is not False
            or receipt.get("inner_models_used_for_outer_prediction") is not False
            or receipt.get("outer_prediction_uses_final_refit") is not True
            or receipt.get("outer_prediction_model")
            != "fixed_epoch_final_outer_training_refit"
            or receipt.get("final_fit_rows") != receipt.get("training_rows")
            or receipt.get("final_fit_donors") != receipt.get("training_donors")
            or receipt.get("final_fit_studies")
            != len(receipt.get("training_studies", []))
            or receipt.get("prediction_tables_contain_observed_labels") is not False
            or receipt.get("metrics_calculated") is not False
            or receipt.get("histology_read") is not False
            or receipt.get("sealed_outcomes_read") is not False
            or receipt.get("development_labels_read") != ["broad_label"]
            or receipt.get("input_artifacts_sha256", {}).get("source")
            != expected_source_artifacts_sha256
            or receipt.get("input_artifacts_sha256", {}).get("split")
            != expected_split_artifacts_sha256
        ):
            raise CommonHeadAuditError("prediction shard receipt differs")
        fold = key[2]
        expected_test_studies = fold_studies[fold]
        expected_training_studies = sorted(
            {study for values in fold_studies.values() for study in values}
            - set(expected_test_studies)
        )
        inner_models = receipt.get("inner_models", [])
        if (
            receipt.get("test_rows") != fold_rows[fold]
            or receipt.get("test_donors") != fold_donors[fold]
            or receipt.get("training_rows") != EXPECTED_ROWS - fold_rows[fold]
            or receipt.get("training_donors") != EXPECTED_DONORS - fold_donors[fold]
            or sorted(receipt.get("test_studies", [])) != expected_test_studies
            or sorted(receipt.get("training_studies", [])) != expected_training_studies
            or set(receipt.get("test_studies", []))
            & set(receipt.get("training_studies", []))
            or receipt.get("oof_exact_one_time_coverage") is not True
            or receipt.get("oof_rows") != receipt.get("training_rows")
            or len(inner_models) != 5
            or {item.get("inner_fold") for item in inner_models} != set(range(5))
            or any(
                item.get("fixed_epochs") != FIXED_EPOCHS
                or item.get("validation_labels_used_for_checkpoint_selection") is not False
                or item.get("fit_rows", 0) + item.get("prediction_rows", 0)
                != receipt.get("training_rows")
                or item.get("fit_donors", 0) + item.get("prediction_donors", 0)
                != receipt.get("training_donors")
                for item in inner_models
            )
        ):
            raise CommonHeadAuditError("shard donor-safe outer-fold receipt differs")
        observed_identity = {
            field: receipt[field]
            for field in (
                "model_id",
                "representation_id",
                "checkpoint_sha256",
                "exposure_status",
                "embedding_width",
                "embedding_relative",
                "sealed_champion_eligible",
            )
        }
        observed_identity["embeddings_artifacts_sha256"] = receipt[
            "input_artifacts_sha256"
        ]["embeddings"]
        if (
            observed_identity["representation_id"] != "common"
            or SHA256_PATTERN.fullmatch(
                str(observed_identity["checkpoint_sha256"])
            )
            is None
            or observed_identity["exposure_status"] not in EXPOSURE_STATUSES
            or type(observed_identity["embedding_width"]) is not int
            or observed_identity["embedding_width"] < 1
            or observed_identity["embedding_relative"]
            != "embeddings/common_embeddings.npz"
            or type(observed_identity["sealed_champion_eligible"]) is not bool
            or SHA256_PATTERN.fullmatch(
                str(observed_identity["embeddings_artifacts_sha256"])
            )
            is None
            or (
                observed_identity["sealed_champion_eligible"]
                and observed_identity["exposure_status"]
                not in {"clean_declared", "target_label_unexposed"}
            )
        ):
            raise CommonHeadAuditError("shard model identity contract differs")
        if identity is None:
            identity = observed_identity
        elif observed_identity != identity:
            raise CommonHeadAuditError("shard model identities differ")

        fields, rows = read_tsv(shard_root / "predictions.tsv")
        if fields != list(PREDICTION_FIELDS):
            raise CommonHeadAuditError("prediction shard schema differs")
        if any("label" in field.lower() or "histolog" in field.lower() for field in fields):
            raise CommonHeadAuditError("prediction shard exposes evaluator-only fields")
        if len(rows) != fold_rows[fold] or len({row["row_id"] for row in rows}) != len(rows):
            raise CommonHeadAuditError("prediction shard row coverage differs")
        for row in rows:
            expected = expected_by_row.get(row["row_id"])
            if expected is None or (
                row["donor_id"], row["dataset"], int(row["outer_fold"])
            ) != expected or expected[2] != fold:
                raise CommonHeadAuditError("prediction row and frozen split differ")
        validate_probabilities(rows)
        observed[key] = rows
        audit_rows.append(
            {
                "head_id": key[0],
                "screen_seed": key[1],
                "outer_fold": fold,
                "rows": len(rows),
                "donors": len({row["donor_id"] for row in rows}),
                "studies": ",".join(sorted({row["dataset"] for row in rows})),
                "artifacts_sha256": sha256_file(shard_root / "ARTIFACTS.json"),
                "predictions_sha256": sha256_file(shard_root / "predictions.tsv"),
            }
        )
    if set(observed) != expected_shard_keys() or identity is None:
        raise CommonHeadAuditError("common-head shard census is incomplete")

    output.mkdir(mode=0o750)
    prediction_root = output / "predictions"
    prediction_root.mkdir()
    for head in HEADS:
        for seed in SCREEN_SEEDS:
            combined = [row for fold in OUTER_FOLDS for row in observed[(head, seed, fold)]]
            if len(combined) != EXPECTED_ROWS or set(row["row_id"] for row in combined) != set(row_ids):
                raise CommonHeadAuditError("aggregated row coverage differs")
            combined.sort(key=lambda row: row_position[row["row_id"]])
            if [row["row_id"] for row in combined] != row_ids:
                raise CommonHeadAuditError("aggregated row order differs")
            write_tsv(
                prediction_root / f"{head}__seed{seed}.tsv",
                PREDICTION_FIELDS,
                combined,
            )
    write_tsv(
        output / "shard_audit.tsv",
        [
            "head_id",
            "screen_seed",
            "outer_fold",
            "rows",
            "donors",
            "studies",
            "artifacts_sha256",
            "predictions_sha256",
        ],
        sorted(audit_rows, key=lambda row: (row["head_id"], row["screen_seed"], row["outer_fold"])),
    )
    receipt = {
        "schema_version": "masld-bench-common-cell-head-audited-predictions-v1",
        "status": "pass_development_predictions_audited",
        **identity,
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "studies": EXPECTED_STUDIES,
        "heads": list(HEADS),
        "screen_seeds": list(SCREEN_SEEDS),
        "outer_folds": list(OUTER_FOLDS),
        "logical_shards": len(expected_shard_keys()),
        "prediction_tables": len(HEADS) * len(SCREEN_SEEDS),
        "fixed_epochs": FIXED_EPOCHS,
        "bundle_artifacts_sha256": expected_bundle_artifacts_sha256,
        "split_artifacts_sha256": expected_split_artifacts_sha256,
        "source_artifacts_sha256_from_shards": expected_source_artifacts_sha256,
        "row_ids_aligned_to_split": True,
        "donor_safe_outer_folds_verified": True,
        "whole_study_outer_folds_verified": True,
        "head_seed_fold_census_complete": True,
        "probability_schema_verified": True,
        "checkpoint_selection_uses_validation_labels": False,
        "metrics_calculated": False,
        "development_labels_read_by_aggregator": [],
        "prediction_tables_contain_observed_labels": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "clinical_claim_allowed": False,
    }
    write_json_exclusive(output / "prediction_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_foundation_common_head_audited_predictions",
            "model_id": expected_model_id,
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": EXPECTED_ROWS,
            "logical_shards": len(expected_shard_keys()),
            "metrics_calculated": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--bundle", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-bundle-artifacts-sha256", required=True)
    value.add_argument("--expected-split-artifacts-sha256", required=True)
    value.add_argument("--expected-source-artifacts-sha256", required=True)
    value.add_argument("--expected-model-id", required=True)
    return value


def main() -> int:
    arguments = parser().parse_args()
    run(
        bundle=arguments.bundle,
        split=arguments.split,
        output=arguments.output,
        expected_bundle_artifacts_sha256=arguments.expected_bundle_artifacts_sha256,
        expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
        expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
        expected_model_id=arguments.expected_model_id,
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
