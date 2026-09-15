#!/usr/bin/env python3
"""Convert frozen raw cell embeddings into the donor-safe common-head requirements.

This preparation step is outcome-blind. It reads raw embeddings, their row
order, the frozen study split, and model authorities. It never opens the source
label table or any held-back outcome.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
OUTER_FOLDS = frozenset(range(5))
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


class CommonEmbeddingPreparationError(ValueError):
    """Raised when raw embeddings cannot enter the common-head requirements."""


def require_sha256(value: str, label: str) -> str:
    if SHA256_PATTERN.fullmatch(value) is None:
        raise CommonEmbeddingPreparationError(f"{label} is not a lowercase SHA-256")
    return value


def resolve_inside(root: Path, relative: Path, label: str) -> Path:
    if relative.is_absolute() or relative == Path("."):
        raise CommonEmbeddingPreparationError(f"{label} must be a relative file")
    try:
        path = (root / relative).resolve(strict=True)
        path.relative_to(root.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise CommonEmbeddingPreparationError(f"{label} escapes its frozen root") from error
    if not path.is_file():
        raise CommonEmbeddingPreparationError(f"{label} is not a regular file")
    return path


def read_split_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["row_id", "donor_id", "dataset", "outer_fold"]:
            raise CommonEmbeddingPreparationError("split table schema differs")
        return list(reader)


def validate_split(
    rows: list[dict[str, str]], *, expected_rows: int, expected_donors: int, expected_studies: int
) -> np.ndarray:
    try:
        folds = np.asarray([int(row["outer_fold"]) for row in rows], dtype=np.int8)
    except (KeyError, TypeError, ValueError) as error:
        raise CommonEmbeddingPreparationError("split fold values are invalid") from error
    donors = np.asarray([row["donor_id"] for row in rows], dtype=str)
    studies = np.asarray([row["dataset"] for row in rows], dtype=str)
    row_ids = [row["row_id"] for row in rows]
    if (
        len(rows) != expected_rows
        or len(set(row_ids)) != expected_rows
        or len(set(donors.tolist())) != expected_donors
        or len(set(studies.tolist())) != expected_studies
        or set(map(int, folds)) != OUTER_FOLDS
        or any(len(set(map(int, folds[donors == donor]))) != 1 for donor in set(donors))
        or any(len(set(map(int, folds[studies == study]))) != 1 for study in set(studies))
    ):
        raise CommonEmbeddingPreparationError("split is not the frozen donor/study-safe partition")
    return folds


def checkpoint_record(
    registry: dict[str, Any], *, model_id: str, artifact_path: str, checkpoint_sha256: str
) -> dict[str, Any]:
    if model_id not in registry.get("architecture_contract", {}):
        raise CommonEmbeddingPreparationError("model_id is absent from checkpoint authority")
    matches = [
        item
        for item in registry.get("artifacts", [])
        if item.get("path") == artifact_path
        and item.get("sha256") == checkpoint_sha256
    ]
    if len(matches) != 1:
        raise CommonEmbeddingPreparationError("checkpoint identity differs from authority")
    return matches[0]


def run(
    *,
    raw_root: Path,
    raw_embeddings_relative: Path,
    raw_row_order_relative: Path,
    raw_receipt_relative: Path,
    split_root: Path,
    checkpoint_registry: Path,
    exposure_audit: Path,
    output: Path,
    expected_raw_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    expected_checkpoint_registry_sha256: str,
    expected_exposure_audit_sha256: str,
    expected_model_id: str,
    expected_checkpoint_artifact_path: str,
    expected_checkpoint_sha256: str,
    expected_exposure_status: str,
    expected_sealed_champion_eligible: bool,
    expected_rows: int,
    expected_width: int,
    expected_donors: int,
    expected_studies: int,
    dataset_view_id: str = DATASET_VIEW_ID,
    split_id: str = SPLIT_ID,
) -> None:
    if output.exists():
        raise CommonEmbeddingPreparationError("refusing to overwrite common embeddings")
    if MODEL_ID_PATTERN.fullmatch(expected_model_id) is None:
        raise CommonEmbeddingPreparationError("expected model_id is invalid")
    if expected_exposure_status not in EXPOSURE_STATUSES:
        raise CommonEmbeddingPreparationError("expected exposure status is invalid")
    if type(expected_sealed_champion_eligible) is not bool:
        raise CommonEmbeddingPreparationError("champion eligibility is not boolean")
    if expected_sealed_champion_eligible and expected_exposure_status not in {
        "clean_declared",
        "target_label_unexposed",
    }:
        raise CommonEmbeddingPreparationError("champion eligibility contradicts exposure")
    if min(expected_rows, expected_width, expected_donors, expected_studies) < 1:
        raise CommonEmbeddingPreparationError("expected dimensions must be positive")
    expected_hashes = {
        raw_root / "ARTIFACTS.json": require_sha256(
            expected_raw_artifacts_sha256, "raw ARTIFACTS hash"
        ),
        split_root / "ARTIFACTS.json": require_sha256(
            expected_split_artifacts_sha256, "split ARTIFACTS hash"
        ),
        checkpoint_registry: require_sha256(
            expected_checkpoint_registry_sha256, "checkpoint authority hash"
        ),
        exposure_audit: require_sha256(
            expected_exposure_audit_sha256, "exposure authority hash"
        ),
    }
    require_sha256(expected_checkpoint_sha256, "checkpoint hash")
    for path, expected in expected_hashes.items():
        if not path.is_file() or sha256_file(path) != expected:
            raise CommonEmbeddingPreparationError(f"pinned input hash differs: {path}")

    raw_manifest = verify_frozen_tree(raw_root)
    split_manifest = verify_frozen_tree(split_root)
    raw_metadata = raw_manifest.get("metadata", {})
    split_metadata = split_manifest.get("metadata", {})
    if (
        raw_metadata.get("artifact_class")
        != "geneformer_v2_316m_frozen_screen_raw_embeddings"
        or raw_metadata.get("model_id") != expected_model_id
        or raw_metadata.get("dataset_view_id") != dataset_view_id
        or raw_metadata.get("rows") != expected_rows
        or raw_metadata.get("raw_unrefined_embeddings_only") is not True
        or raw_metadata.get("downstream_head_fit") is not False
        or raw_metadata.get("evaluation_labels_read") is not False
        or raw_metadata.get("histology_read") is not False
        or raw_metadata.get("sealed_outcomes_read") is not False
        or raw_metadata.get("status") != "passed"
    ):
        raise CommonEmbeddingPreparationError("raw embedding firewall metadata differs")
    if (
        split_metadata.get("artifact_class") != "cell_state_study_outer_split"
        or split_metadata.get("dataset_view_id") != dataset_view_id
        or split_metadata.get("split_id") != split_id
        or split_metadata.get("rows") != expected_rows
        or split_metadata.get("donors") != expected_donors
        or split_metadata.get("studies") != expected_studies
        or split_metadata.get("target_labels_used_for_assignment") is not False
        or split_metadata.get("sealed_outcomes_read") is not False
        or split_metadata.get("histology_read") is not False
        or split_metadata.get("status") != "passed"
    ):
        raise CommonEmbeddingPreparationError("split firewall metadata differs")

    raw_embeddings = resolve_inside(raw_root, raw_embeddings_relative, "raw embeddings")
    raw_row_order = resolve_inside(raw_root, raw_row_order_relative, "raw row order")
    raw_receipt_path = resolve_inside(raw_root, raw_receipt_relative, "raw receipt")
    raw_receipt = json.loads(raw_receipt_path.read_text(encoding="utf-8"))
    if (
        raw_receipt.get("model_id") != expected_model_id
        or raw_receipt.get("dataset_view_id") != dataset_view_id
        or raw_receipt.get("rows") != expected_rows
        or raw_receipt.get("embedding_shape") != [expected_rows, expected_width]
        or raw_receipt.get("embedding_dtype") != "float32"
        or raw_receipt.get("embedding_sha256") != sha256_file(raw_embeddings)
        or raw_receipt.get("row_order_sha256") != sha256_file(raw_row_order)
        or raw_receipt.get("fixture_row_order_sha256") != sha256_file(raw_row_order)
        or raw_receipt.get("raw_unrefined_embeddings_only") is not True
        or raw_receipt.get("downstream_head_fit") is not False
        or raw_receipt.get("evaluation_labels_read") is not False
        or raw_receipt.get("histology_read") is not False
        or raw_receipt.get("sealed_outcomes_read") is not False
        or raw_receipt.get("status") != "pass_raw_outcome_blind_extraction"
    ):
        raise CommonEmbeddingPreparationError("raw extraction receipt differs")

    registry = json.loads(checkpoint_registry.read_text(encoding="utf-8"))
    checkpoint = checkpoint_record(
        registry,
        model_id=expected_model_id,
        artifact_path=expected_checkpoint_artifact_path,
        checkpoint_sha256=expected_checkpoint_sha256,
    )
    exposure = json.loads(exposure_audit.read_text(encoding="utf-8"))
    finding = exposure.get("checkpoint_findings", {}).get(expected_model_id, {})
    if finding.get("exposure_state") != expected_exposure_status:
        raise CommonEmbeddingPreparationError("model exposure differs from authority")
    eligibility = exposure.get("sealed_champion_eligibility")
    if expected_sealed_champion_eligible and eligibility != (
        "eligible_on_exposure_only_other_admission_gates_still_apply"
    ):
        raise CommonEmbeddingPreparationError("sealed eligibility differs from authority")

    split_rows = read_split_rows(split_root / "row_outer_folds.tsv")
    outer_folds = validate_split(
        split_rows,
        expected_rows=expected_rows,
        expected_donors=expected_donors,
        expected_studies=expected_studies,
    )
    row_ids = raw_row_order.read_text(encoding="utf-8").splitlines()
    split_row_ids = [row["row_id"] for row in split_rows]
    if (
        len(row_ids) != expected_rows
        or len(set(row_ids)) != expected_rows
        or any(not row_id for row_id in row_ids)
        or row_ids != split_row_ids
    ):
        raise CommonEmbeddingPreparationError("raw embeddings and split row order differ")

    embeddings = np.load(raw_embeddings, mmap_mode="r", allow_pickle=False)
    if embeddings.shape != (expected_rows, expected_width) or embeddings.dtype != np.float32:
        raise CommonEmbeddingPreparationError("raw embedding array shape or dtype differs")
    for start in range(0, expected_rows, 2048):
        if not np.all(np.isfinite(embeddings[start : start + 2048])):
            raise CommonEmbeddingPreparationError("raw embeddings contain non-finite values")

    output.mkdir(mode=0o750)
    embedding_root = output / "embeddings"
    embedding_root.mkdir()
    representation = embedding_root / "common_embeddings.npz"
    row_dtype = f"<U{max(len(value) for value in row_ids)}"
    with representation.open("xb") as handle:
        np.savez(
            handle,
            embeddings=embeddings,
            outer_folds=outer_folds,
            row_ids=np.asarray(row_ids, dtype=row_dtype),
        )
    with np.load(representation, allow_pickle=False) as prepared:
        if (
            set(prepared.files) != {"embeddings", "outer_folds", "row_ids"}
            or prepared["embeddings"].shape != (expected_rows, expected_width)
            or prepared["embeddings"].dtype != np.float32
            or prepared["outer_folds"].tolist() != outer_folds.tolist()
            or prepared["row_ids"].astype(str).tolist() != row_ids
        ):
            raise CommonEmbeddingPreparationError("written common embedding bundle differs")

    receipt = {
        "schema_version": "masld-bench-common-cell-embeddings-v1",
        "status": "pass_outcome_blind_common_embedding_preparation",
        "model_id": expected_model_id,
        "checkpoint_artifact_path": expected_checkpoint_artifact_path,
        "checkpoint_sha256": expected_checkpoint_sha256,
        "checkpoint_training_cutoff": checkpoint.get("training_cutoff"),
        "exposure_status": expected_exposure_status,
        "sealed_champion_eligible": expected_sealed_champion_eligible,
        "dataset_view_id": dataset_view_id,
        "split_id": split_id,
        "rows": expected_rows,
        "embedding_width": expected_width,
        "embedding_dtype": "float32",
        "policies": {
            "common": {
                "transformation": "identity_from_frozen_raw_embedding",
                "source_embedding_policy": raw_receipt.get("embedding_policy"),
            }
        },
        "raw_artifacts_sha256": expected_raw_artifacts_sha256,
        "raw_embeddings_sha256": sha256_file(raw_embeddings),
        "raw_row_order_sha256": sha256_file(raw_row_order),
        "split_artifacts_sha256": expected_split_artifacts_sha256,
        "checkpoint_registry_sha256": expected_checkpoint_registry_sha256,
        "exposure_audit_sha256": expected_exposure_audit_sha256,
        "row_order_aligned_to_split": True,
        "donor_safe_outer_folds_verified": True,
        "whole_study_outer_folds_verified": True,
        "evaluation_label_columns_read": [],
        "source_label_table_read": False,
        "downstream_head_fit": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(embedding_root / "receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_foundation_common_embeddings",
            "model_id": expected_model_id,
            "dataset_view_id": dataset_view_id,
            "development_rows": expected_rows,
            "exposure_status": expected_exposure_status,
            "sealed_champion_eligible": expected_sealed_champion_eligible,
            "evaluation_labels_read": False,
            "source_label_table_read": False,
            "downstream_head_fit": False,
            "histology_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--raw-root", required=True, type=Path)
    value.add_argument("--raw-embeddings-relative", required=True, type=Path)
    value.add_argument("--raw-row-order-relative", required=True, type=Path)
    value.add_argument("--raw-receipt-relative", required=True, type=Path)
    value.add_argument("--split-root", required=True, type=Path)
    value.add_argument("--checkpoint-registry", required=True, type=Path)
    value.add_argument("--exposure-audit", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-raw-artifacts-sha256", required=True)
    value.add_argument("--expected-split-artifacts-sha256", required=True)
    value.add_argument("--expected-checkpoint-registry-sha256", required=True)
    value.add_argument("--expected-exposure-audit-sha256", required=True)
    value.add_argument("--expected-model-id", required=True)
    value.add_argument("--expected-checkpoint-artifact-path", required=True)
    value.add_argument("--expected-checkpoint-sha256", required=True)
    value.add_argument(
        "--expected-exposure-status", required=True, choices=sorted(EXPOSURE_STATUSES)
    )
    value.add_argument(
        "--expected-sealed-champion-eligible",
        required=True,
        choices=("true", "false"),
    )
    value.add_argument("--expected-rows", required=True, type=int)
    value.add_argument("--expected-width", required=True, type=int)
    value.add_argument("--expected-donors", required=True, type=int)
    value.add_argument("--expected-studies", required=True, type=int)
    value.add_argument("--dataset-view-id", default=DATASET_VIEW_ID)
    value.add_argument("--split-id", default=SPLIT_ID)
    return value


def main() -> int:
    arguments = parser().parse_args()
    run(
        raw_root=arguments.raw_root,
        raw_embeddings_relative=arguments.raw_embeddings_relative,
        raw_row_order_relative=arguments.raw_row_order_relative,
        raw_receipt_relative=arguments.raw_receipt_relative,
        split_root=arguments.split_root,
        checkpoint_registry=arguments.checkpoint_registry,
        exposure_audit=arguments.exposure_audit,
        output=arguments.output,
        expected_raw_artifacts_sha256=arguments.expected_raw_artifacts_sha256,
        expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
        expected_checkpoint_registry_sha256=arguments.expected_checkpoint_registry_sha256,
        expected_exposure_audit_sha256=arguments.expected_exposure_audit_sha256,
        expected_model_id=arguments.expected_model_id,
        expected_checkpoint_artifact_path=arguments.expected_checkpoint_artifact_path,
        expected_checkpoint_sha256=arguments.expected_checkpoint_sha256,
        expected_exposure_status=arguments.expected_exposure_status,
        expected_sealed_champion_eligible=(
            arguments.expected_sealed_champion_eligible == "true"
        ),
        expected_rows=arguments.expected_rows,
        expected_width=arguments.expected_width,
        expected_donors=arguments.expected_donors,
        expected_studies=arguments.expected_studies,
        dataset_view_id=arguments.dataset_view_id,
        split_id=arguments.split_id,
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
