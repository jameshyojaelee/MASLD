#!/usr/bin/env python3
"""Reproject outcome-blind GSE281364 allele embeddings under 524-kb blocks."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import verify_frozen_tree
from scripts.gse281364_dna_lm_native_contract import (
    ALLELES,
    apply_projection,
    fit_projection,
    head_features,
)


SCHEMA = "masld-bench-gse281364-long-range-projected-embeddings-v1"
MODELS = ("caduceus", "dnabert2", "hyenadna", "nucleotide_transformer")
ROW_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
)


class LongRangeProjectionError(RuntimeError):
    """Raised when raw features cannot be safely moved to the long-range folds."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.tobytes())
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    try:
        config = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise LongRangeProjectionError(f"invalid projection config: {error}") from error
    if (
        config.get("schema_version")
        != "masld-bench-gse281364-sequence-context-refit-dag-v1"
        or config.get("outcome_access_authorized") is not False
        or config.get("prediction_generation_authorized") is not False
        or config.get("model_fit_authorized") is not False
    ):
        raise LongRangeProjectionError("projection firewall differs")
    return config


def safe_path(root: Path, relative_value: object, *, label: str) -> Path:
    if not isinstance(relative_value, str) or not relative_value:
        raise LongRangeProjectionError(f"{label} path is missing")
    relative = Path(relative_value)
    if relative.is_absolute() or ".." in relative.parts:
        raise LongRangeProjectionError(f"unsafe {label} path")
    try:
        resolved = (root / relative).resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise LongRangeProjectionError(f"{label} path is unavailable") from error
    if resolved.is_symlink():
        raise LongRangeProjectionError(f"{label} path is symlinked")
    return resolved


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise LongRangeProjectionError(f"cannot read TSV {path}: {error}") from error


def build_long_range_index(rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, str]]:
    if not rows or tuple(rows[0]) != ROW_FIELDS:
        raise LongRangeProjectionError("row-universe header differs")
    result: dict[str, dict[str, str]] = {}
    identities: set[tuple[str, str, str]] = set()
    for row in rows:
        if (
            row["study_id"] != "gse281364"
            or row["stratum"] != "all"
            or row["outer_fold"] not in {f"fold-{index}" for index in range(5)}
            or row["assay_context_id"] not in {"HepG2_control", "HepG2_PAOA"}
        ):
            raise LongRangeProjectionError("row-universe role or fold differs")
        identity = (row["seed"], row["row_hash"], row["assay_context_id"])
        if identity in identities:
            raise LongRangeProjectionError("duplicate row-universe identity")
        identities.add(identity)
        metadata = {
            "element_id": row["element_id"],
            "source_locus_group_id": row["source_locus_group_id"],
            "long_range_block_id": row["long_range_block_id"],
            "outer_fold": row["outer_fold"],
            "unit_hash": row["unit_hash"],
            "block_hash": row["block_hash"],
        }
        previous = result.setdefault(row["element_id"], metadata)
        if previous != metadata:
            raise LongRangeProjectionError("element metadata differs across seeds or contexts")
    if len(result) != 1033 or len(rows) != 10330:
        raise LongRangeProjectionError("row-universe denominator differs")
    return result


def _fixture_for_model(config: Mapping[str, Any], model_id: str) -> Mapping[str, Any]:
    matches = [
        fixture
        for fixture in config.get("projection_fixture", ())
        if model_id in fixture.get("model_ids", ())
    ]
    if len(matches) != 1:
        raise LongRangeProjectionError(f"{model_id} fixture binding differs")
    return matches[0]


def _inventory_for_model(config: Mapping[str, Any], model_id: str) -> Mapping[str, Any]:
    matches = [
        item
        for item in config.get("artifact_inventory", ())
        if item.get("model_id") == model_id
        and item.get("action") == "reindex_raw_then_recompute_projection_and_heads"
    ]
    if len(matches) != 1:
        raise LongRangeProjectionError(f"{model_id} raw inventory binding differs")
    return matches[0]


def project_model(
    *,
    root: Path,
    config: Mapping[str, Any],
    model_id: str,
    long_range: Mapping[str, Mapping[str, str]],
    output: Path,
) -> dict[str, Any]:
    inventory = _inventory_for_model(config, model_id)
    fixture = _fixture_for_model(config, model_id)
    source_root = safe_path(root, inventory["tree_path"], label=f"{model_id} source")
    fixture_root = safe_path(root, fixture["tree_path"], label=f"{model_id} fixture")
    if file_sha256(source_root / "ARTIFACTS.json") != inventory["artifacts_sha256"]:
        raise LongRangeProjectionError(f"{model_id} source ARTIFACTS hash differs")
    if file_sha256(fixture_root / "ARTIFACTS.json") != fixture["artifacts_sha256"]:
        raise LongRangeProjectionError(f"{model_id} fixture ARTIFACTS hash differs")
    verify_frozen_tree(source_root)
    verify_frozen_tree(fixture_root)
    fixture_path = fixture_root / str(fixture["manifest_member"])
    fixture_rows = read_tsv(fixture_path)
    raw_path = source_root / f"raw/{model_id}/allele_embeddings.npz"
    with np.load(raw_path, allow_pickle=False) as source:
        fixture_ids = source["fixture_ids"].astype(str)
        old_groups = source["outer_locus_sequence_group_ids"].astype(str)
        old_folds = source["outer_folds"].astype(np.int64)
        allele_order = tuple(source["allele_order"].astype(str).tolist())
        embeddings = source["embeddings"]
    if (
        len(fixture_rows) != 1033
        or fixture_ids.shape != (1033,)
        or old_groups.shape != (1033,)
        or old_folds.shape != (1033,)
        or allele_order != ALLELES
        or embeddings.ndim != 3
        or embeddings.shape[:2] != (1033, 4)
        or embeddings.shape[2] < 256
        or not np.isfinite(embeddings).all()
    ):
        raise LongRangeProjectionError(f"{model_id} raw embedding geometry differs")
    fixture_lookup = {row["fixture_id"]: row for row in fixture_rows}
    if len(fixture_lookup) != 1033 or set(fixture_ids) != set(fixture_lookup):
        raise LongRangeProjectionError(f"{model_id} fixture identity differs")
    element_ids: list[str] = []
    new_groups: list[str] = []
    block_ids: list[str] = []
    folds: list[int] = []
    for fixture_id in fixture_ids:
        fixture_row = fixture_lookup[fixture_id]
        element = fixture_row["element_id"]
        metadata = long_range.get(element)
        if metadata is None:
            raise LongRangeProjectionError(f"{model_id} element missing from row universe")
        if fixture_row["outer_locus_sequence_group_id"] != old_groups[len(element_ids)]:
            raise LongRangeProjectionError(f"{model_id} raw group does not match fixture")
        element_ids.append(element)
        new_groups.append(metadata["source_locus_group_id"])
        block_ids.append(metadata["long_range_block_id"])
        folds.append(int(metadata["outer_fold"].removeprefix("fold-")))
    fold_array = np.asarray(folds, dtype=np.int64)
    if set(fold_array.tolist()) != set(range(5)) or len(set(block_ids)) != 239:
        raise LongRangeProjectionError(f"{model_id} long-range fold coverage differs")
    reassigned = int(np.sum(old_folds != fold_array))
    if reassigned != 830:
        raise LongRangeProjectionError(f"{model_id} fold reassignment count differs")
    model_output = output / model_id
    model_output.mkdir(parents=True, exist_ok=False)
    fold_receipts = []
    for held_fold in range(5):
        parameters = fit_projection(embeddings, fold_array, held_fold, width=256)
        features = head_features(apply_projection(embeddings, parameters))
        fold_output = model_output / f"heldout_fold{held_fold}"
        fold_output.mkdir()
        np.savez_compressed(fold_output / "projection_parameters.npz", **parameters)
        np.savez_compressed(
            fold_output / "head_features.npz",
            fixture_ids=fixture_ids,
            element_ids=np.asarray(element_ids),
            source_locus_group_ids=np.asarray(new_groups),
            long_range_block_ids=np.asarray(block_ids),
            outer_folds=fold_array,
            feature_block_order=np.asarray(
                ["REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF"]
            ),
            features=features,
        )
        fold_receipts.append(
            {
                "held_out_fold": held_fold,
                "training_elements": int(np.sum(fold_array != held_fold)),
                "held_out_elements": int(np.sum(fold_array == held_fold)),
                "feature_sha256": array_sha256(features),
                "projection_components_sha256": array_sha256(parameters["components"]),
            }
        )
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_outcome_blind_long_range_projection",
        "dataset_id": "gse281364",
        "model_id": model_id,
        "elements": 1033,
        "source_locus_groups": 1033,
        "long_range_blocks": 239,
        "outer_folds": 5,
        "largest_receptive_field_buffer_bp": 524288,
        "source_groups_reassigned_from_original_fold_map": reassigned,
        "projection_width": 256,
        "head_input_width": 1024,
        "raw_embedding_sha256": file_sha256(raw_path),
        "raw_artifacts_sha256": inventory["artifacts_sha256"],
        "fixture_artifacts_sha256": fixture["artifacts_sha256"],
        "row_universe_sha256": config["row_binding"]["row_tsv_sha256"],
        "folds": fold_receipts,
        "outcomes_read": False,
        "prediction_values_read": False,
        "downstream_head_fit": False,
        "sealed_assets_read": False,
    }
    (model_output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config)
    if arguments.output.exists():
        raise LongRangeProjectionError("projection output exists")
    row_path = arguments.row_universe.resolve(strict=True)
    if file_sha256(row_path) != config["row_binding"]["row_tsv_sha256"]:
        raise LongRangeProjectionError("row-universe checksum differs")
    rows = read_tsv(row_path)
    long_range = build_long_range_index(rows)
    requested = tuple(value for value in arguments.model_roster.split(",") if value)
    if not requested or len(set(requested)) != len(requested) or any(
        model not in MODELS for model in requested
    ):
        raise LongRangeProjectionError("model roster differs")
    arguments.output.mkdir(parents=True, mode=0o750)
    model_receipts = {
        model: project_model(
            root=root,
            config=config,
            model_id=model,
            long_range=long_range,
            output=arguments.output,
        )
        for model in requested
    }
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_outcome_blind_long_range_projection",
        "dataset_id": "gse281364",
        "models": list(requested),
        "model_receipts": model_receipts,
        "elements": 1033,
        "long_range_blocks": 239,
        "fixed_future_head_seeds": config["row_binding"]["fixed_seeds"],
        "outcomes_read": False,
        "prediction_values_read": False,
        "downstream_head_fit": False,
        "sealed_assets_read": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--row-universe", type=Path, required=True)
    parser.add_argument("--model-roster", required=True)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
