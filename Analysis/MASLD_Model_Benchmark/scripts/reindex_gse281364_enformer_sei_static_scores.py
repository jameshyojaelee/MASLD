#!/usr/bin/env python3
"""Bind outcome-blind Enformer and Sei features to long-range-safe folds."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA = "masld-bench-gse281364-enformer-sei-static-reindex-v1"
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
ENFORMER_FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "hepg2_accessibility_sad",
    "hepg2_accessibility_sar",
    "liver_accessibility_sad",
    "liver_accessibility_sar",
    "all_accessibility_mean_sad",
    "all_accessibility_mean_sar",
)
SEI_FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "max_abs_sequence_class_index0",
    "signed_max_abs_score",
    "max_abs_score",
    "mean_signed_score",
)
TRACK_FIELDS = (
    "position",
    "canonical_index",
    "identifier",
    "description",
    "assay_family",
    "clip",
    "scale",
    "sum_stat",
)
BINDING_FIELDS = ROW_FIELDS + (
    "enformer_feature_index0",
    "sei_feature_index0",
    "legacy_outer_fold",
    "fold_reassigned",
)
BASE_FIELDS = (
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "outer_fold",
    "enformer_feature_index0",
    "sei_feature_index0",
    "legacy_outer_fold",
    "fold_reassigned",
    "enformer_hepg2_accessibility_sad",
    "enformer_hepg2_accessibility_sar",
    "enformer_liver_accessibility_sad",
    "enformer_liver_accessibility_sar",
    "enformer_all_accessibility_mean_sad",
    "enformer_all_accessibility_mean_sar",
    "sei_max_abs_sequence_class_index0",
    "sei_signed_max_abs_score",
    "sei_max_abs_score",
    "sei_mean_signed_score",
)
HEPG2_TRACKS = (27, 91, 234)
LIVER_TRACKS = (26, 448)


class StaticReindexError(RuntimeError):
    """Raised when an authority, feature matrix, or row binding differs."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise StaticReindexError(f"invalid JSON: {path}: {error}") from error
    if not isinstance(value, dict):
        raise StaticReindexError(f"JSON object required: {path}")
    return value


def safe_path(root: Path, value: object, *, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise StaticReindexError(f"{label} path is missing")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise StaticReindexError(f"unsafe {label} path")
    try:
        resolved = (root / relative).resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise StaticReindexError(f"{label} path is missing or escapes root") from error
    return resolved


def load_config(path: Path) -> dict[str, Any]:
    config = load_json(path)
    firewall = config.get("firewall", {})
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != "prepared_outcome_blind_static_reindex"
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or any(
            firewall.get(field) is not False
            for field in (
                "outcomes_may_be_read",
                "sealed_assets_may_be_read",
                "model_fit_may_run",
                "calibration_may_run",
                "models_may_be_ranked",
                "native_outputs_are_interchangeable",
                "global_frozen_census_may_change",
            )
        )
        or firewall.get("static_features_may_be_read") is not True
    ):
        raise StaticReindexError("top-level identity or outcome firewall differs")
    return config


def verify_tree(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    tree = safe_path(root, binding.get("tree_path"), label=label)
    expected = binding.get("artifacts_sha256")
    if not isinstance(expected, str) or file_sha256(tree / "ARTIFACTS.json") != expected:
        raise StaticReindexError(f"{label} ARTIFACTS identity differs")
    try:
        verify_frozen_tree(tree)
    except ArtifactError as error:
        raise StaticReindexError(f"{label} frozen tree differs: {error}") from error
    return tree


def verify_member(
    tree: Path,
    member: object,
    expected_sha256: object,
    *,
    label: str,
) -> Path:
    if not isinstance(member, str) or not isinstance(expected_sha256, str):
        raise StaticReindexError(f"{label} member contract is missing")
    relative = Path(member)
    if relative.is_absolute() or ".." in relative.parts:
        raise StaticReindexError(f"unsafe {label} member")
    path = tree / relative
    if path.is_symlink() or not path.is_file() or file_sha256(path) != expected_sha256:
        raise StaticReindexError(f"{label} member identity differs")
    manifest = load_json(tree / "ARTIFACTS.json")
    registered = {
        str(item.get("path")): item
        for item in manifest.get("artifacts", ())
        if isinstance(item, dict)
    }
    item = registered.get(member)
    if item is None or item.get("sha256") != expected_sha256 or item.get("size_bytes") != path.stat().st_size:
        raise StaticReindexError(f"{label} member is not bound by ARTIFACTS")
    return path


def verify_bounded_member(
    root: Path, binding: Mapping[str, Any], *, label: str
) -> Path:
    """Verify one small member of a frozen tree without rehashing huge weights."""
    tree = safe_path(root, binding.get("tree_path"), label=label)
    if file_sha256(tree / "ARTIFACTS.json") != binding.get("artifacts_sha256"):
        raise StaticReindexError(f"{label} ARTIFACTS identity differs")
    path = verify_member(
        tree,
        binding.get("member"),
        binding.get("member_sha256"),
        label=label,
    )
    if path.stat().st_size != binding.get("member_size_bytes"):
        raise StaticReindexError(f"{label} member size differs")
    return path


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise StaticReindexError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def collapse_row_universe(
    rows: Sequence[Mapping[str, str]],
    *,
    expected_seeds: Sequence[int],
    expected_contexts: Sequence[str],
) -> list[dict[str, str]]:
    seeds = tuple(str(value) for value in expected_seeds)
    contexts = tuple(expected_contexts)
    if not rows or any(tuple(row) != ROW_FIELDS for row in rows):
        raise StaticReindexError("row-universe schema differs")
    if len({(row["seed"], row["row_hash"]) for row in rows}) != len(rows):
        raise StaticReindexError("seed-by-row hashes are not unique")
    grouped: dict[str, list[Mapping[str, str]]] = {}
    order: list[str] = []
    for row in rows:
        element = row["element_id"]
        if element not in grouped:
            grouped[element] = []
            order.append(element)
        grouped[element].append(row)
    output: list[dict[str, str]] = []
    expected_pairs = {(seed, context) for seed in seeds for context in contexts}
    for element in order:
        values = grouped[element]
        pairs = {(row["seed"], row["assay_context_id"]) for row in values}
        invariant_fields = (
            "unit_hash",
            "block_hash",
            "stratum",
            "outer_fold",
            "study_id",
            "source_locus_group_id",
            "long_range_block_id",
        )
        if (
            len(values) != len(expected_pairs)
            or pairs != expected_pairs
            or any(len({row[field] for row in values}) != 1 for field in invariant_fields)
            or values[0]["stratum"] != "all"
            or values[0]["study_id"] != "gse281364"
            or values[0]["outer_fold"] not in {f"fold-{index}" for index in range(5)}
        ):
            raise StaticReindexError(f"row-universe topology differs for {element}")
        output.append({field: values[0][field] for field in ROW_FIELDS[2:]})
    if len({row["source_locus_group_id"] for row in output}) != len(output):
        raise StaticReindexError("source locus groups are not one-to-one with elements")
    return output


def indexed_source_rows(
    rows: Sequence[Mapping[str, str]],
    fields: Sequence[str],
    *,
    label: str,
) -> tuple[dict[str, tuple[int, dict[str, str]]], list[dict[str, str]]]:
    if any(tuple(row) != tuple(fields) for row in rows):
        raise StaticReindexError(f"{label} table schema differs")
    result: dict[str, tuple[int, dict[str, str]]] = {}
    ordered: list[dict[str, str]] = []
    for index, raw in enumerate(rows):
        row = dict(raw)
        element = row["element_id"]
        if element in result:
            raise StaticReindexError(f"duplicate {label} element")
        try:
            fold = int(row["outer_fold"])
        except ValueError as error:
            raise StaticReindexError(f"invalid {label} legacy fold") from error
        if fold not in range(5):
            raise StaticReindexError(f"invalid {label} legacy fold")
        result[element] = (index, row)
        ordered.append(row)
    return result, ordered


def validate_matrix(
    path: Path,
    *,
    shape: Sequence[int],
    dtype: str,
    label: str,
) -> np.ndarray:
    try:
        matrix = np.load(path, allow_pickle=False, mmap_mode="r")
    except (OSError, ValueError) as error:
        raise StaticReindexError(f"cannot load {label}: {error}") from error
    if tuple(matrix.shape) != tuple(shape) or str(matrix.dtype) != dtype or not np.isfinite(matrix).all():
        raise StaticReindexError(f"{label} matrix geometry or values differ")
    return matrix


def validate_sei_scalarization(
    matrix: np.ndarray, rows: Sequence[Mapping[str, str]]
) -> None:
    values = np.asarray(matrix, dtype=np.float64)
    indices = np.argmax(np.abs(values), axis=1)
    signed = values[np.arange(values.shape[0]), indices]
    maximum = np.abs(signed)
    means = values.mean(axis=1)
    for index, row in enumerate(rows):
        expected = (
            int(row["max_abs_sequence_class_index0"]),
            float(row["signed_max_abs_score"]),
            float(row["max_abs_score"]),
            float(row["mean_signed_score"]),
        )
        observed = (int(indices[index]), signed[index], maximum[index], means[index])
        if expected[0] != observed[0] or not all(
            math.isclose(float(left), float(right), rel_tol=2e-9, abs_tol=2e-10)
            for left, right in zip(expected[1:], observed[1:])
        ):
            raise StaticReindexError("Sei scalar summary does not match native class matrix")


def validate_enformer_scalarization(
    sad: np.ndarray,
    sar: np.ndarray,
    tracks: Sequence[Mapping[str, str]],
    rows: Sequence[Mapping[str, str]],
) -> None:
    access = [index for index, row in enumerate(tracks) if row["assay_family"] == "accessibility"]
    if len(access) != 684:
        raise StaticReindexError("Enformer accessibility-track census differs")
    # Reproduce the original writer's one-row-at-a-time float32 reductions.
    # A matrix-wide reduction can change cancellation rounding for signed SAD.
    def row_means(matrix: np.ndarray, indices: Sequence[int]) -> np.ndarray:
        return np.asarray(
            [matrix[index, indices].mean() for index in range(matrix.shape[0])],
            dtype=np.float32,
        )

    scalar = {
        "hepg2_accessibility_sad": row_means(sad, HEPG2_TRACKS),
        "hepg2_accessibility_sar": row_means(sar, HEPG2_TRACKS),
        "liver_accessibility_sad": row_means(sad, LIVER_TRACKS),
        "liver_accessibility_sar": row_means(sar, LIVER_TRACKS),
        "all_accessibility_mean_sad": row_means(sad, access),
        "all_accessibility_mean_sar": row_means(sar, access),
    }
    for index, row in enumerate(rows):
        for field, values in scalar.items():
            if not math.isclose(float(row[field]), float(values[index]), rel_tol=2e-6, abs_tol=2e-8):
                raise StaticReindexError(f"Enformer scalar summary differs: {field}")


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="", write_through=True) as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def reindex(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise StaticReindexError("output already exists")
    config = load_config(config_path)
    row_binding = config["row_universe"]
    row_tree = verify_tree(root, row_binding, label="row universe")
    row_path = verify_member(
        row_tree,
        row_binding["member"],
        row_binding["member_sha256"],
        label="row universe",
    )
    row_fields, row_rows = read_tsv(row_path)
    if row_fields != ROW_FIELDS or len(row_rows) != row_binding["rows"]:
        raise StaticReindexError("row-universe fields or row count differ")
    base = collapse_row_universe(
        row_rows,
        expected_seeds=row_binding["seeds"],
        expected_contexts=row_binding["assay_contexts"],
    )
    if (
        len(base) != row_binding["elements"]
        or len({row["long_range_block_id"] for row in base}) != row_binding["long_range_blocks"]
    ):
        raise StaticReindexError("row-universe element or long-range-block count differs")

    authorities: dict[str, dict[str, Any]] = {}
    for authority_id, authority in config["model_authorities"].items():
        path = safe_path(root, authority["path"], label=authority_id)
        if path.is_symlink() or not path.is_file() or file_sha256(path) != authority["sha256"]:
            raise StaticReindexError(f"{authority_id} identity differs")
        authorities[authority_id] = load_json(path)
    if (
        authorities["enformer_checkpoint"]["native_vs_converted_parity"]["parity_status"] != "UNRESOLVED"
        or authorities["enformer_exposure"]["checkpoint_findings"]["official_native_sonnet_human"]["exposure_state"]
        != "target_label_unexposed"
        or authorities["enformer_exposure"]["checkpoint_findings"]["registered_crested_enformer_human"]["exposure_state"]
        != "target_label_unexposed_with_conversion_uncertainty"
        or authorities["sei_exposure"]["checkpoint_findings"]["sei_zenodo_4906997"]["exposure_state"]
        != "target_label_unexposed"
        or authorities["sei_checkpoint"]["terms_audit"]["result"] != "restricted_comparator"
    ):
        raise StaticReindexError("license, parity, or exposure authority differs")

    enformer_config = config["enformer"]
    enformer_tree = verify_tree(root, enformer_config, label="Enformer predictions")
    enformer_receipt_path = verify_member(
        enformer_tree,
        enformer_config["receipt_member"],
        enformer_config["receipt_sha256"],
        label="Enformer receipt",
    )
    enformer_receipt = load_json(enformer_receipt_path)
    if (
        enformer_receipt.get("status") != "pass_outcome_blind_restricted_prediction"
        or enformer_receipt.get("registered_scientific_identity") != enformer_config["scientific_identity"]
        or enformer_receipt.get("elements") != row_binding["elements"]
        or enformer_receipt.get("outcomes_read")
        or enformer_receipt.get("sealed_outcomes_read")
        or enformer_receipt.get("model_fitted_or_adapted")
        or enformer_receipt.get("native_sonnet_parity_established")
        or enformer_receipt.get("champion_eligible")
    ):
        raise StaticReindexError("Enformer prediction receipt differs")
    enformer_table_path = verify_member(
        enformer_tree,
        enformer_config["table_member"],
        enformer_config["table_sha256"],
        label="Enformer prediction table",
    )
    enformer_fields, enformer_rows = read_tsv(enformer_table_path)
    if enformer_fields != ENFORMER_FIELDS or len(enformer_rows) != row_binding["elements"]:
        raise StaticReindexError("Enformer prediction table differs")
    enformer_index, enformer_ordered = indexed_source_rows(
        enformer_rows, ENFORMER_FIELDS, label="Enformer"
    )

    crosswalk_binding = enformer_config["track_crosswalk"]
    crosswalk_tree = verify_tree(root, crosswalk_binding, label="Enformer track crosswalk")
    crosswalk_path = verify_member(
        crosswalk_tree,
        crosswalk_binding["member"],
        crosswalk_binding["member_sha256"],
        label="Enformer track crosswalk",
    )
    track_fields, tracks = read_tsv(crosswalk_path)
    if (
        track_fields != TRACK_FIELDS
        or len(tracks) != 5313
        or [int(row["position"]) for row in tracks] != list(range(5313))
        or [int(row["canonical_index"]) for row in tracks] != list(range(5313))
    ):
        raise StaticReindexError("Enformer track crosswalk geometry differs")
    enformer_matrices: dict[str, np.ndarray] = {}
    matrix_contract: dict[str, dict[str, Any]] = {}
    for binding in enformer_config["feature_matrices"]:
        path = verify_member(
            enformer_tree,
            binding["member"],
            binding["sha256"],
            label=f"Enformer {binding['feature_id']}",
        )
        matrix = validate_matrix(
            path,
            shape=binding["shape"],
            dtype=binding["dtype"],
            label=f"Enformer {binding['feature_id']}",
        )
        enformer_matrices[binding["feature_id"]] = matrix
        matrix_contract[binding["feature_id"]] = {
            "tree_path": enformer_config["tree_path"],
            "artifacts_sha256": enformer_config["artifacts_sha256"],
            "member": binding["member"],
            "sha256": binding["sha256"],
            "shape": binding["shape"],
            "dtype": binding["dtype"],
        }
    if set(enformer_matrices) != {"forward_sad", "forward_sar", "rc_ensemble_sad", "rc_ensemble_sar"}:
        raise StaticReindexError("Enformer feature-matrix roster differs")
    validate_enformer_scalarization(
        enformer_matrices["rc_ensemble_sad"],
        enformer_matrices["rc_ensemble_sar"],
        tracks,
        enformer_ordered,
    )

    sei_config = config["sei"]
    sei_tree = verify_tree(root, sei_config, label="Sei predictions")
    sei_receipt_path = verify_member(
        sei_tree,
        sei_config["receipt_member"],
        sei_config["receipt_sha256"],
        label="Sei receipt",
    )
    sei_receipt = load_json(sei_receipt_path)
    if (
        sei_receipt.get("status") != "pass_outcome_blind_native_prediction"
        or sei_receipt.get("elements") != row_binding["elements"]
        or sei_receipt.get("native_sequence_classes") != 40
        or sei_receipt.get("outcomes_read")
        or sei_receipt.get("sealed_outcomes_read")
        or sei_receipt.get("model_fitted_or_adapted")
        or sei_receipt.get("champion_eligible")
    ):
        raise StaticReindexError("Sei prediction receipt differs")
    sei_table_path = verify_member(
        sei_tree,
        sei_config["table_member"],
        sei_config["table_sha256"],
        label="Sei prediction table",
    )
    sei_fields, sei_rows = read_tsv(sei_table_path)
    if sei_fields != SEI_FIELDS or len(sei_rows) != row_binding["elements"]:
        raise StaticReindexError("Sei prediction table differs")
    sei_index, sei_ordered = indexed_source_rows(sei_rows, SEI_FIELDS, label="Sei")
    sei_matrix_path = verify_member(
        sei_tree,
        sei_config["class_matrix_member"],
        sei_config["class_matrix_sha256"],
        label="Sei native class matrix",
    )
    sei_matrix = validate_matrix(
        sei_matrix_path,
        shape=sei_config["class_matrix_shape"],
        dtype=sei_config["class_matrix_dtype"],
        label="Sei native class",
    )
    validate_sei_scalarization(sei_matrix, sei_ordered)
    class_names_path = verify_bounded_member(
        root, sei_config["class_manifest"], label="Sei class manifest"
    )
    class_names = class_names_path.read_text(encoding="utf-8").splitlines()
    if len(class_names) != sei_config["class_manifest"]["rows"] or any(not name for name in class_names):
        raise StaticReindexError("Sei class-name manifest differs")

    base_by_element = {row["element_id"]: row for row in base}
    if set(base_by_element) != set(enformer_index) or set(base_by_element) != set(sei_index):
        raise StaticReindexError("static feature and long-range row universes differ")
    base_output: list[dict[str, object]] = []
    base_binding: dict[str, dict[str, object]] = {}
    reassigned = 0
    for row in base:
        element = row["element_id"]
        enformer_offset, enformer_row = enformer_index[element]
        sei_offset, sei_row = sei_index[element]
        if (
            enformer_row["outer_locus_sequence_group_id"] != row["source_locus_group_id"]
            or sei_row["outer_locus_sequence_group_id"] != row["source_locus_group_id"]
            or enformer_row["outer_fold"] != sei_row["outer_fold"]
        ):
            raise StaticReindexError("source group or legacy fold differs across models")
        changed = int(enformer_row["outer_fold"]) != int(row["outer_fold"].split("-")[1])
        reassigned += int(changed)
        binding = {
            "enformer_feature_index0": enformer_offset,
            "sei_feature_index0": sei_offset,
            "legacy_outer_fold": enformer_row["outer_fold"],
            "fold_reassigned": str(changed).lower(),
        }
        base_binding[element] = binding
        base_output.append(
            {
                "element_id": element,
                "source_locus_group_id": row["source_locus_group_id"],
                "long_range_block_id": row["long_range_block_id"],
                "outer_fold": row["outer_fold"],
                **binding,
                **{
                    f"enformer_{field}": enformer_row[field]
                    for field in enformer_config["scalar_features"]
                },
                "sei_max_abs_sequence_class_index0": sei_row["max_abs_sequence_class_index0"],
                **{
                    f"sei_{field}": sei_row[field]
                    for field in sei_config["scalar_features"]
                },
            }
        )
    if reassigned != row_binding["source_groups_reassigned_from_legacy_folds"]:
        raise StaticReindexError("legacy-to-long-range fold reassignment count differs")

    binding_rows = [
        {**row, **base_binding[row["element_id"]]}
        for row in row_rows
    ]
    output.mkdir(parents=True)
    write_gzip_tsv(output / "row_binding.tsv.gz", BINDING_FIELDS, binding_rows)
    write_tsv(output / "base_static_scores.tsv", BASE_FIELDS, base_output)
    write_tsv(output / "enformer_track_manifest.tsv", TRACK_FIELDS, tracks)
    write_tsv(
        output / "sei_sequence_class_manifest.tsv",
        ("class_index0", "class_name"),
        (
            {"class_index0": index, "class_name": name}
            for index, name in enumerate(class_names)
        ),
    )
    source_contract = {
        "schema_version": SCHEMA,
        "enformer": {
            "model_id": enformer_config["model_id"],
            "scientific_identity": enformer_config["scientific_identity"],
            "native_sonnet_parity_established": False,
            "open_champion_eligible": False,
            "feature_matrices": matrix_contract,
            "track_manifest_sha256": file_sha256(output / "enformer_track_manifest.tsv"),
            "track_count": 5313,
            "allowed_future_use": "restricted_development_head_fit_inside_long_range_outer_training_blocks_only",
        },
        "sei": {
            "model_id": sei_config["model_id"],
            "open_champion_eligible": False,
            "class_matrix": {
                "tree_path": sei_config["tree_path"],
                "artifacts_sha256": sei_config["artifacts_sha256"],
                "member": sei_config["class_matrix_member"],
                "sha256": sei_config["class_matrix_sha256"],
                "shape": sei_config["class_matrix_shape"],
                "dtype": sei_config["class_matrix_dtype"],
            },
            "class_manifest_sha256": file_sha256(output / "sei_sequence_class_manifest.tsv"),
            "class_count": 40,
            "allowed_future_use": "restricted_development_head_fit_inside_long_range_outer_training_blocks_only",
        },
        "outputs_are_interchangeable": False,
        "models_ranked": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "model_fit": False,
    }
    (output / "feature_source_contract.json").write_text(
        json.dumps(source_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt: dict[str, Any] = {
        "schema_version": SCHEMA,
        "status": "pass_outcome_blind_static_reindex",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "models": [enformer_config["model_id"], sei_config["model_id"]],
        "elements": len(base_output),
        "row_bindings": len(binding_rows),
        "long_range_blocks": len({row["long_range_block_id"] for row in base_output}),
        "outer_folds": 5,
        "fixed_seeds": row_binding["seeds"],
        "assay_contexts": row_binding["assay_contexts"],
        "source_groups_reassigned_from_legacy_folds": reassigned,
        "enformer_track_features_per_delta": 5313,
        "enformer_delta_feature_matrices": 4,
        "enformer_scalar_features": len(enformer_config["scalar_features"]),
        "sei_native_sequence_classes": 40,
        "sei_scalar_features": len(sei_config["scalar_features"]),
        "row_binding_sha256": file_sha256(output / "row_binding.tsv.gz"),
        "base_static_scores_sha256": file_sha256(output / "base_static_scores.tsv"),
        "feature_source_contract_sha256": file_sha256(output / "feature_source_contract.json"),
        "prediction_values_read": True,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "model_fit": False,
        "calibration_fit": False,
        "models_ranked": False,
        "native_outputs_are_interchangeable": False,
        "global_frozen_census_changed": False,
        "open_champion_eligible": False,
        "next_gate": "separately_freeze_five_seed_outer_training_only_heads_and_task_native_controls",
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = reindex(
        args.root.resolve(strict=True),
        args.config.resolve(strict=True),
        args.output,
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
