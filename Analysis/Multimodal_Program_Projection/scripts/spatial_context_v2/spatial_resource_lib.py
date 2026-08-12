#!/usr/bin/env python3
"""Fail-closed helpers for the spatial Resource candidate release."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


RESOURCE_RELEASE_ID = "spatial-resource-candidate-2026-08-11"
PROGRAM_RELEASE_ID = "program-context-v2-candidate-2026-08-07"
EXPECTED_REGISTRY_SHA256 = "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
EXPECTED_N_PROGRAMS = 117
EXPECTED_N_CONFIRMATORY = 2
OBSERVABLE_MIN_GENES = 8
OBSERVABLE_MIN_WEIGHT = 0.20
DETECTION_MIN_FRACTION = 0.01

ALLOWED_DATASET_GATES = {"pass", "source_dependent", "metadata_pending", "skipped", "dropped"}
ALLOWED_COVERAGE_STATES = {"observable", "partial", "untestable"}
ALLOWED_EVIDENCE_STATES = {
    "supported", "indeterminate", "untestable", "not_applicable", "source_dependent", "skipped", "dropped"
}
HMSMA_JOIN_COLUMNS = (
    "array_id", "donor_id", "array_relationship", "run_aggregation_id",
    "clinical_label", "nas_score", "fibrosis_score", "age", "sex", "bmi",
    "acquisition_batch", "histology_registration_status", "histology_image_id",
    "histology_orientation", "scrna_overlap_status", "maldi_sample_id",
    "maldi_registration_status",
)

REGISTRY_COLUMNS = (
    "dataset_id", "assay_id", "record_type", "decision", "paper_role", "dataset_gate",
    "source_dependence", "coverage_required", "biological_unit", "biological_unit_resolution",
    "n_biological", "technical_unit", "n_technical", "donor_join", "phenotype_join",
    "histology_join", "spatial_coordinate_join", "metabolite_join", "gene_axis_source", "detection_matrix",
    "abundance_unit", "lineage_context", "permitted_estimands", "prohibited_claims",
    "source_version", "license", "processing_provenance",
)

COVERAGE_COLUMNS = (
    "release_id", "program_release_id", "registry_sha256", "dataset_id", "assay_id",
    "program_uid", "program_label", "cell_type", "membership_sha256", "n_program_genes",
    "n_genes_on_axis", "n_genes_measured", "n_genes_detected", "retained_l1_weight", "detection_fraction",
    "gene_universe_denominator", "lineage_context", "coverage_status", "testability_reason",
    "dataset_gate", "biological_unit", "biological_unit_resolution", "source_dependence",
)

EFFECT_COLUMNS = (
    "release_id", "program_release_id", "registry_sha256", "dataset_id", "assay_id",
    "program_uid", "program_label", "membership_sha256", "biological_unit",
    "biological_unit_resolution", "n_biological", "technical_unit", "n_technical",
    "source_dependence", "dataset_gate", "estimand", "effect_unit", "estimate",
    "matched_null_sd", "pvalue", "qvalue", "multiplicity_family", "n_null_draws",
    "n_genes_measured", "retained_l1_weight", "equal_weight_sign_agree",
    "leave_top_weighted_gene_sign_agree", "within_source_result", "evidence_state", "testability_reason",
    "uncertainty_semantics", "source_release_id", "source_row_sha256",
)


class SpatialResourceError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_row_sha256(row: Mapping[str, object], columns: Sequence[str]) -> str:
    payload = "\t".join(str(row.get(column, "")) for column in columns) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def read_tsv(path: Path, required: Sequence[str] = ()) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise SpatialResourceError(f"missing required table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [field for field in required if field not in fields]
        if missing:
            raise SpatialResourceError(f"{path} missing columns: {missing}")
        return fields, [dict(row) for row in reader]


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def parse_bool(value: object, label: str) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise SpatialResourceError(f"invalid boolean for {label}: {value!r}")


def load_dataset_registry(path: Path) -> list[dict[str, str]]:
    _, rows = read_tsv(path, REGISTRY_COLUMNS)
    keys = [(row["dataset_id"], row["assay_id"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise SpatialResourceError("dataset registry contains duplicate dataset/assay IDs")
    for row in rows:
        if row["dataset_gate"] not in ALLOWED_DATASET_GATES:
            raise SpatialResourceError(f"invalid dataset_gate for {row['dataset_id']}: {row['dataset_gate']}")
        resolution = row["biological_unit_resolution"]
        if resolution == "resolved" and row["biological_unit"] == "unknown_public_biological_unit":
            raise SpatialResourceError(f"resolved biological unit cannot be unknown: {row['dataset_id']}")
        if row["dataset_gate"] in {"metadata_pending", "dropped", "skipped"} and not row["prohibited_claims"].strip():
            raise SpatialResourceError(f"closed source lacks prohibited-claim contract: {row['dataset_id']}")
        parse_bool(row["coverage_required"], f"coverage_required[{row['dataset_id']}]")
    return rows


def assert_authorized_inference(row: Mapping[str, str], inference_kind: str) -> None:
    if inference_kind == "donor" and row["biological_unit_resolution"] != "resolved":
        raise SpatialResourceError(f"donor inference refused for unresolved unit: {row['dataset_id']}")
    if inference_kind == "clinical" and row["phenotype_join"] != "resolved":
        raise SpatialResourceError(f"clinical inference refused for unresolved join: {row['dataset_id']}")
    if inference_kind == "histology" and row["histology_join"] != "resolved":
        raise SpatialResourceError(f"histology localization refused for unresolved join: {row['dataset_id']}")
    if inference_kind == "metabolite" and row["metabolite_join"] != "resolved":
        raise SpatialResourceError(f"metabolite localization refused for unresolved join: {row['dataset_id']}")


def validate_hmsma_metadata_join(
    rows: Sequence[Mapping[str, str]],
    expected_array_ids: Sequence[str],
    inference_kind: str = "donor",
) -> None:
    """Validate an authoritative HMSMA join without inferring missing identities.

    Repeated donors are permitted only when every repeated array explicitly
    declares its primary/section/technical relationship. Clinical, histology,
    and metabolite gates add requirements; they never impute missing fields.
    """
    if inference_kind not in {"donor", "clinical", "histology", "metabolite"}:
        raise SpatialResourceError(f"unknown HMSMA inference kind: {inference_kind}")
    if not rows:
        raise SpatialResourceError("HMSMA metadata join is empty")

    def text_value(row: Mapping[str, str], field: str) -> str:
        value = row.get(field, "")
        return "" if value is None else str(value).strip()

    for index, row in enumerate(rows, start=1):
        missing_columns = [column for column in HMSMA_JOIN_COLUMNS if column not in row]
        if missing_columns:
            raise SpatialResourceError(f"HMSMA metadata join row {index} missing columns: {missing_columns}")

    expected = {str(value).strip() for value in expected_array_ids if str(value).strip()}
    array_ids = [text_value(row, "array_id") for row in rows]
    if len(array_ids) != len(set(array_ids)):
        raise SpatialResourceError("HMSMA metadata join contains duplicated or ambiguous array IDs")
    if set(array_ids) != expected:
        missing = sorted(expected - set(array_ids))
        extra = sorted(set(array_ids) - expected)
        raise SpatialResourceError(f"HMSMA array join mismatch: missing={missing} extra={extra}")

    forbidden_identity_tokens = {"", "unknown", "ambiguous", "inferred", "provisional", "tbd", "na", "nan"}
    allowed_relationships = {"primary", "repeated_section", "technical_replicate"}
    donor_groups: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        donor_id = text_value(row, "donor_id")
        donor_lower = donor_id.lower()
        if donor_lower in forbidden_identity_tokens or any(
            token in donor_lower for token in ("unknown", "ambiguous", "inferred", "provisional", "tbd")
        ):
            raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks an authoritative donor ID")
        relationship = text_value(row, "array_relationship").lower()
        if relationship not in allowed_relationships:
            raise SpatialResourceError(f"HMSMA array {row['array_id']} has unresolved array relationship")
        if not text_value(row, "run_aggregation_id"):
            raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks a run aggregation ID")
        donor_groups[donor_id].append(row)

    for donor_id, donor_rows in donor_groups.items():
        if len(donor_rows) == 1:
            continue
        relationships = [text_value(row, "array_relationship").lower() for row in donor_rows]
        # Exactly one array anchors each repeated-donor group; row order has no
        # inferential meaning.
        if relationships.count("primary") != 1:
            raise SpatialResourceError(f"repeated HMSMA donor {donor_id} lacks exactly one primary array")
        if any(value not in {"primary", "repeated_section", "technical_replicate"} for value in relationships):
            raise SpatialResourceError(f"repeated HMSMA donor {donor_id} has ambiguous replicate relationships")

    if inference_kind in {"clinical", "histology", "metabolite"}:
        allowed_labels = {"control", "masld", "mash"}
        allowed_sex = {"female", "male", "intersex", "not_reported"}
        for row in rows:
            if text_value(row, "clinical_label").lower() not in allowed_labels:
                raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks an exact clinical label")
            if text_value(row, "sex").lower() not in allowed_sex:
                raise SpatialResourceError(f"HMSMA array {row['array_id']} has unresolved sex metadata")
            numeric_values: dict[str, float] = {}
            for field in ("nas_score", "fibrosis_score", "age", "bmi"):
                try:
                    value = float(text_value(row, field))
                except ValueError as exc:
                    raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks numeric {field}") from exc
                if not math.isfinite(value):
                    raise SpatialResourceError(f"HMSMA array {row['array_id']} has non-finite {field}")
                numeric_values[field] = value
            bounds = {
                "nas_score": (0, 8),
                "fibrosis_score": (0, 4),
                "age": (0, 120),
                "bmi": (10, 100),
            }
            for field, (lower, upper) in bounds.items():
                if not lower <= numeric_values[field] <= upper:
                    raise SpatialResourceError(f"HMSMA array {row['array_id']} has implausible {field}")
            if not text_value(row, "acquisition_batch"):
                raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks acquisition batch")
            if text_value(row, "scrna_overlap_status").lower() not in {"overlap", "no_overlap"}:
                raise SpatialResourceError(f"HMSMA array {row['array_id']} has unresolved spatial-to-scRNA overlap")

    if inference_kind == "histology":
        for row in rows:
            if text_value(row, "histology_registration_status").lower() != "registered":
                raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks registered histology")
            if not text_value(row, "histology_image_id") or not text_value(row, "histology_orientation"):
                raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks histology image/orientation")

    if inference_kind == "metabolite":
        for row in rows:
            if not text_value(row, "maldi_sample_id"):
                raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks a MALDI sample ID")
            if text_value(row, "maldi_registration_status").lower() != "registered":
                raise SpatialResourceError(f"HMSMA array {row['array_id']} lacks registered MALDI coordinates")


def load_frozen_programs(hotspot_root: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, float]]]:
    registry_path = hotspot_root / "program_registry_v2.tsv"
    if sha256_file(registry_path) != EXPECTED_REGISTRY_SHA256:
        raise SpatialResourceError("frozen program registry hash drift")
    _, rows = read_tsv(
        registry_path,
        ("program_uid", "cell_type", "module_name", "membership_sha256", "external_test_eligible"),
    )
    if len(rows) != EXPECTED_N_PROGRAMS or len({row["program_uid"] for row in rows}) != EXPECTED_N_PROGRAMS:
        raise SpatialResourceError(f"expected {EXPECTED_N_PROGRAMS} unique frozen programs, found {len(rows)}")
    selected = [row for row in rows if parse_bool(row["external_test_eligible"], row["program_uid"])]
    if len(selected) != EXPECTED_N_CONFIRMATORY:
        raise SpatialResourceError(f"expected exactly two confirmatory programs, found {len(selected)}")

    _, memberships = read_tsv(
        hotspot_root / "program_membership_v2.tsv",
        (
            "program_uid", "canonical_gene", "mapped_symbol", "mapped_symbol_status",
            "original_l1_weight", "membership_sha256",
        ),
    )
    weights: dict[str, dict[str, float]] = defaultdict(dict)
    source_weight_sums: dict[str, float] = defaultdict(float)
    membership_hashes: dict[str, set[str]] = defaultdict(set)
    for row in memberships:
        uid = row["program_uid"]
        if uid not in {item["program_uid"] for item in rows}:
            raise SpatialResourceError(f"membership references unknown program: {uid}")
        weight = float(row["original_l1_weight"])
        if not math.isfinite(weight) or weight <= 0:
            raise SpatialResourceError(f"invalid positive L1 weight for {uid}/{row['canonical_gene']}")
        source_weight_sums[uid] += weight
        if row["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed":
            gene = row["mapped_symbol"].strip()
            if not gene:
                raise SpatialResourceError(f"confirmed membership has no mapped symbol: {uid}")
            weights[uid][gene] = weights[uid].get(gene, 0.0) + weight
        membership_hashes[uid].add(row["membership_sha256"])
    by_uid = {row["program_uid"]: row for row in rows}
    for uid, row in by_uid.items():
        if membership_hashes[uid] != {row["membership_sha256"]}:
            raise SpatialResourceError(f"membership hash mismatch for {uid}")
        if not math.isclose(source_weight_sums[uid], 1.0, rel_tol=0, abs_tol=1e-9):
            raise SpatialResourceError(f"source L1 weights do not sum to one for {uid}")
        if not weights[uid]:
            raise SpatialResourceError(f"program has no confirmed assay-facing symbols: {uid}")
    return by_uid, dict(weights)


def coverage_state(n_measured: int, retained_weight: float) -> tuple[str, str]:
    if n_measured >= OBSERVABLE_MIN_GENES and retained_weight >= OBSERVABLE_MIN_WEIGHT:
        return "observable", "at_least_8_genes_and_20pct_frozen_l1_weight"
    if n_measured > 0:
        return "partial", "some_program_genes_measured_but_observability_threshold_failed"
    return "untestable", "no_program_genes_measured_or_gene_axis_unavailable"


def require_complete_coverage(rows: Sequence[Mapping[str, str]], registry: Sequence[Mapping[str, str]]) -> None:
    required = {(row["dataset_id"], row["assay_id"]) for row in registry if parse_bool(row["coverage_required"], row["dataset_id"])}
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in rows:
        groups[(row["dataset_id"], row["assay_id"])].append(row["program_uid"])
        if row["coverage_status"] not in ALLOWED_COVERAGE_STATES:
            raise SpatialResourceError(f"invalid coverage state: {row['coverage_status']}")
    if set(groups) != required:
        raise SpatialResourceError(f"coverage dataset groups differ: observed={sorted(groups)} expected={sorted(required)}")
    for key, uids in groups.items():
        if len(uids) != EXPECTED_N_PROGRAMS or len(set(uids)) != EXPECTED_N_PROGRAMS:
            raise SpatialResourceError(f"{key} does not have exactly {EXPECTED_N_PROGRAMS} unique coverage rows")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def dataframe_index_values(group) -> list[str]:
    """Read an AnnData dataframe index from an open h5py group."""
    index_name = group.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode("utf-8")
    if index_name not in group:
        raise SpatialResourceError(f"AnnData group has no index dataset {index_name!r}")
    values = group[index_name][...]
    output = []
    for value in values:
        output.append(value.decode("utf-8") if isinstance(value, bytes) else str(value))
    return output


def h5ad_column_values(node) -> list[str]:
    """Read a plain or categorical AnnData dataframe column."""
    if hasattr(node, "shape"):
        values = node[...]
    elif "codes" in node and "categories" in node:
        codes = node["codes"][...]
        categories = node["categories"][...]
        values = [categories[int(code)] if int(code) >= 0 else "" for code in codes]
    else:
        raise SpatialResourceError("unsupported AnnData dataframe column encoding")
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def h5ad_gene_sets(
    path: Path,
    matrix_path: str,
    min_detection_fraction: float = DETECTION_MIN_FRACTION,
) -> tuple[set[str], set[str]]:
    """Return gene-axis and adequately detected symbols without loading spots."""
    try:
        import h5py
        import numpy as np
    except ImportError as exc:
        raise SpatialResourceError("h5py and numpy are required for the coverage branch") from exc

    with h5py.File(path, "r") as handle:
        var = handle["var"]
        names = dataframe_index_values(var)
        if "gene_symbol" in var:
            symbols = h5ad_column_values(var["gene_symbol"])
        else:
            symbols = names
        symbols = [symbol.strip() for symbol in symbols]
        measured = {symbol for symbol in symbols if symbol}
        nonzero_counts = np.zeros(len(symbols), dtype=np.int64)
        if not matrix_path:
            raise SpatialResourceError(f"registered H5AD source lacks an explicit detection matrix: {path}")
        if matrix_path not in handle:
            raise SpatialResourceError(f"registered detection matrix {matrix_path!r} is absent: {path}")
        matrix = handle[matrix_path]
        if isinstance(matrix, h5py.Dataset):
            n_rows, n_columns = matrix.shape
            for start in range(0, n_rows, 2048):
                chunk = matrix[start : min(start + 2048, n_rows), :]
                nonzero_counts += np.count_nonzero(chunk, axis=0)
        else:
            encoding = matrix.attrs.get("encoding-type", "")
            if isinstance(encoding, bytes):
                encoding = encoding.decode("utf-8")
            data = matrix["data"]
            indices = matrix["indices"]
            indptr = matrix["indptr"]
            if encoding == "csr_matrix":
                n_rows = len(indptr) - 1
                for start in range(0, len(data), 5_000_000):
                    stop = min(start + 5_000_000, len(data))
                    values = data[start:stop]
                    cols = indices[start:stop]
                    present_columns = cols[np.asarray(values) != 0].astype(int)
                    nonzero_counts += np.bincount(present_columns, minlength=len(symbols))
            elif encoding == "csc_matrix":
                pointers = indptr[...]
                shape = matrix.attrs.get("shape")
                n_rows = int(shape[0]) if shape is not None else int(indices[...].max()) + 1
                for column in range(len(symbols)):
                    start, stop = int(pointers[column]), int(pointers[column + 1])
                    if stop > start:
                        nonzero_counts[column] = int(np.count_nonzero(data[start:stop]))
            else:
                raise SpatialResourceError(f"unsupported AnnData X encoding {encoding!r}: {path}")
        threshold = max(1, int(math.ceil(min_detection_fraction * n_rows)))
        detected_indices = np.flatnonzero(nonzero_counts >= threshold)
        detected = {symbols[index] for index in detected_indices if symbols[index]}
        return measured, detected
