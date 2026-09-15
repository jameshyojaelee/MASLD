#!/usr/bin/env python3
"""Build a label-free GSE281367 fixed-window ATAC exchange axis."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, write_json_exclusive


SCHEMA = "masld-bench-gse281367-label-free-atac-exchange-v1"
EXCHANGE_ID = "gse296875_to_gse281367_fixed_window_atac_transport_20260825"
PRIMARY_LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
CANONICAL_LINEAGES = (*PRIMARY_LINEAGES, "t_cell")
LABEL_CONTRACT = {
    "B_Cell": ("b_cell", "secondary"),
    "Cholangiocyte": ("cholangiocyte", "primary"),
    "Endothelial": ("endothelial_cell", "secondary"),
    "Hepatocyte": ("hepatocyte", "primary"),
    "Kupffer_Cell": ("macrophage", "primary"),
    "LSEC": ("endothelial_cell", "secondary"),
    "Macrophage": ("macrophage", "primary"),
    "NK_T_Cell": ("t_nk_cell", "secondary"),
    "Plasma_Cell": ("b_cell", "secondary"),
    "Stellate_Cell": ("fibroblast", "primary"),
}
MEMBERSHIP_FIELDS = (
    "dataset_id",
    "donor_id",
    "raw_barcode",
    "lineage_id",
    "analysis_role",
    "outer_fold",
)
WINDOW_FIELDS = ("window_index", "window_id", "role", "contig", "start", "end")
FORBIDDEN_TOKENS = ("condition", "diagnosis", "mash", "masld", "fibrosis", "nas", "sex", "age", "bmi")


class LabelFreeATACExchangeError(ValueError):
    """Raised when a label-free exchange authority does not meet its requirements."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise LabelFreeATACExchangeError("exchange contract is not an object")
    return value


def _read_h5_column(group: Any, key: str) -> list[str]:
    import numpy as np

    if key not in group:
        raise LabelFreeATACExchangeError(f"H5AD obs column missing: {key}")
    value = group[key]
    if hasattr(value, "shape"):
        values = value[...]
        return [item.decode() if isinstance(item, bytes) else str(item) for item in values]
    if "categories" not in value or "codes" not in value:
        raise LabelFreeATACExchangeError(f"unsupported H5AD obs encoding: {key}")
    categories = [
        item.decode() if isinstance(item, bytes) else str(item)
        for item in value["categories"][...]
    ]
    codes = np.asarray(value["codes"][...], dtype=int)
    return ["" if code < 0 else categories[code] for code in codes]


def read_label_free_obs(path: Path, fields: Sequence[str]) -> dict[str, list[str]]:
    """Read only named identifier/annotation axes; never access condition."""

    if any(any(token in field.lower() for token in FORBIDDEN_TOKENS) for field in fields):
        raise LabelFreeATACExchangeError("forbidden phenotype-like H5AD field requested")
    import h5py

    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        index_key = obs.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        result = {"source_cell_id": _read_h5_column(obs, str(index_key))}
        result.update({field: _read_h5_column(obs, field) for field in fields})
    if len({len(values) for values in result.values()}) != 1:
        raise LabelFreeATACExchangeError("H5AD label-free axes differ")
    return result


def resolve_raw_barcodes(merged_ids: Sequence[str], raw_ids: Sequence[str]) -> list[str]:
    raw = set(raw_ids)
    if len(raw) != len(raw_ids):
        raise LabelFreeATACExchangeError("per-donor raw barcodes are duplicated")
    used: set[str] = set()
    result: list[str] = []
    for merged in merged_ids:
        candidate = merged
        matches: list[str] = []
        while candidate:
            if candidate in raw:
                matches.append(candidate)
            suffix = re.fullmatch(r"(.+)-([0-9]+)", candidate)
            if suffix is None:
                break
            candidate = suffix.group(1)
        available = [item for item in matches if item not in used]
        if len(available) != 1:
            raise LabelFreeATACExchangeError("merged/raw barcode join is not bijective")
        used.add(available[0])
        result.append(available[0])
    if used != raw:
        raise LabelFreeATACExchangeError("merged and per-donor cell universes differ")
    return result


def fold_index(dataset_id: str, donor_id: str, seed: int = 20260821) -> int:
    digest = sha256(f"{seed}\0{dataset_id}\0{donor_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % 5


def _read_inventory(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise LabelFreeATACExchangeError(f"inventory fields differ: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise LabelFreeATACExchangeError("inventory is empty")
    return rows


def _verify_tree(root: Path, expected: str, artifact_class: str) -> None:
    manifest_path = root / "ARTIFACTS.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        not isinstance(manifest, dict)
        or _digest(manifest_path) != expected
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise LabelFreeATACExchangeError(f"artifact authority differs: {root}")


def validate_contract(root: Path, contract: Mapping[str, Any]) -> dict[str, Path]:
    if contract.get("schema_version") != SCHEMA or contract.get("exchange_id") != EXCHANGE_ID:
        raise LabelFreeATACExchangeError("exchange identity differs")
    for field in (
        "condition_labels_available_to_model_jobs",
        "condition_labels_used_for_axis_selection",
        "condition_labels_used_for_model_selection",
        "phenotype_values_available_to_model_jobs",
        "development_outcomes_available_to_model_jobs",
        "sealed_data_available",
        "external_or_sealed_evaluation",
        "champion_claim_allowed",
        "universal_claim_allowed",
    ):
        if contract.get(field) is not False:
            raise LabelFreeATACExchangeError(f"label or claim firewall opened: {field}")
    resolved: dict[str, Path] = {}
    authorities = (
        ("parent_readiness", "gse281367_atac_transport_production_readiness_audit"),
        ("activation_authority", "gse281367_activation_readiness_audit"),
        ("source_authority", "atac_transport_source_admission_audit"),
        ("window_authority", "atac_transport_evaluator_outcomes"),
        ("census_authority", "observed_multiome_additive_census_revision"),
    )
    for section, artifact_class in authorities:
        spec = contract.get(section)
        if not isinstance(spec, dict):
            raise LabelFreeATACExchangeError(f"missing authority: {section}")
        path = (root / str(spec.get("path", ""))).resolve(strict=True)
        path.relative_to(root)
        _verify_tree(path, str(spec.get("artifacts_sha256", "")), artifact_class)
        resolved[section] = path
    source = contract["source_authority"]
    window = contract["window_authority"]
    source_inventory = resolved["source_authority"] / str(source["inventory_path"])
    donor_inventory_root = (root / str(source["per_donor_inventory_authority"])).resolve(strict=True)
    _verify_tree(
        donor_inventory_root,
        str(source["per_donor_inventory_artifacts_sha256"]),
        "atac_transport_cell_membership",
    )
    donor_inventory = donor_inventory_root / str(source["per_donor_inventory_path"])
    windows = resolved["window_authority"] / str(window["coordinate_only_path"])
    for path, expected in (
        (source_inventory, source["inventory_sha256"]),
        (donor_inventory, source["per_donor_inventory_sha256"]),
        (windows, window["coordinate_only_sha256"]),
    ):
        if _digest(path) != expected:
            raise LabelFreeATACExchangeError(f"bound file differs: {path}")
    resolved.update(
        source_inventory=source_inventory,
        donor_inventory=donor_inventory,
        windows=windows,
    )
    axis = contract.get("label_free_axis", {})
    exchange = contract.get("window_exchange", {})
    rotations = contract.get("rotation", [])
    if (
        axis.get("donor_count") != 12
        or tuple(axis.get("primary_lineages", ())) != PRIMARY_LINEAGES
        or tuple(axis.get("canonical_lineages", ())) != CANONICAL_LINEAGES
        or axis.get("minimum_cells") != 50
        or axis.get("outer_fold_uses_condition") is not False
        or axis.get("missing_as_zero") is not False
        or exchange.get("window_count") != 32000
        or exchange.get("largest_admitted_receptive_field_bp") != 524288
        or exchange.get("target_role_contigs_removed_from_query_atac") is not True
        or [(row.get("observed_atac_context_role"), row.get("scored_target_role")) for row in rotations]
        != [("valid", "test"), ("test", "valid")]
    ):
        raise LabelFreeATACExchangeError("axis, mask, or rotation contract differs")
    if contract.get("claim_boundary", {}).get("gse281367_is_non_champion") is not True:
        raise LabelFreeATACExchangeError("GSE281367 non-champion boundary differs")
    implementation = contract.get("implementation")
    if not isinstance(implementation, dict):
        raise LabelFreeATACExchangeError("implementation authority is absent")
    for prefix in ("builder", "builder_test", "commit", "commit_test", "scorer", "scorer_test"):
        path = (root / str(implementation.get(f"{prefix}_path", ""))).resolve(strict=True)
        path.relative_to(root)
        if path.is_symlink() or not path.is_file() or _digest(path) != implementation.get(f"{prefix}_sha256"):
            raise LabelFreeATACExchangeError(f"implementation authority differs: {prefix}")
    family = contract.get("family_eligibility")
    if not isinstance(family, list) or len(family) != 25:
        raise LabelFreeATACExchangeError("family-native eligibility roster differs")
    model_ids = [row.get("model_id") for row in family]
    if len(set(model_ids)) != len(model_ids) or set(model_ids) != {
        "alphagenome", "borzoi_ensemble", "bpnet", "chrombpnet", "context_borzoi",
        "corgi_plus", "corgi_regular", "dna_language_models", "enformer", "epiagent",
        "epibert", "epcotv2", "get", "lsi", "multivi", "observed_atac_glm",
        "observed_atac_only", "peakvi", "scbasset", "scbasset_observed_cell_embedding",
        "scooby_released", "sei", "sequence_cnn_control", "sequence_transformer_control",
        "seurat_wnn",
    }:
        raise LabelFreeATACExchangeError("family-native model identities differ")
    missingness = contract.get("missingness_contract", {})
    if (
        missingness.get("missing_prediction_value") != "nan"
        or missingness.get("missing_prediction_value_may_be_zero") is not False
        or sorted(
            value for key, value in missingness.items() if key.endswith("_code")
        )
        != list(range(8))
    ):
        raise LabelFreeATACExchangeError("missingness encoding differs")
    decision = contract.get("decision", {})
    if (
        decision.get("label_free_exchange_contract_frozen") is not True
        or decision.get("condition_blind_scorer_registered") is not True
        or decision.get("prediction_hash_commit_registered") is not True
        or decision.get("biological_scoring_authorized_before_prediction_commit") is not False
        or decision.get("outcome_selected_model_repair_allowed") is not False
    ):
        raise LabelFreeATACExchangeError("exchange decision differs")
    return resolved


def _hash_source(record: Mapping[str, str]) -> tuple[Path, str]:
    path = Path(record["path"]).resolve(strict=True)
    if path.is_symlink() or path.stat().st_size != int(record["size_bytes"]):
        raise LabelFreeATACExchangeError(f"source path or size differs: {path}")
    observed = _digest(path)
    if observed != record["sha256"]:
        raise LabelFreeATACExchangeError(f"source hash differs: {path}")
    return path, observed


def build_membership(
    merged_path: Path,
    donor_records: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    merged = read_label_free_obs(merged_path, ("donor_id", "cell_type"))
    donor_ids = sorted(set(merged["donor_id"]))
    if donor_ids != [f"Z{index:02d}" for index in range(1, 13)]:
        raise LabelFreeATACExchangeError("target donor axis differs")
    by_donor = {Path(record["path"]).stem: record for record in donor_records}
    if set(by_donor) != set(donor_ids):
        raise LabelFreeATACExchangeError("per-donor H5AD roster differs")
    rows: list[dict[str, Any]] = []
    counts: Counter[tuple[str, str]] = Counter()
    for donor in donor_ids:
        indexes = [index for index, value in enumerate(merged["donor_id"]) if value == donor]
        if not indexes:
            raise LabelFreeATACExchangeError(f"donor absent from merged H5AD: {donor}")
        raw = read_label_free_obs(Path(by_donor[donor]["path"]), ())
        merged_ids = [merged["source_cell_id"][index] for index in indexes]
        raw_barcodes = resolve_raw_barcodes(merged_ids, raw["source_cell_id"])
        for index, barcode in zip(indexes, raw_barcodes, strict=True):
            source_label = merged["cell_type"][index]
            if source_label not in LABEL_CONTRACT:
                raise LabelFreeATACExchangeError(f"unregistered cell label: {source_label}")
            lineage, role = LABEL_CONTRACT[source_label]
            rows.append(
                {
                    "dataset_id": "gse281367",
                    "donor_id": donor,
                    "raw_barcode": barcode,
                    "lineage_id": lineage,
                    "analysis_role": role,
                    "outer_fold": fold_index("gse281367", donor),
                }
            )
            if role == "primary":
                counts[(donor, lineage)] += 1
    if len(rows) != 226224 or len({(row["donor_id"], row["raw_barcode"]) for row in rows}) != len(rows):
        raise LabelFreeATACExchangeError("project-frozen cell census or barcode uniqueness differs")
    axis_rows: list[dict[str, Any]] = []
    for donor_index, donor in enumerate(donor_ids):
        for lineage_index, lineage in enumerate(CANONICAL_LINEAGES):
            cells = counts[(donor, lineage)] if lineage in PRIMARY_LINEAGES else None
            if lineage == "t_cell":
                state, eligible = "not_applicable", False
            elif cells == 0:
                state, eligible = "structurally_missing", False
            elif cells < 50:
                state, eligible = "below_qc", False
            else:
                state, eligible = "observed", True
            axis_rows.append(
                {
                    "dataset_id": "gse281367",
                    "donor_index": donor_index,
                    "donor_id": donor,
                    "lineage_index": lineage_index,
                    "lineage_id": lineage,
                    "outer_fold": fold_index("gse281367", donor),
                    "cells": "" if cells is None else cells,
                    "evidence_state": state,
                    "eligible_min_50_cells": str(eligible).lower(),
                }
            )
    return rows, axis_rows


def read_windows(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != WINDOW_FIELDS:
            raise LabelFreeATACExchangeError("window fields differ")
        rows = [dict(row) for row in reader]
    if len(rows) != 32000 or [int(row["window_index"]) for row in rows] != list(range(32000)):
        raise LabelFreeATACExchangeError("window index axis differs")
    role_counts = Counter(row["role"] for row in rows)
    if role_counts != {"valid": 16000, "test": 16000}:
        raise LabelFreeATACExchangeError("window role counts differ")
    contigs = {role: {row["contig"] for row in rows if row["role"] == role} for role in role_counts}
    if contigs["valid"] & contigs["test"]:
        raise LabelFreeATACExchangeError("context and target chromosome roles overlap")
    for row in rows:
        start, end = int(row["start"]), int(row["end"])
        if start < 0 or end - start != 1000:
            raise LabelFreeATACExchangeError("window geometry differs")
    return rows


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_membership(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=MEMBERSHIP_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def build(root: Path, contract_path: Path, output: Path, hash_workers: int) -> dict[str, Any]:
    if output.exists() or not 1 <= hash_workers <= 4:
        raise LabelFreeATACExchangeError("output or hash-worker contract differs")
    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    contract_path.relative_to(root)
    contract = _load_toml(contract_path)
    resolved = validate_contract(root, contract)
    source_rows = _read_inventory(resolved["source_inventory"], ("path", "size_bytes", "sha256"))
    merged_records = [row for row in source_rows if row["path"].endswith(contract["source_authority"]["merged_h5ad_suffix"])]
    donor_rows = [
        row
        for row in _read_inventory(
            resolved["donor_inventory"],
            ("dataset_id", "path", "size_bytes", "sha256"),
        )
        if row["dataset_id"] == "gse281367"
    ]
    if len(merged_records) != 1 or len(donor_rows) != 12:
        raise LabelFreeATACExchangeError("merged or per-donor H5AD roster differs")
    with ThreadPoolExecutor(max_workers=hash_workers) as executor:
        verified = list(executor.map(_hash_source, [*merged_records, *donor_rows]))
    membership, axis = build_membership(verified[0][0], donor_rows)
    windows = read_windows(resolved["windows"])

    output.mkdir(parents=True, mode=0o750)
    _write_membership(output / "label_free_cell_membership.tsv.gz", membership)
    _write_tsv(
        output / "donor_lineage_axis.tsv",
        (
            "dataset_id", "donor_index", "donor_id", "lineage_index", "lineage_id",
            "outer_fold", "cells", "evidence_state", "eligible_min_50_cells",
        ),
        axis,
    )
    _write_tsv(output / "windows.tsv", WINDOW_FIELDS, windows)
    _write_tsv(
        output / "rotations.tsv",
        ("rotation_id", "observed_atac_context_role", "scored_target_role"),
        contract["rotation"],
    )
    _write_tsv(
        output / "family_eligibility.tsv",
        ("model_id", "input_regime", "output_family", "eligibility", "query_atac_allowed"),
        contract["family_eligibility"],
    )
    write_json_exclusive(output / "missingness_contract.json", contract["missingness_contract"])
    write_json_exclusive(output / "overlap_contract.json", contract["overlap_contract"])
    summary = {
        "schema_version": "masld-bench-gse281367-label-free-atac-exchange-summary-v1",
        "exchange_id": EXCHANGE_ID,
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "donors": 12,
        "membership_cells": len(membership),
        "donor_lineage_rows": len(axis),
        "eligible_donor_lineage_rows": sum(row["eligible_min_50_cells"] == "true" for row in axis),
        "windows": len(windows),
        "rotations": 2,
        "condition_h5_dataset_opened": False,
        "condition_values_read": False,
        "condition_values_used": False,
        "biological_outcome_arrays_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
        "gse281367_champion_claim_allowed": False,
        "source_hashes_reverified": len(verified),
    }
    write_json_exclusive(output / "exchange_summary.json", summary)
    (output / "source.sha256").write_text(
        f"{_digest(contract_path)}  {contract_path}\n{_digest(Path(__file__))}  {Path(__file__).resolve()}\n",
        encoding="utf-8",
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse281367_label_free_atac_exchange_axis",
            "exchange_id": EXCHANGE_ID,
            "source_dataset_id": "gse296875",
            "target_dataset_id": "gse281367",
            "condition_values_read": False,
            "biological_outcome_arrays_read": False,
            "metrics_calculated": False,
            "champion_claim_allowed": False,
            "status": "passed",
        },
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hash-workers", type=int, default=2)
    arguments = parser.parse_args()
    print(json.dumps(build(arguments.root, arguments.contract, arguments.output, arguments.hash_workers), sort_keys=True))


if __name__ == "__main__":
    main()
