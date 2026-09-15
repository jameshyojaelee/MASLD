#!/usr/bin/env python3
"""Build five target-value-free GSE296875 observed-multiome model inputs."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)
from masld_bench.observed_multiome_materialization import (
    ObservedMultiomeMaterializationError,
    decode_strings,
    donor_lineage_aggregation,
    read_full_csr,
    salted_hash,
    selected_value_positions,
    selectively_read_csr_columns,
    write_csr,
    write_string_dataset,
)


SCHEMA = "masld-bench-observed-multiome-model-input-materialization-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_model_inputs_20260825"


class ModelInputMaterializationError(ValueError):
    """Raised when a model-visible output file can contain held target values."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ModelInputMaterializationError("JSON object required")
    return value


def _read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ModelInputMaterializationError(f"TSV schema differs: {path.name}")
        return [dict(row) for row in reader]


def _bound_tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(
        root / str(record.get("path", "")), label=label
    )
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise ModelInputMaterializationError(f"{label} artifact drifted")
    return path


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("campaign_id") != CAMPAIGN_ID
        or config.get("dataset_id") != "gse296875"
        or config.get("dataset_view_id") != "gse296875_rna_atac_smoke_1000_v1"
        or config.get("stage") != "smoke"
    ):
        raise ModelInputMaterializationError("campaign identity differs")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {
        "materialization_contract",
        "mask_plan",
    }:
        raise ModelInputMaterializationError("parent roster differs")
    resolved = {
        "materialization_contract": _bound_tree(
            root, parents["materialization_contract"], "materialization contract"
        ),
        "mask_plan": _bound_tree(root, parents["mask_plan"], "mask plan"),
    }
    h5_record = config.get("input_h5")
    if not isinstance(h5_record, dict):
        raise ModelInputMaterializationError("input HDF5 record is absent")
    h5_lexical = reject_symlink_components(
        root / str(h5_record.get("path", "")), label="input HDF5"
    )
    h5_path = h5_lexical.resolve(strict=True)
    h5_path.relative_to(root)
    if not h5_path.is_file() or _digest(h5_path) != h5_record.get("sha256"):
        raise ModelInputMaterializationError("input HDF5 drifted")
    resolved["input_h5"] = h5_path
    if config.get("expected") != {
        "nuclei": 1000,
        "donors": 39,
        "donor_lineage_units": 195,
        "rna_features": 36601,
        "atac_features": 306706,
        "genomic_folds": 5,
        "input_peaks_per_fold": 2000,
        "target_peaks_per_fold": 1000,
    }:
        raise ModelInputMaterializationError("expected census differs")
    namespaces = config.get("hash_namespaces")
    if namespaces != {
        "donor": "gse296875:observed_multiome:donor:v1",
        "row": "gse296875:observed_multiome:row:v1",
        "peak": "gse296875:observed_multiome:peak:v1",
    }:
        raise ModelInputMaterializationError("hash namespaces differ")
    artifact = config.get("model_artifact")
    if (
        not isinstance(artifact, dict)
        or artifact.get("one_separately_frozen_child_per_genomic_fold") is not True
        or artifact.get("rna_values")
        != "raw_integer_donor_by_lineage_sum_all_features"
        or artifact.get("observed_atac_values")
        != "raw_integer_donor_by_lineage_sum_selected_input_peaks_only"
        or artifact.get("target_values") != "forbidden"
        or artifact.get("target_identifiers_and_coordinates")
        != "required_without_values"
        or artifact.get("normalization") != "none"
        or artifact.get("feature_selection") != "none"
        or artifact.get("raw_donor_ids_exported") is not False
        or artifact.get("explicit_rna_and_atac_observed_masks") is not True
    ):
        raise ModelInputMaterializationError("model artifact contract differs")
    selective = config.get("selective_read")
    if (
        not isinstance(selective, dict)
        or selective.get("atac_indices_and_indptr_read") is not True
        or selective.get("atac_data_values_read_by_monotone_selected_position_index")
        is not True
        or selective.get("full_atac_data_array_read") is not False
        or selective.get("same_fold_selected_input_and_target_value_positions_disjoint")
        is not True
    ):
        raise ModelInputMaterializationError("selective read contract differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise ModelInputMaterializationError("model materialization firewall is open")
    return resolved


def _rows_and_masks(
    mask_root: Path, cell_donors: Sequence[str], cell_lineages: Sequence[str], namespaces: Mapping[str, str]
) -> tuple[list[dict[str, str]], Any, list[str], list[str]]:
    planned = _read_tsv(
        mask_root / "donor_lineage_units.tsv",
        ("donor_hash", "lineage", "donor_fold", "nuclei"),
    )
    cell_unit_keys = [
        (salted_hash(namespaces["donor"], donor), lineage)
        for donor, lineage in zip(cell_donors, cell_lineages, strict=True)
    ]
    aggregation, _ = donor_lineage_aggregation(
        cell_unit_keys=cell_unit_keys, planned_units=planned
    )
    row_hashes = [
        salted_hash(namespaces["row"], row["donor_hash"], row["lineage"])
        for row in planned
    ]
    raw_donor_ids = sorted(set(cell_donors))
    if len(set(row_hashes)) != len(planned) or len(raw_donor_ids) != 39:
        raise ModelInputMaterializationError("row hashes or donor census differ")
    return planned, aggregation, row_hashes, raw_donor_ids


def _peak_index(
    peak_ids: Sequence[str], chromosomes: Sequence[str], starts: Sequence[int], ends: Sequence[int]
) -> dict[str, tuple[int, str, int, int]]:
    if not (
        len(peak_ids) == len(chromosomes) == len(starts) == len(ends)
        and len(set(peak_ids)) == len(peak_ids)
    ):
        raise ModelInputMaterializationError("peak axes differ")
    return {
        peak_id: (index, chromosome, int(start), int(end))
        for index, (peak_id, chromosome, start, end) in enumerate(
            zip(peak_ids, chromosomes, starts, ends, strict=True)
        )
    }


def _resolve_mask_rows(
    *,
    peak_index: Mapping[str, tuple[int, str, int, int]],
    rows: Sequence[Mapping[str, str]],
    held_fold: int,
    peak_namespace: str,
    kind: str,
) -> tuple[list[int], list[dict[str, str]]]:
    expected_hash_field = "input_hash" if kind == "input" else "target_hash"
    by_hash: dict[str, tuple[int, str, int, int]] = {}
    for peak_id, metadata in peak_index.items():
        hash_value = (
            salted_hash(peak_namespace, "input", held_fold, peak_id)
            if kind == "input"
            else salted_hash(peak_namespace, "target", peak_id)
        )
        by_hash[hash_value] = metadata
    columns: list[int] = []
    output_rows: list[dict[str, str]] = []
    for row in rows:
        metadata = by_hash.get(row[expected_hash_field])
        if metadata is None:
            raise ModelInputMaterializationError(f"{kind} hash is absent from HDF5")
        column, chromosome, start, end = metadata
        if (
            chromosome != row["chromosome"]
            or start != int(row["bed_start_0based"])
            or end != int(row["bed_end_half_open"])
        ):
            raise ModelInputMaterializationError(f"{kind} coordinate differs")
        columns.append(column)
        output_rows.append(dict(row))
    if len(columns) != len(set(columns)):
        raise ModelInputMaterializationError(f"{kind} columns are duplicated")
    return columns, output_rows


def _write_fold_h5(
    *,
    path: Path,
    held_fold: int,
    planned_units: Sequence[Mapping[str, str]],
    row_hashes: Sequence[str],
    rna_counts: Any,
    ensembl_ids: Sequence[str],
    gene_names: Sequence[str],
    atac_counts: Any,
    input_rows: Sequence[Mapping[str, str]],
    target_rows: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    import h5py

    with h5py.File(path, "x") as handle:
        handle.attrs["schema_version"] = "masld-bench-observed-multiome-model-input-h5-v1"
        handle.attrs["dataset_id"] = "gse296875"
        handle.attrs["stage"] = "smoke"
        handle.attrs["held_genomic_fold"] = held_fold
        handle.attrs["biological_unit"] = "donor_by_lineage"
        handle.attrs["target_atac_values_present"] = False
        rows = handle.create_group("rows")
        write_string_dataset(rows, "row_hash", row_hashes)
        write_string_dataset(rows, "unit_hash", [row["donor_hash"] for row in planned_units])
        write_string_dataset(rows, "lineage", [row["lineage"] for row in planned_units])
        rows.create_dataset(
            "donor_fold", data=np.asarray([int(row["donor_fold"]) for row in planned_units], dtype=np.int8)
        )
        rows.create_dataset(
            "nuclei", data=np.asarray([int(row["nuclei"]) for row in planned_units], dtype=np.int16)
        )
        rows.create_dataset("rna_observed_mask", data=np.ones(len(planned_units), dtype=np.uint8))
        rows.create_dataset("atac_observed_mask", data=np.ones(len(planned_units), dtype=np.uint8))
        rna = handle.create_group("rna")
        rna_identity = write_csr(rna.create_group("counts_csr"), rna_counts)
        write_string_dataset(rna, "ensembl_id", ensembl_ids)
        write_string_dataset(rna, "gene_name", gene_names)
        rna.attrs["scale"] = "raw_integer_donor_by_lineage_sum"
        atac = handle.create_group("observed_atac_input")
        atac_identity = write_csr(atac.create_group("counts_csr"), atac_counts)
        write_string_dataset(atac, "input_hash", [row["input_hash"] for row in input_rows])
        write_string_dataset(atac, "chromosome", [row["chromosome"] for row in input_rows])
        atac.create_dataset("bed_start_0based", data=np.asarray([int(row["bed_start_0based"]) for row in input_rows], dtype=np.int64))
        atac.create_dataset("bed_end_half_open", data=np.asarray([int(row["bed_end_half_open"]) for row in input_rows], dtype=np.int64))
        atac.attrs["scale"] = "raw_integer_donor_by_lineage_sum"
        targets = handle.create_group("targets_without_values")
        write_string_dataset(targets, "target_hash", [row["target_hash"] for row in target_rows])
        write_string_dataset(targets, "chromosome", [row["chromosome"] for row in target_rows])
        targets.create_dataset("bed_start_0based", data=np.asarray([int(row["bed_start_0based"]) for row in target_rows], dtype=np.int64))
        targets.create_dataset("bed_end_half_open", data=np.asarray([int(row["bed_end_half_open"]) for row in target_rows], dtype=np.int64))
        targets.attrs["atac_values_present"] = False
    return {"rna": rna_identity, "observed_atac_input": atac_identity}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    config = _load_json(config_path)
    resolved = validate_config(root, config)
    output = reject_symlink_components(args.output, label="materialization output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    mask_root = resolved["mask_plan"]
    input_rows_all = _read_tsv(
        mask_root / "observed_atac_input_peaks.tsv",
        (
            "input_hash",
            "chromosome",
            "bed_start_0based",
            "bed_end_half_open",
            "source_genomic_fold",
            "held_genomic_fold",
        ),
    )
    target_rows_all = _read_tsv(
        mask_root / "target_peaks.tsv",
        (
            "target_hash",
            "chromosome",
            "bed_start_0based",
            "bed_end_half_open",
            "genomic_fold",
        ),
    )
    namespaces = config["hash_namespaces"]
    expected = config["expected"]
    child_receipts: list[dict[str, Any]] = []
    import h5py

    with h5py.File(resolved["input_h5"], "r") as handle:
        if handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1":
            raise ModelInputMaterializationError("input HDF5 schema differs")
        cell_donors = decode_strings(handle["obs/donor_id"][:], label="donor")
        cell_lineages = decode_strings(handle["obs/broad_label"][:], label="lineage")
        planned_units, aggregation, row_hashes, raw_donor_ids = _rows_and_masks(
            mask_root, cell_donors, cell_lineages, namespaces
        )
        if len(cell_donors) != expected["nuclei"] or len(planned_units) != expected["donor_lineage_units"]:
            raise ModelInputMaterializationError("row census differs")
        rna_cell = read_full_csr(handle["rna/counts_csr"])
        if rna_cell.shape != (expected["nuclei"], expected["rna_features"]):
            raise ModelInputMaterializationError("RNA matrix shape differs")
        rna_counts = (aggregation @ rna_cell).tocsr()
        ensembl_ids = decode_strings(handle["rna/ensembl_id"][:], label="Ensembl ID")
        gene_names = decode_strings(handle["rna/gene_name"][:], label="gene name")
        peak_ids = decode_strings(handle["atac/peak_id"][:], label="peak ID")
        chromosomes = decode_strings(handle["atac/chromosome"][:], label="chromosome")
        starts = [int(value) for value in handle["atac/bed_start_0based"][:]]
        ends = [int(value) for value in handle["atac/bed_end_half_open"][:]]
        peak_index = _peak_index(peak_ids, chromosomes, starts, ends)
        atac_group = handle["atac/counts_csr"]
        atac_shape = tuple(int(value) for value in atac_group["shape"][:])
        if atac_shape != (expected["nuclei"], expected["atac_features"]):
            raise ModelInputMaterializationError("ATAC matrix shape differs")
        atac_indices = np.asarray(atac_group["indices"][:], dtype=np.int64)
        for held_fold in range(expected["genomic_folds"]):
            input_plan = [row for row in input_rows_all if int(row["held_genomic_fold"]) == held_fold]
            target_plan = [row for row in target_rows_all if int(row["genomic_fold"]) == held_fold]
            input_columns, input_rows = _resolve_mask_rows(
                peak_index=peak_index,
                rows=input_plan,
                held_fold=held_fold,
                peak_namespace=namespaces["peak"],
                kind="input",
            )
            target_columns, target_rows = _resolve_mask_rows(
                peak_index=peak_index,
                rows=target_plan,
                held_fold=held_fold,
                peak_namespace=namespaces["peak"],
                kind="target",
            )
            if len(input_columns) != expected["input_peaks_per_fold"] or len(target_columns) != expected["target_peaks_per_fold"]:
                raise ModelInputMaterializationError("fold peak budget differs")
            if set(input_columns).intersection(target_columns):
                raise ModelInputMaterializationError("held target column entered model input")
            input_positions_expected = selected_value_positions(
                atac_indices, input_columns, total_columns=expected["atac_features"]
            )
            target_positions = selected_value_positions(
                atac_indices, target_columns, total_columns=expected["atac_features"]
            )
            if np.intersect1d(input_positions_expected, target_positions).size:
                raise ModelInputMaterializationError("held target value position entered model input")
            atac_cell, input_positions_read = selectively_read_csr_columns(
                atac_group, input_columns
            )
            if not np.array_equal(input_positions_read, input_positions_expected):
                raise ModelInputMaterializationError("selective ATAC value positions differ")
            atac_counts = (aggregation @ atac_cell).tocsr()
            child = output / f"genomic_fold_{held_fold}"
            child.mkdir()
            identities = _write_fold_h5(
                path=child / "model_input.h5",
                held_fold=held_fold,
                planned_units=planned_units,
                row_hashes=row_hashes,
                rna_counts=rna_counts,
                ensembl_ids=ensembl_ids,
                gene_names=gene_names,
                atac_counts=atac_counts,
                input_rows=input_rows,
                target_rows=target_rows,
            )
            child_receipt = {
                "schema_version": "masld-bench-observed-multiome-model-input-receipt-v1",
                "dataset_id": "gse296875",
                "stage": "smoke",
                "held_genomic_fold": held_fold,
                "donor_lineage_units": len(planned_units),
                "rna_features": len(ensembl_ids),
                "observed_atac_input_peaks": len(input_columns),
                "target_identifiers_without_values": len(target_columns),
                "rna_identity": identities["rna"],
                "observed_atac_input_identity": identities["observed_atac_input"],
                "atac_sparse_value_positions_read": int(input_positions_read.size),
                "held_target_sparse_value_positions": int(target_positions.size),
                "same_fold_input_target_value_position_overlap": 0,
                "target_atac_values_read": False,
                "target_atac_values_written": False,
                "full_atac_data_array_read": False,
                "normalization_applied": False,
                "feature_selection_applied": False,
                "raw_donor_ids_exported": False,
                "model_fit": False,
                "metric_calculated": False,
                "sealed_outcomes_read": False,
            }
            write_json_exclusive(child / "receipt.json", child_receipt)
            child_hash = freeze_tree(
                child,
                metadata={
                    "artifact_class": "gse296875_observed_multiome_model_input",
                    "held_genomic_fold": held_fold,
                    "target_atac_values_present": False,
                    "sealed_outcomes_accessed": False,
                },
            )
            child_receipts.append(
                {
                    "held_genomic_fold": held_fold,
                    "path": child.relative_to(output).as_posix(),
                    "artifacts_sha256": child_hash,
                }
            )
    exported_tokens = {
        value
        for child in output.glob("genomic_fold_*")
        for path in child.glob("*.json")
        for line in path.read_text(encoding="utf-8").splitlines()
        for value in line.replace('"', " ").replace(":", " ").replace(",", " ").split()
    }
    if set(raw_donor_ids).intersection(exported_tokens):
        raise ModelInputMaterializationError("raw donor ID entered JSON output")
    parent_receipt = {
        "schema_version": "masld-bench-observed-multiome-model-input-campaign-receipt-v1",
        "campaign_id": CAMPAIGN_ID,
        "dataset_id": "gse296875",
        "stage": "smoke",
        "child_artifacts": child_receipts,
        "child_artifact_count": len(child_receipts),
        "donors": 39,
        "donor_lineage_units": 195,
        "target_atac_values_read": False,
        "target_atac_values_written": False,
        "evaluator_artifact_read": False,
        "full_atac_data_array_read": False,
        "development_metric_calculated": False,
        "sealed_outcomes_read": False,
        "model_fit": False,
        "ranking_authorized": False,
        "champion_claim_allowed": False,
        "next_gate": "separate_evaluator_only_target_count_materialization",
    }
    write_json_exclusive(output / "receipt.json", parent_receipt)
    sources = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "src/masld_bench/observed_multiome_materialization.py",
        root / "tests/unit/test_observed_multiome_model_input_materialization.py",
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8"
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "gse296875_observed_multiome_model_input_campaign",
            "campaign_id": CAMPAIGN_ID,
            "child_artifact_count": 5,
            "target_atac_values_present": False,
            "sealed_outcomes_accessed": False,
            "model_fit": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
