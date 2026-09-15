#!/usr/bin/env python3
"""Build a condition-blind GSE244832 fixed-window ATAC exchange axis."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from scripts.build_gse281367_label_free_atac_exchange import (
    CANONICAL_LINEAGES,
    MEMBERSHIP_FIELDS,
    PRIMARY_LINEAGES,
    WINDOW_FIELDS,
    LabelFreeATACExchangeError,
    _digest,
    _hash_source,
    _load_toml,
    _read_inventory,
    _verify_tree,
    _write_membership,
    _write_tsv,
    fold_index,
    read_label_free_obs,
    read_windows,
    resolve_raw_barcodes,
)


SCHEMA = "masld-bench-gse244832-label-free-atac-exchange-v1"
EXCHANGE_ID = "gse296875_to_gse244832_fixed_window_atac_transport_20260825"
DATASET_ID = "gse244832"
DONORS = tuple(f"D{index:02d}" for index in range(1, 19))
LABEL_CONTRACT = {
    "Cholangiocytes": ("cholangiocyte", "primary"),
    "Circulating_NK_NKT": ("t_nk_cell", "secondary"),
    "Endothelial_cells": ("endothelial_cell", "secondary"),
    "Fibroblasts": ("fibroblast", "primary"),
    "Hepatocytes": ("hepatocyte", "primary"),
    "Low_confidence": ("low_confidence", "excluded"),
    "Macrophages": ("macrophage", "primary"),
    "Plasma_cells": ("b_cell", "secondary"),
    "Resident_NK": ("t_nk_cell", "secondary"),
    "T_cells": ("t_nk_cell", "secondary"),
}


def _resolve_authority(
    root: Path,
    contract: Mapping[str, Any],
    section: str,
    artifact_class: str,
) -> Path:
    spec = contract.get(section)
    if not isinstance(spec, dict):
        raise LabelFreeATACExchangeError(f"missing authority: {section}")
    path = (root / str(spec.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    _verify_tree(path, str(spec.get("artifacts_sha256", "")), artifact_class)
    return path


def validate_contract(root: Path, contract: Mapping[str, Any]) -> dict[str, Any]:
    if contract.get("schema_version") != SCHEMA or contract.get("exchange_id") != EXCHANGE_ID:
        raise LabelFreeATACExchangeError("exchange identity differs")
    if contract.get("source_dataset_id") != "gse296875" or contract.get("target_dataset_id") != DATASET_ID:
        raise LabelFreeATACExchangeError("source or target dataset differs")
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
        "rna_assay_available_to_atac_model_jobs",
        "cross_assay_participant_join_inferred",
    ):
        if contract.get(field) is not False:
            raise LabelFreeATACExchangeError(f"label, pairing, or claim firewall opened: {field}")

    activation_root = _resolve_authority(
        root, contract, "activation_authority", "gse244832_activation_readiness_audit"
    )
    source_root = _resolve_authority(
        root, contract, "source_authority", "atac_transport_source_admission_audit"
    )
    window_root = _resolve_authority(
        root, contract, "window_authority", "atac_transport_evaluator_outcomes"
    )
    census_root = _resolve_authority(
        root, contract, "census_authority", "observed_multiome_additive_census_revision"
    )

    source = contract["source_authority"]
    window = contract["window_authority"]
    source_inventory = source_root / str(source["inventory_path"])
    donor_inventory_root = (root / str(source["per_donor_inventory_authority"])).resolve(strict=True)
    donor_inventory_root.relative_to(root)
    _verify_tree(
        donor_inventory_root,
        str(source["per_donor_inventory_artifacts_sha256"]),
        "atac_transport_cell_membership",
    )
    donor_inventory = donor_inventory_root / str(source["per_donor_inventory_path"])
    windows = window_root / str(window["coordinate_only_path"])
    for path, expected in (
        (source_inventory, source["inventory_sha256"]),
        (donor_inventory, source["per_donor_inventory_sha256"]),
        (windows, window["coordinate_only_sha256"]),
    ):
        if _digest(path) != expected:
            raise LabelFreeATACExchangeError(f"bound file differs: {path}")

    family_spec = contract.get("family_eligibility_authority", {})
    family_path = (root / str(family_spec.get("path", ""))).resolve(strict=True)
    family_path.relative_to(root)
    if family_path.is_symlink() or _digest(family_path) != family_spec.get("sha256"):
        raise LabelFreeATACExchangeError("family-eligibility authority differs")
    with family_path.open("rb") as handle:
        family_contract = tomllib.load(handle)
    family = family_contract.get("family_eligibility")
    if not isinstance(family, list) or len(family) != 25:
        raise LabelFreeATACExchangeError("family-native eligibility roster differs")
    by_id = {row.get("model_id"): row for row in family}
    if len(by_id) != 25 or by_id.get("epibert", {}).get("query_atac_allowed") is not True:
        raise LabelFreeATACExchangeError("observed-ATAC family eligibility differs")
    if by_id.get("scbasset", {}).get("query_atac_allowed") is not False:
        raise LabelFreeATACExchangeError("sequence-only family eligibility differs")

    axis = contract.get("label_free_axis", {})
    exchange = contract.get("window_exchange", {})
    rotations = contract.get("rotation", [])
    if (
        axis.get("donor_count") != 18
        or tuple(axis.get("primary_lineages", ())) != PRIMARY_LINEAGES
        or tuple(axis.get("canonical_lineages", ())) != CANONICAL_LINEAGES
        or axis.get("minimum_cells") != 50
        or axis.get("outer_fold_uses_condition") is not False
        or axis.get("missing_as_zero") is not False
        or axis.get("rna_state") != "structurally_missing_for_atac_exchange"
        or exchange.get("window_count") != 32000
        or exchange.get("query_fragment_cut_sites") != "start_plus_4_and_end_minus_5"
        or exchange.get("target_role_contigs_removed_from_query_atac") is not True
        or [(row.get("observed_atac_context_role"), row.get("scored_target_role")) for row in rotations]
        != [("valid", "test"), ("test", "valid")]
    ):
        raise LabelFreeATACExchangeError("axis, fragment geometry, mask, or rotation contract differs")
    if contract.get("claim_boundary", {}).get("gse244832_is_non_champion") is not True:
        raise LabelFreeATACExchangeError("GSE244832 non-champion boundary differs")
    if contract.get("claim_boundary", {}).get("cross_assay_join_unresolved") is not True:
        raise LabelFreeATACExchangeError("GSE244832 cross-assay boundary differs")

    implementation = contract.get("implementation", {})
    for prefix in ("builder", "builder_test"):
        path = (root / str(implementation.get(f"{prefix}_path", ""))).resolve(strict=True)
        path.relative_to(root)
        if path.is_symlink() or not path.is_file() or _digest(path) != implementation.get(f"{prefix}_sha256"):
            raise LabelFreeATACExchangeError(f"implementation authority differs: {prefix}")
    decision = contract.get("decision", {})
    if (
        decision.get("label_free_axis_build_authorized") is not True
        or decision.get("biological_scoring_authorized_before_prediction_commit") is not False
        or decision.get("cross_assay_fusion_authorized") is not False
        or decision.get("outcome_selected_model_repair_allowed") is not False
    ):
        raise LabelFreeATACExchangeError("exchange decision differs")
    return {
        "activation_root": activation_root,
        "source_inventory": source_inventory,
        "donor_inventory": donor_inventory,
        "windows": windows,
        "census_root": census_root,
        "family": family,
    }


def build_membership(
    merged_path: Path,
    donor_records: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Read identifiers and annotations only; never open the condition column."""

    merged = read_label_free_obs(merged_path, ("donor_id", "cell_type"))
    if sorted(set(merged["donor_id"])) != list(DONORS):
        raise LabelFreeATACExchangeError("target donor axis differs")
    by_donor = {Path(record["path"]).stem: record for record in donor_records}
    if set(by_donor) != set(DONORS):
        raise LabelFreeATACExchangeError("per-donor H5AD roster differs")
    rows: list[dict[str, Any]] = []
    counts: Counter[tuple[str, str]] = Counter()
    for donor in DONORS:
        indexes = [index for index, value in enumerate(merged["donor_id"]) if value == donor]
        raw = read_label_free_obs(Path(by_donor[donor]["path"]), ())
        raw_barcodes = resolve_raw_barcodes(
            [merged["source_cell_id"][index] for index in indexes],
            raw["source_cell_id"],
        )
        for index, barcode in zip(indexes, raw_barcodes, strict=True):
            source_label = merged["cell_type"][index]
            if source_label not in LABEL_CONTRACT:
                raise LabelFreeATACExchangeError(f"unregistered cell label: {source_label}")
            lineage, role = LABEL_CONTRACT[source_label]
            rows.append(
                {
                    "dataset_id": DATASET_ID,
                    "donor_id": donor,
                    "raw_barcode": barcode,
                    "lineage_id": lineage,
                    "analysis_role": role,
                    "outer_fold": fold_index(DATASET_ID, donor),
                }
            )
            if role == "primary":
                counts[(donor, lineage)] += 1
    if len(rows) != 88814 or len({(row["donor_id"], row["raw_barcode"]) for row in rows}) != len(rows):
        raise LabelFreeATACExchangeError("project-frozen cell census or barcode uniqueness differs")

    axis_rows: list[dict[str, Any]] = []
    for donor_index, donor in enumerate(DONORS):
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
                    "dataset_id": DATASET_ID,
                    "donor_index": donor_index,
                    "donor_id": donor,
                    "lineage_index": lineage_index,
                    "lineage_id": lineage,
                    "outer_fold": fold_index(DATASET_ID, donor),
                    "cells": "" if cells is None else cells,
                    "evidence_state": state,
                    "eligible_min_50_cells": str(eligible).lower(),
                }
            )
    return rows, axis_rows


def build(root: Path, contract_path: Path, output: Path, hash_workers: int) -> dict[str, Any]:
    if output.exists() or not 1 <= hash_workers <= 4:
        raise LabelFreeATACExchangeError("output or hash-worker contract differs")
    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    contract_path.relative_to(root)
    contract = _load_toml(contract_path)
    resolved = validate_contract(root, contract)
    source_rows = _read_inventory(resolved["source_inventory"], ("path", "size_bytes", "sha256"))
    merged_records = [
        row for row in source_rows
        if row["path"].endswith(contract["source_authority"]["merged_h5ad_suffix"])
    ]
    donor_rows = [
        row
        for row in _read_inventory(
            resolved["donor_inventory"],
            ("dataset_id", "path", "size_bytes", "sha256"),
        )
        if row["dataset_id"] == DATASET_ID
    ]
    if len(merged_records) != 1 or len(donor_rows) != 18:
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
        resolved["family"],
    )
    write_json_exclusive(output / "missingness_contract.json", contract["missingness_contract"])
    write_json_exclusive(output / "overlap_contract.json", contract["overlap_contract"])
    summary = {
        "schema_version": "masld-bench-gse244832-label-free-atac-exchange-summary-v1",
        "exchange_id": EXCHANGE_ID,
        "source_dataset_id": "gse296875",
        "target_dataset_id": DATASET_ID,
        "donors": 18,
        "membership_cells": len(membership),
        "donor_lineage_rows": len(axis),
        "eligible_donor_lineage_rows": sum(row["eligible_min_50_cells"] == "true" for row in axis),
        "windows": len(windows),
        "rotations": 2,
        "condition_h5_dataset_opened": False,
        "condition_values_read": False,
        "condition_values_used": False,
        "rna_assay_opened": False,
        "cross_assay_join_inferred": False,
        "biological_outcome_arrays_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
        "gse244832_champion_claim_allowed": False,
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
            "artifact_class": "gse244832_label_free_atac_exchange_axis",
            "exchange_id": EXCHANGE_ID,
            "source_dataset_id": "gse296875",
            "target_dataset_id": DATASET_ID,
            "condition_values_read": False,
            "rna_assay_opened": False,
            "cross_assay_join_inferred": False,
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
