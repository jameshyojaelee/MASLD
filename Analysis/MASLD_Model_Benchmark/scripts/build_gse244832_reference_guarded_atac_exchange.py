#!/usr/bin/env python3
"""Rebuild the GSE244832 ATAC exchange with target-reference family guards."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from scripts.build_gse244832_label_free_atac_exchange import (
    DATASET_ID,
    DONORS,
    EXCHANGE_ID,
    build_membership,
    validate_contract,
)
from scripts.build_gse281367_label_free_atac_exchange import (
    WINDOW_FIELDS,
    LabelFreeATACExchangeError,
    _digest,
    _hash_source,
    _load_toml,
    _read_inventory,
    _write_membership,
    _write_tsv,
    read_windows,
)


GUARD_SCHEMA = "masld-bench-gse244832-atac-reference-guard-v1"
REVISION_ID = "gse244832_label_free_atac_exchange_reference_guard_20260825"


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise LabelFreeATACExchangeError("reference guard must be an object")
    return value


def target_family_eligibility(
    family: list[dict[str, Any]], guard: Mapping[str, Any]
) -> list[dict[str, Any]]:
    blocked = set(guard.get("blocked_input_regimes", ()))
    allowed = set(guard.get("allowed_input_regimes", ()))
    if allowed != {"observed_atac"} or "sequence_only" not in blocked:
        raise LabelFreeATACExchangeError("reference-guard input regimes differ")
    result: list[dict[str, Any]] = []
    for source in family:
        row = dict(source)
        regime = str(row.get("input_regime"))
        if regime in blocked:
            row["eligibility"] = guard["blocked_eligibility"]
            row["query_atac_allowed"] = False
        elif regime not in allowed:
            raise LabelFreeATACExchangeError(f"unclassified target input regime: {regime}")
        result.append(row)
    if len(result) != 25 or len({row["model_id"] for row in result}) != 25:
        raise LabelFreeATACExchangeError("target family roster differs")
    if any(
        row["input_regime"] != "observed_atac" and not row["eligibility"].startswith("ineligible_")
        for row in result
    ):
        raise LabelFreeATACExchangeError("reference- or RNA-dependent family remained executable")
    if not any(row["input_regime"] == "observed_atac" for row in result):
        raise LabelFreeATACExchangeError("observed-ATAC family roster is empty")
    return result


def validate_guard(
    root: Path,
    contract_path: Path,
    guard_path: Path,
    implementation_path: Path,
    implementation_test_path: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    contract = _load_toml(contract_path)
    resolved = validate_contract(root, contract)
    guard = _load_json(guard_path)
    if (
        guard.get("schema_version") != GUARD_SCHEMA
        or guard.get("revision_id") != REVISION_ID
        or guard.get("dataset_id") != DATASET_ID
        or (root / guard.get("base_contract_path", "")).resolve() != contract_path
        or guard.get("base_contract_sha256") != _digest(contract_path)
        or guard.get("source_native_reference_bundle_resolved") is not False
        or guard.get("source_native_exact_annotation_resolved") is not False
        or guard.get("sequence_or_coordinate_family_execution_authorized") is not False
        or guard.get("observed_atac_query_materialization_authorized") is not True
        or guard.get("rna_or_cross_assay_input_authorized") is not False
        or guard.get("cross_assay_participant_join_inferred") is not False
        or guard.get("condition_values_read") is not False
        or guard.get("rna_assay_opened") is not False
        or guard.get("biological_outcome_arrays_read") is not False
        or guard.get("prediction_or_metric_files_produced") is not False
        or guard.get("champion_claim_allowed") is not False
    ):
        raise LabelFreeATACExchangeError("reference guard differs")
    activation = (root / str(guard.get("activation_authority_path", ""))).resolve(strict=True)
    if _digest(activation / "ARTIFACTS.json") != guard.get("activation_authority_artifacts_sha256"):
        raise LabelFreeATACExchangeError("reference-guard activation authority differs")
    implementation = guard.get("implementation", {})
    if (
        _digest(implementation_path) != implementation.get("builder_sha256")
        or _digest(implementation_test_path) != implementation.get("builder_test_sha256")
    ):
        raise LabelFreeATACExchangeError("reference-guard implementation differs")
    resolved["family"] = target_family_eligibility(resolved["family"], guard)
    return contract, guard, resolved


def build(
    root: Path,
    contract_path: Path,
    guard_path: Path,
    output: Path,
    hash_workers: int,
) -> dict[str, Any]:
    if output.exists() or not 1 <= hash_workers <= 4:
        raise LabelFreeATACExchangeError("output or hash-worker contract differs")
    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    guard_path = guard_path.resolve(strict=True)
    implementation_path = Path(__file__).resolve(strict=True)
    implementation_test_path = (
        root / "tests/unit/test_build_gse244832_reference_guarded_atac_exchange.py"
    ).resolve(strict=True)
    for path in (contract_path, guard_path, implementation_path, implementation_test_path):
        path.relative_to(root)
    contract, guard, resolved = validate_guard(
        root, contract_path, guard_path, implementation_path, implementation_test_path
    )
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
    if len(merged_records) != 1 or len(donor_rows) != len(DONORS):
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
    write_json_exclusive(output / "reference_guard.json", guard)
    write_json_exclusive(output / "missingness_contract.json", contract["missingness_contract"])
    write_json_exclusive(output / "overlap_contract.json", contract["overlap_contract"])
    summary = {
        "schema_version": "masld-bench-gse244832-reference-guarded-atac-exchange-summary-v1",
        "revision_id": REVISION_ID,
        "exchange_id": EXCHANGE_ID,
        "source_dataset_id": "gse296875",
        "target_dataset_id": DATASET_ID,
        "donors": len(DONORS),
        "membership_cells": len(membership),
        "donor_lineage_rows": len(axis),
        "eligible_donor_lineage_rows": sum(row["eligible_min_50_cells"] == "true" for row in axis),
        "family_rows": len(resolved["family"]),
        "observed_atac_family_rows": sum(
            row["input_regime"] == "observed_atac" for row in resolved["family"]
        ),
        "non_observed_atac_executable_rows": sum(
            row["input_regime"] != "observed_atac"
            and not row["eligibility"].startswith("ineligible_")
            for row in resolved["family"]
        ),
        "windows": len(windows),
        "rotations": 2,
        "source_native_reference_bundle_resolved": False,
        "sequence_or_coordinate_family_execution_authorized": False,
        "observed_atac_query_materialization_authorized": True,
        "condition_values_read": False,
        "rna_assay_opened": False,
        "cross_assay_join_inferred": False,
        "biological_outcome_arrays_read": False,
        "metrics_calculated": False,
        "gse244832_champion_claim_allowed": False,
        "source_hashes_reverified": len(verified),
    }
    write_json_exclusive(output / "exchange_summary.json", summary)
    (output / "source.sha256").write_text(
        "".join(
            f"{_digest(path)}  {path}\n"
            for path in (contract_path, guard_path, implementation_path, implementation_test_path)
        ),
        encoding="utf-8",
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse244832_reference_guarded_atac_exchange_axis",
            "revision_id": REVISION_ID,
            "source_dataset_id": "gse296875",
            "target_dataset_id": DATASET_ID,
            "condition_values_read": False,
            "rna_assay_opened": False,
            "cross_assay_join_inferred": False,
            "biological_outcome_arrays_read": False,
            "sequence_execution_authorized": False,
            "observed_atac_query_materialization_authorized": True,
            "champion_claim_allowed": False,
            "status": "passed",
        },
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--guard", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hash-workers", type=int, default=2)
    arguments = parser.parse_args()
    print(
        json.dumps(
            build(
                arguments.root,
                arguments.contract,
                arguments.guard,
                arguments.output,
                arguments.hash_workers,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
