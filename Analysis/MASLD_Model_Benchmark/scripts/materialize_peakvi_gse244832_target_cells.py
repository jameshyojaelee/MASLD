#!/usr/bin/env python3
"""Materialize outcome-blind GSE244832 cell contexts for source-fit PeakVI."""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import shutil
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.materialize_peakvi_cross_cohort_cells import (
    LINEAGES,
    ROLES,
    PeakVICellMaterializationError,
    WindowIndex,
    digest,
    merge_shards,
    read_json,
    read_reference_lengths,
    read_windows,
    windows_at,
    write_axis,
    write_windows,
)
from scripts.materialize_gse244832_label_free_query_atac import tn5_positions


SCHEMA = "masld-bench-peakvi-gse244832-target-cell-materialization-v1"
DATASET_ID = "gse244832"
DONORS = tuple(f"D{index:02d}" for index in range(1, 19))
MASK_COUNTS = {0: 48, 1: 4, 3: 20}


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def resolve_inside(root: Path, value: str, label: str) -> Path:
    path = (root / value).resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as error:
        raise PeakVICellMaterializationError(f"{label} escapes benchmark root") from error
    if path.is_symlink():
        raise PeakVICellMaterializationError(f"{label} is a symlink")
    return path


def verify_tree(root: Path, record: Mapping[str, Any], artifact_class: str) -> Path:
    path = resolve_inside(root, str(record.get("path", "")), artifact_class)
    manifest = verify_frozen_tree(path)
    if (
        digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise PeakVICellMaterializationError(f"{artifact_class} authority differs")
    return path


def verify_manifest(root: Path, record: Mapping[str, Any], artifact_class: str) -> Path:
    path = resolve_inside(root, str(record.get("path", "")), artifact_class)
    manifest = read_json(path / "ARTIFACTS.json")
    if (
        digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise PeakVICellMaterializationError(f"{artifact_class} manifest differs")
    return path


def read_unit_axis(path: Path) -> tuple[list[dict[str, str]], Any]:
    import numpy as np

    fields = (
        "donor_index", "donor_id", "lineage_index", "lineage_id", "outer_fold",
        "cells", "evidence_state", "missing_state_code", "eligible_min_50_cells",
    )
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise PeakVICellMaterializationError("registered donor-lineage fields differ")
        rows = list(reader)
    if len(rows) != 72:
        raise PeakVICellMaterializationError("registered donor-lineage census differs")
    states = np.empty((18, 4), dtype=np.uint8)
    expected_state = {0: "observed", 1: "structurally_missing", 3: "below_qc"}
    for offset, row in enumerate(rows):
        donor_index, lineage_index = divmod(offset, 4)
        state = int(row["missing_state_code"])
        if (
            int(row["donor_index"]) != donor_index
            or row["donor_id"] != DONORS[donor_index]
            or int(row["lineage_index"]) != lineage_index
            or row["lineage_id"] != LINEAGES[lineage_index]
            or state not in expected_state
            or row["evidence_state"] != expected_state[state]
            or (row["eligible_min_50_cells"] == "true") != (state == 0)
            or (int(row["cells"]) >= 50) != (state == 0)
        ):
            raise PeakVICellMaterializationError("registered donor-lineage row differs")
        states[donor_index, lineage_index] = state
    if {state: int((states == state).sum()) for state in MASK_COUNTS} != MASK_COUNTS:
        raise PeakVICellMaterializationError("registered missingness counts differ")
    return rows, states


def read_membership(path: Path) -> list[dict[str, str]]:
    fields = (
        "dataset_id", "donor_id", "raw_barcode", "lineage_id", "analysis_role", "outer_fold"
    )
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise PeakVICellMaterializationError("GSE244832 membership fields differ")
        for row in reader:
            identity = (row["donor_id"], row["raw_barcode"])
            if (
                row["dataset_id"] != DATASET_ID
                or row["donor_id"] not in DONORS
                or identity in seen
            ):
                raise PeakVICellMaterializationError("GSE244832 membership identity differs")
            seen.add(identity)
            rows.append(row)
    if {row["donor_id"] for row in rows} != set(DONORS):
        raise PeakVICellMaterializationError("GSE244832 membership donor roster differs")
    return rows


def select_cells(
    membership: Sequence[Mapping[str, str]],
    unit_rows: Sequence[Mapping[str, str]],
    unit_states: Any,
    *,
    cap: int,
    seed: int,
) -> list[dict[str, Any]]:
    """Apply the registered unit mask before deterministic ATAC-cell selection."""

    if cap != 64 or seed != 20260825:
        raise PeakVICellMaterializationError("cell-selection contract differs")
    grouped: dict[tuple[str, str], list[Mapping[str, str]]] = defaultdict(list)
    for row in membership:
        if row["analysis_role"] == "primary" and row["lineage_id"] in LINEAGES:
            grouped[(row["donor_id"], row["lineage_id"])].append(row)
    selected: list[dict[str, Any]] = []
    for unit_offset, unit in enumerate(unit_rows):
        donor_index, lineage_index = divmod(unit_offset, 4)
        donor, lineage = unit["donor_id"], unit["lineage_id"]
        candidates = grouped.get((donor, lineage), [])
        if len(candidates) != int(unit["cells"]):
            raise PeakVICellMaterializationError("membership and registered unit count differ")
        if int(unit_states[donor_index, lineage_index]) != 0:
            continue
        ordered = sorted(
            candidates,
            key=lambda row: sha256(
                f"{seed}\0{DATASET_ID}\0{donor}\0{row['raw_barcode']}".encode()
            ).hexdigest(),
        )
        for rank, row in enumerate(ordered[:cap]):
            barcode = row["raw_barcode"]
            selected.append(
                {
                    "dataset_id": DATASET_ID,
                    "cell_id": barcode,
                    "cell_hash": sha256(
                        f"peakvi-cross-cohort-v2\0{DATASET_ID}\0{donor}\0{barcode}".encode()
                    ).hexdigest(),
                    "raw_barcode": barcode,
                    "well_id": donor,
                    "donor_id": donor,
                    "lineage_id": lineage,
                    "outer_fold": row["outer_fold"],
                    "selection_rank": rank,
                    "available_nuclei": len(ordered),
                }
            )
    observed_units = {
        (row["donor_id"], row["lineage_id"])
        for row in selected
    }
    expected_units = {
        (row["donor_id"], row["lineage_id"])
        for offset, row in enumerate(unit_rows)
        if int(unit_states.reshape(-1)[offset]) == 0
    }
    if (
        observed_units != expected_units
        or len(observed_units) != 48
        or len({row["cell_hash"] for row in selected}) != len(selected)
        or any(int(row["available_nuclei"]) < 50 for row in selected)
    ):
        raise PeakVICellMaterializationError("eligible-cell selection differs")
    return selected


def read_inventory(path: Path, fragment_root: Path) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("path", "size_bytes", "sha256"):
            raise PeakVICellMaterializationError("GSE244832 source inventory fields differ")
        for row in reader:
            candidate = Path(row["path"])
            if candidate.parent != fragment_root:
                continue
            match = re.fullmatch(r"(D[0-9]{2})_fragments[.]tsv[.]gz", candidate.name)
            if match is not None:
                output[match.group(1)] = row
    if tuple(sorted(output)) != DONORS:
        raise PeakVICellMaterializationError("GSE244832 fragment roster differs")
    return output


def process_fragment_file(
    *,
    donor_id: str,
    path: str,
    expected_size: int,
    expected_sha256: str,
    selected_lookup: Mapping[str, int],
    windows: Mapping[str, WindowIndex],
    reference_lengths: Mapping[str, int],
    output: str,
    n_rows: int,
) -> dict[str, Any]:
    """Write one target shard using custom-Bowtie +4/-5 cut geometry."""

    import numpy as np
    from scipy import sparse

    source = Path(path).resolve(strict=True)
    if source.is_symlink() or source.stat().st_size != expected_size:
        raise PeakVICellMaterializationError(f"fragment path or size differs: {donor_id}")
    role_rows: dict[str, list[int]] = {role: [] for role in ROLES}
    role_columns: dict[str, list[int]] = {role: [] for role in ROLES}
    seen: set[str] = set()
    prior_contig: str | None = None
    prior_start = -1
    closed: set[str] = set()
    selected_fragment_rows = 0
    compressed_digest = sha256()

    class DigestReader(io.RawIOBase):
        def __init__(self, handle: Any) -> None:
            self.handle = handle

        def readable(self) -> bool:
            return True

        def readinto(self, buffer: Any) -> int:
            data = self.handle.read(len(buffer))
            if not data:
                return 0
            compressed_digest.update(data)
            buffer[: len(data)] = data
            return len(data)

    with source.open("rb") as raw:
        with gzip.GzipFile(fileobj=DigestReader(raw), mode="rb") as compressed:
            handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
            for line_number, line in enumerate(handle, start=1):
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise PeakVICellMaterializationError(
                        f"fragment width differs: {donor_id}:{line_number}"
                    )
                contig, start_text, end_text, barcode = fields[:4]
                try:
                    start, end = int(start_text), int(end_text)
                except ValueError as error:
                    raise PeakVICellMaterializationError("fragment coordinate differs") from error
                if (
                    contig not in reference_lengths
                    or start < 0
                    or end <= start
                    or end > reference_lengths[contig]
                ):
                    raise PeakVICellMaterializationError("fragment geometry differs")
                if contig != prior_contig:
                    if contig in closed:
                        raise PeakVICellMaterializationError("fragment contig reappears")
                    if prior_contig is not None:
                        closed.add(prior_contig)
                    prior_contig, prior_start = contig, -1
                if start < prior_start:
                    raise PeakVICellMaterializationError("fragments are unsorted")
                prior_start = start
                row_index = selected_lookup.get(barcode)
                if row_index is None:
                    continue
                seen.add(barcode)
                selected_fragment_rows += 1
                for position in tn5_positions(start, end):
                    for role in ROLES:
                        for column in windows_at(windows[role], contig, position):
                            role_rows[role].append(row_index)
                            role_columns[role].append(column)
    if compressed_digest.hexdigest() != expected_sha256:
        raise PeakVICellMaterializationError(f"fragment SHA-256 differs: {donor_id}")
    if len(seen) != len(selected_lookup):
        raise PeakVICellMaterializationError(
            f"selected barcodes absent from fragment file: {donor_id}"
        )
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=False)
    identities: dict[str, Any] = {}
    for role in ROLES:
        matrix = sparse.coo_matrix(
            (
                np.ones(len(role_rows[role]), dtype=np.uint16),
                (
                    np.asarray(role_rows[role], dtype=np.int32),
                    np.asarray(role_columns[role], dtype=np.int32),
                ),
            ),
            shape=(n_rows, 16_000),
            dtype=np.uint16,
        ).tocsr()
        matrix.sum_duplicates()
        matrix.sort_indices()
        sparse.save_npz(output_path / f"{role}.counts.npz", matrix, compressed=True)
        identities[role] = {
            "shape": list(matrix.shape),
            "nnz": int(matrix.nnz),
            "sum": int(matrix.sum(dtype=np.uint64)),
            "sha256": digest(output_path / f"{role}.counts.npz"),
        }
    return {
        "donor_id": donor_id,
        "path": str(source),
        "selected_barcodes": len(selected_lookup),
        "selected_barcodes_seen": len(seen),
        "selected_fragment_rows": selected_fragment_rows,
        "gzip_crc_verified_to_eof": True,
        "compressed_sha256_verified": True,
        "fragment_cut_sites": "start_plus_4_and_end_minus_5",
        "matrices": identities,
    }


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_target_output(
    output: Path,
    role: str,
    cells: Sequence[Mapping[str, Any]],
    windows: WindowIndex,
    matrix: Any,
    aggregate_query: Path,
    unit_rows: Sequence[Mapping[str, str]],
    unit_states: Any,
    mask_authority: Path,
    authorities: Mapping[str, Any],
) -> None:
    import numpy as np
    from scipy import sparse

    output.mkdir(parents=True, mode=0o750)
    write_axis(output / "cell_axis.tsv", cells)
    sparse.save_npz(output / "context.counts.npz", matrix, compressed=True)
    write_windows(output / "context.windows.tsv", windows)
    write_tsv(output / "donor_lineage_axis.tsv", tuple(unit_rows[0]), unit_rows)
    shutil.copyfile(mask_authority, output / "donor_lineage_missing_state.uint8.npy")
    if digest(output / "donor_lineage_missing_state.uint8.npy") != digest(mask_authority):
        raise PeakVICellMaterializationError("registered missingness copy differs")
    aggregate = np.load(
        aggregate_query / "query_unit_fragment_total.uint32.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    if aggregate.shape != (18, 4, 16_000) or aggregate.dtype != np.uint32:
        raise PeakVICellMaterializationError("aggregate query geometry differs")
    donor_index = {donor: index for index, donor in enumerate(DONORS)}
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    selected_total = np.zeros((18, 4, 16_000), dtype=np.uint64)
    selected_groups = np.zeros((18, 4), dtype=np.int64)
    for row_index, row in enumerate(cells):
        group = (donor_index[row["donor_id"]], lineage_index[row["lineage_id"]])
        selected_total[group] += matrix.getrow(row_index).toarray()[0]
        selected_groups[group] += 1
    if (
        np.any(selected_total > np.asarray(aggregate, dtype=np.uint64))
        or np.any(selected_groups[unit_states != 0] != 0)
        or np.any(selected_groups[unit_states == 0] == 0)
    ):
        raise PeakVICellMaterializationError("selected cells violate aggregate or mask boundary")
    write_json_exclusive(output / "input_authorities.json", authorities)
    write_json_exclusive(
        output / "contract.json",
        {
            "schema_version": "masld-bench-peakvi-gse244832-target-cell-context-v1",
            "dataset_id": DATASET_ID,
            "context_role": role,
            "scored_target_role": "test" if role == "valid" else "valid",
            "cells": len(cells),
            "donors": 18,
            "lineages": list(LINEAGES),
            "eligible_donor_lineage_units": 48,
            "structurally_missing_units": 4,
            "below_qc_units": 20,
            "features": 16_000,
            "matrix_unit": "deduplicated_custom_Bowtie2_Tn5_cut_count_per_selected_cell",
            "fragment_cut_sites": "start_plus_4_and_end_minus_5",
            "native_peakvi_transform": "binary_x_greater_than_zero_applied_only_in_model_adapter",
            "registered_mask_applied_before_cell_selection": True,
            "ineligible_cells_entered_transform": False,
            "selected_cell_counts_bounded_by_frozen_aggregate_query": True,
            "cell_selection_assay": "atac_only",
            "other_rotation_values_present": False,
            "scored_target_contigs_present": False,
            "condition_values_read": False,
            "phenotype_values_read": False,
            "rna_assay_read": False,
            "cross_assay_join_consumed": False,
            "evaluator_outcomes_read": False,
            "metrics_calculated": False,
            "missing_as_zero": False,
            "rna_state": "structurally_missing",
            "champion_claim_allowed": False,
        },
    )
    freeze_tree(
        output,
        {
            "artifact_class": "peakvi_gse244832_target_cell_context",
            "dataset_id": DATASET_ID,
            "context_role": role,
            "native_output_family": "joint_representation",
            "eligible_donor_lineage_units": 48,
            "fragment_cut_sites": "start_plus_4_and_end_minus_5",
            "rna_assay_read": False,
            "cross_assay_join_consumed": False,
            "condition_values_read": False,
            "evaluator_outcomes_accessed": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )


def validate_config(root: Path, config_path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    config = read_json(config_path.resolve(strict=True))
    if (
        config.get("schema_version") != SCHEMA
        or config.get("source_dataset_id") != "gse296875"
        or config.get("target_dataset_id") != DATASET_ID
        or config.get("cell_selection", {}).get("maximum_per_donor_lineage") != 64
        or config.get("cell_selection", {}).get("seed") != 20260825
        or config.get("cell_selection", {}).get("assay") != "atac_only"
        or config.get("cell_selection", {}).get("mask_before_selection") is not True
        or config.get("native_family_contract", {}).get(
            "source_input_reused_without_rematerialization"
        )
        is not True
        or config.get("native_family_contract", {}).get(
            "source_fit_states_reused_without_refit"
        )
        is not True
    ):
        raise PeakVICellMaterializationError("campaign identity or cell selection differs")
    for field in (
        "condition_values_read", "phenotype_values_read", "rna_assay_read",
        "cross_assay_join_consumed", "evaluator_outcomes_read", "sealed_data_read",
        "metrics_calculated", "missing_as_zero",
    ):
        if config.get("firewall", {}).get(field) is not False:
            raise PeakVICellMaterializationError(f"campaign firewall opened: {field}")
    if (
        config.get("firewall", {}).get("ineligible_cells_enter_transform") is not False
        or config.get("firewall", {}).get("scoring_authorized") is not False
    ):
        raise PeakVICellMaterializationError("campaign execution firewall differs")
    paths = {
        "registration": verify_tree(
            root, config["execution_registration"],
            "gse244832_observed_atac_execution_registration",
        ),
        "axis": verify_tree(
            root, config["exchange_axis"],
            "gse244832_reference_guarded_atac_exchange_axis",
        ),
        "source": verify_tree(
            root, config["target_fragment_source"],
            "atac_transport_source_admission_audit",
        ),
        "query_valid": verify_tree(
            root, config["target_queries"]["valid"],
            "gse244832_label_free_query_atac",
        ),
        "query_test": verify_tree(
            root, config["target_queries"]["test"],
            "gse244832_label_free_query_atac",
        ),
        "source_cell_input": verify_manifest(
            root, config["source_cell_input"],
            "peakvi_gse296875_source_cell_input",
        ),
    }
    registration = read_json(paths["registration"] / "execution_contract.json")
    if (
        registration.get("axis_artifacts_sha256")
        != config["exchange_axis"]["artifacts_sha256"]
        or registration.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
        or registration.get("missingness", {}).get("mask_before_any_query_transform")
        is not True
        or registration.get("rna_assay_opened") is not False
        or registration.get("sequence_execution_authorized") is not False
    ):
        raise PeakVICellMaterializationError("execution registration differs")
    if (
        digest(paths["axis"] / "windows.tsv")
        != config["exchange_axis"]["windows_sha256"]
        or digest(paths["axis"] / "label_free_cell_membership.tsv.gz")
        != config["exchange_axis"]["membership_sha256"]
        or digest(
            paths["source"] / config["target_fragment_source"]["inventory_path"]
        )
        != config["target_fragment_source"]["inventory_sha256"]
    ):
        raise PeakVICellMaterializationError("exchange windows differ")
    for role in ROLES:
        query = read_json(paths[f"query_{role}"] / "query_contract.json")
        if (
            query.get("context_role") != role
            or query.get("scored_target_role") != ("test" if role == "valid" else "valid")
            or query.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
            or query.get("rna_assay_opened") is not False
            or query.get("metrics_calculated") is not False
        ):
            raise PeakVICellMaterializationError(f"aggregate query contract differs: {role}")
    source_contract = read_json(paths["source_cell_input"] / "contract.json")
    if (
        source_contract.get("dataset_id") != "gse296875"
        or source_contract.get("cells") != 8121
        or source_contract.get("target_dataset_values_read") is not False
        or source_contract.get("metrics_calculated") is not False
    ):
        raise PeakVICellMaterializationError("exact source revision2 cell input differs")
    reference = Path(config["reference"]["fai_path"]).resolve(strict=True)
    if digest(reference) != config["reference"]["fai_sha256"]:
        raise PeakVICellMaterializationError("reference FAI differs")
    paths["reference"] = reference
    for name, record in config.get("implementation", {}).items():
        path = resolve_inside(root, str(record.get("path", "")), f"implementation {name}")
        if not path.is_file() or digest(path) != record.get("sha256"):
            raise PeakVICellMaterializationError(f"implementation differs: {name}")
    return config, paths


def build(
    *,
    root: Path,
    config_path: Path,
    target_valid_output: Path,
    target_test_output: Path,
    work: Path,
    workers: int,
) -> dict[str, Any]:
    import numpy as np

    if any(path.exists() for path in (target_valid_output, target_test_output, work)):
        raise PeakVICellMaterializationError("output or work path already exists")
    if not 1 <= workers <= 8:
        raise PeakVICellMaterializationError("worker count differs")
    root = root.resolve(strict=True)
    config, paths = validate_config(root, config_path)
    windows = read_windows(paths["axis"] / "windows.tsv")
    reference = read_reference_lengths(paths["reference"])
    registration = read_json(paths["registration"] / "execution_contract.json")
    mask_path = paths["registration"] / registration["query_lineage_missing_state_path"]
    if digest(mask_path) != registration["query_lineage_missing_state_sha256"]:
        raise PeakVICellMaterializationError("registered missingness hash differs")
    unit_rows, unit_states = read_unit_axis(
        paths["registration"] / "query_axis_with_missingness.tsv"
    )
    registered_states = np.load(mask_path, mmap_mode="r", allow_pickle=False)
    if not np.array_equal(unit_states, registered_states):
        raise PeakVICellMaterializationError("registered TSV and binary masks differ")
    cells = select_cells(
        read_membership(paths["axis"] / "label_free_cell_membership.tsv.gz"),
        unit_rows,
        unit_states,
        cap=int(config["cell_selection"]["maximum_per_donor_lineage"]),
        seed=int(config["cell_selection"]["seed"]),
    )
    fragment_root = Path(config["fragment_root"]).resolve(strict=True)
    inventory = read_inventory(
        paths["source"] / config["target_fragment_source"]["inventory_path"],
        fragment_root,
    )
    by_donor: dict[str, dict[str, int]] = defaultdict(dict)
    for row_index, row in enumerate(cells):
        by_donor[row["donor_id"]][row["raw_barcode"]] = row_index
    work.mkdir(parents=True, mode=0o750)
    jobs: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for donor in DONORS:
            record = inventory[donor]
            shard = work / f"target-{donor}"
            future = executor.submit(
                process_fragment_file,
                donor_id=donor,
                path=str(Path(record["path"]).resolve(strict=True)),
                expected_size=int(record["size_bytes"]),
                expected_sha256=record["sha256"],
                selected_lookup=by_donor[donor],
                windows=windows,
                reference_lengths=reference,
                output=str(shard),
                n_rows=len(cells),
            )
            futures[future] = donor
        for future in as_completed(futures):
            jobs.append(future.result())
    shards = [work / f"target-{donor}" for donor in DONORS]
    matrices = {role: merge_shards(shards, role, len(cells)) for role in ROLES}
    if any((matrix.getnnz(axis=1) == 0).any() for matrix in matrices.values()):
        raise PeakVICellMaterializationError("selected target cell has zero context support")
    authorities = {
        "campaign_sha256": digest(config_path.resolve(strict=True)),
        "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
        "exchange_axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
        "target_fragment_source_artifacts_sha256": config["target_fragment_source"]["artifacts_sha256"],
        "source_cell_input_artifacts_sha256": config["source_cell_input"]["artifacts_sha256"],
        "query_lineage_missing_state_sha256": registration["query_lineage_missing_state_sha256"],
        "reference_fai_sha256": config["reference"]["fai_sha256"],
        "fragment_sha256_by_donor": {
            row["donor_id"]: inventory[row["donor_id"]]["sha256"] for row in jobs
        },
        "condition_or_phenotype_authority_included": False,
        "rna_or_cross_assay_authority_included": False,
        "evaluator_outcome_authority_included": False,
    }
    write_target_output(
        target_valid_output,
        "valid",
        cells,
        windows["valid"],
        matrices["valid"],
        paths["query_valid"],
        unit_rows,
        unit_states,
        mask_path,
        {**authorities, "aggregate_query_artifacts_sha256": config["target_queries"]["valid"]["artifacts_sha256"]},
    )
    write_target_output(
        target_test_output,
        "test",
        cells,
        windows["test"],
        matrices["test"],
        paths["query_test"],
        unit_rows,
        unit_states,
        mask_path,
        {**authorities, "aggregate_query_artifacts_sha256": config["target_queries"]["test"]["artifacts_sha256"]},
    )
    receipt = {
        "schema_version": SCHEMA,
        "status": "passed",
        "target_cells": len(cells),
        "target_donors": 18,
        "eligible_donor_lineage_units": 48,
        "structurally_missing_units": 4,
        "below_qc_units": 20,
        "features_per_rotation": 16_000,
        "target_roles_frozen_in_separate_artifacts": True,
        "registered_mask_applied_before_cell_selection": True,
        "ineligible_cells_entered_transform": False,
        "fragment_cut_sites": "start_plus_4_and_end_minus_5",
        "cell_hashes_globally_unique": True,
        "native_output_family": "joint_representation",
        "condition_values_read": False,
        "phenotype_values_read": False,
        "rna_assay_read": False,
        "cross_assay_join_consumed": False,
        "evaluator_outcomes_read": False,
        "metrics_calculated": False,
        "target_valid_artifacts_sha256": digest(target_valid_output / "ARTIFACTS.json"),
        "target_test_artifacts_sha256": digest(target_test_output / "ARTIFACTS.json"),
        "fragment_shards": {row["donor_id"]: row for row in sorted(jobs, key=lambda value: value["donor_id"])},
    }
    write_json_exclusive(work / "receipt.json", receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--target-valid-output", required=True, type=Path)
    parser.add_argument("--target-test-output", required=True, type=Path)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    arguments = parser.parse_args()
    print(
        canonical_json(
            build(
                root=arguments.root,
                config_path=arguments.config,
                target_valid_output=arguments.target_valid_output,
                target_test_output=arguments.target_test_output,
                work=arguments.work,
                workers=arguments.workers,
            )
        )
    )


if __name__ == "__main__":
    main()
