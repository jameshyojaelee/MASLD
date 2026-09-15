#!/usr/bin/env python3
"""Build outcome-blind cell-by-window inputs for PeakVI cross-cohort transfer."""

from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import dataclass
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-peakvi-cross-cohort-cell-materialization-v1"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
ROLES = ("valid", "test")


class PeakVICellMaterializationError(ValueError):
    """Raised when PeakVI cell inputs does not meet the frozen exchange requirements."""


@dataclass(frozen=True)
class WindowIndex:
    rows: tuple[tuple[str, int, int, str], ...]
    by_contig: Mapping[str, tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]]


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def resolve_inside(root: Path, value: str, label: str, *, file: bool = False) -> Path:
    path = (root / value).resolve(strict=True)
    try:
        path.relative_to(root)
    except ValueError as error:
        raise PeakVICellMaterializationError(f"{label} escapes benchmark root") from error
    if file and (path.is_symlink() or not path.is_file()):
        raise PeakVICellMaterializationError(f"{label} is not a regular file")
    return path


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise PeakVICellMaterializationError(f"JSON object required: {path}")
    return value


def verify_authority(root: Path, record: Mapping[str, Any], artifact_class: str) -> Path:
    path = resolve_inside(root, str(record.get("path", "")), artifact_class)
    manifest = verify_frozen_tree(path)
    if (
        digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise PeakVICellMaterializationError(f"{artifact_class} authority differs")
    return path


def verify_manifest_authority(
    root: Path, record: Mapping[str, Any], artifact_class: str
) -> Path:
    """Verify a large frozen tree's read-only manifest without rereading its payload."""

    path = resolve_inside(root, str(record.get("path", "")), artifact_class)
    manifest_path = path / "ARTIFACTS.json"
    manifest = read_json(manifest_path)
    if (
        digest(manifest_path) != record.get("artifacts_sha256")
        or manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise PeakVICellMaterializationError(f"{artifact_class} manifest authority differs")
    return path


def read_windows(path: Path) -> dict[str, WindowIndex]:
    by_role: dict[str, list[tuple[str, int, int, str]]] = {role: [] for role in ROLES}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "window_index", "window_id", "role", "contig", "start", "end"
        ):
            raise PeakVICellMaterializationError("window fields differ")
        seen: set[str] = set()
        for row in reader:
            role = row["role"]
            start, end = int(row["start"]), int(row["end"])
            if (
                role not in by_role
                or row["window_id"] in seen
                or start < 0
                or end - start != 1000
            ):
                raise PeakVICellMaterializationError("window identity or geometry differs")
            seen.add(row["window_id"])
            by_role[role].append((row["contig"], start, end, row["window_id"]))
    if any(len(rows) != 16_000 for rows in by_role.values()):
        raise PeakVICellMaterializationError("role-specific window count differs")
    if {row[0] for row in by_role["valid"]} & {row[0] for row in by_role["test"]}:
        raise PeakVICellMaterializationError("rotation contigs overlap")
    result: dict[str, WindowIndex] = {}
    for role, rows in by_role.items():
        grouped: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
        for offset, (contig, start, end, _identifier) in enumerate(rows):
            grouped[contig].append((start, end, offset))
        packed: dict[str, tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = {}
        for contig, values in grouped.items():
            values.sort()
            packed[contig] = tuple(
                tuple(value[column] for value in values) for column in range(3)
            )  # type: ignore[assignment]
        result[role] = WindowIndex(tuple(rows), packed)
    return result


def windows_at(index: WindowIndex, contig: str, position: int) -> Iterable[int]:
    packed = index.by_contig.get(contig)
    if packed is None:
        return ()
    starts, ends, offsets = packed
    cursor = bisect_right(starts, position) - 1
    hits: list[int] = []
    while cursor >= 0 and ends[cursor] > position:
        hits.append(offsets[cursor])
        cursor -= 1
    return hits


def selected_cells(
    rows: Sequence[Mapping[str, str]], *, dataset_id: str, cap: int, seed: int
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    if cap < 8 or seed < 0:
        raise PeakVICellMaterializationError("cell-selection cap or seed differs")
    grouped: dict[tuple[str, str], list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        lineage = row["lineage_id"]
        if row["analysis_role"] == "primary" and lineage in LINEAGES:
            grouped[(row["donor_id"], lineage)].append(row)
    expected_units = 39 * 4 if dataset_id == "gse296875" else 12 * 4
    if len(grouped) != expected_units:
        raise PeakVICellMaterializationError(f"{dataset_id} donor-lineage census differs")
    selected: list[dict[str, Any]] = []
    available: dict[str, int] = {}
    for (donor, lineage), candidates in sorted(grouped.items()):
        ordered = sorted(
            candidates,
            key=lambda row: sha256(
                f"{seed}\0{dataset_id}\0{row['cell_id']}".encode()
            ).hexdigest(),
        )
        key = f"{donor}\0{lineage}"
        available[key] = len(ordered)
        for rank, row in enumerate(ordered[:cap]):
            selected.append(
                {
                    "dataset_id": dataset_id,
                    "cell_id": row["cell_id"],
                    "cell_hash": sha256(
                        (
                            f"peakvi-cross-cohort-v2\0{dataset_id}\0{donor}\0"
                            f"{row['cell_id']}"
                        ).encode()
                    ).hexdigest(),
                    "raw_barcode": row["raw_barcode"],
                    "well_id": row.get("well_id", ""),
                    "donor_id": donor,
                    "lineage_id": lineage,
                    "outer_fold": row["outer_fold"],
                    "selection_rank": rank,
                    "available_nuclei": len(ordered),
                }
            )
    if not selected:
        raise PeakVICellMaterializationError("cell selection is empty")
    return selected, available


def read_source_membership(path: Path) -> list[dict[str, str]]:
    expected = (
        "well_id", "raw_barcode", "cell_id", "donor_id", "source_label",
        "lineage_id", "analysis_role", "outer_fold",
    )
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != expected:
            raise PeakVICellMaterializationError("GSE296875 membership fields differ")
        rows = list(reader)
    if len(rows) != 68_398:
        raise PeakVICellMaterializationError("GSE296875 membership census differs")
    return rows


def read_target_membership(path: Path) -> list[dict[str, str]]:
    expected = (
        "dataset_id", "donor_id", "raw_barcode", "lineage_id", "analysis_role", "outer_fold"
    )
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != expected:
            raise PeakVICellMaterializationError("GSE281367 membership fields differ")
        for row in reader:
            if row["dataset_id"] != "gse281367":
                raise PeakVICellMaterializationError("GSE281367 membership dataset differs")
            row["cell_id"] = row["raw_barcode"]
            row["well_id"] = row["donor_id"]
            rows.append(row)
    if len({row["donor_id"] for row in rows}) != 12:
        raise PeakVICellMaterializationError("GSE281367 donor census differs")
    return rows


def read_reference_lengths(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0] in result or int(fields[1]) <= 0:
                raise PeakVICellMaterializationError("reference FAI differs")
            result[fields[0]] = int(fields[1])
    if not all(f"chr{value}" in result for value in [*range(1, 23), "X", "Y"]):
        raise PeakVICellMaterializationError("reference primary contigs differ")
    return result


def process_fragment_file(
    *,
    path: str,
    expected_size: int,
    expected_sha256: str,
    selected_lookup: Mapping[str, int],
    windows: Mapping[str, WindowIndex],
    reference_lengths: Mapping[str, int],
    output: str,
    n_rows: int,
) -> dict[str, Any]:
    """Write one sparse shard while checking compressed SHA-256 and gzip CRC in one pass."""

    import numpy as np
    from scipy import sparse

    source = Path(path).resolve(strict=True)
    if source.is_symlink() or source.stat().st_size != expected_size:
        raise PeakVICellMaterializationError(f"fragment path or size differs: {source}")
    role_rows: dict[str, list[int]] = {role: [] for role in ROLES}
    role_columns: dict[str, list[int]] = {role: [] for role in ROLES}
    role_values: dict[str, list[int]] = {role: [] for role in ROLES}
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
        reader = DigestReader(raw)
        with gzip.GzipFile(fileobj=reader, mode="rb") as compressed:
            handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
            for line_number, line in enumerate(handle, start=1):
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise PeakVICellMaterializationError(f"fragment width differs: {source}:{line_number}")
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
                for position in (start, end - 1):
                    for role in ROLES:
                        for column in windows_at(windows[role], contig, position):
                            role_rows[role].append(row_index)
                            role_columns[role].append(column)
                            role_values[role].append(1)
    if compressed_digest.hexdigest() != expected_sha256:
        raise PeakVICellMaterializationError(f"fragment SHA-256 differs: {source}")
    if len(seen) != len(selected_lookup):
        raise PeakVICellMaterializationError(
            f"selected barcodes absent from fragment file: {source} "
            f"expected={len(selected_lookup)} observed={len(seen)}"
        )
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=False)
    identities: dict[str, Any] = {}
    for role in ROLES:
        matrix = sparse.coo_matrix(
            (
                np.asarray(role_values[role], dtype=np.uint16),
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
        "path": str(source),
        "selected_barcodes": len(selected_lookup),
        "selected_barcodes_seen": len(seen),
        "selected_fragment_rows": selected_fragment_rows,
        "gzip_crc_verified_to_eof": True,
        "compressed_sha256_verified": True,
        "matrices": identities,
    }


def merge_shards(shards: Sequence[Path], role: str, n_rows: int) -> Any:
    from scipy import sparse

    result = sparse.csr_matrix((n_rows, 16_000), dtype="uint32")
    for shard in shards:
        matrix = sparse.load_npz(shard / f"{role}.counts.npz").tocsr().astype("uint32")
        if matrix.shape != result.shape:
            raise PeakVICellMaterializationError("sparse shard shape differs")
        result = result + matrix
    result.sum_duplicates()
    result.sort_indices()
    return result


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_axis(path: Path, cells: Sequence[Mapping[str, Any]]) -> None:
    write_tsv(
        path,
        (
            "cell_index", "cell_hash", "donor_id", "lineage_id", "outer_fold",
            "selection_rank", "available_nuclei",
        ),
        (
            {
                "cell_index": index,
                "cell_hash": row["cell_hash"],
                "donor_id": row["donor_id"],
                "lineage_id": row["lineage_id"],
                "outer_fold": row["outer_fold"],
                "selection_rank": row["selection_rank"],
                "available_nuclei": row["available_nuclei"],
            }
            for index, row in enumerate(cells)
        ),
    )


def write_windows(path: Path, index: WindowIndex) -> None:
    write_tsv(
        path,
        ("feature_index", "window_id", "contig", "start", "end"),
        (
            {
                "feature_index": offset,
                "window_id": row[3],
                "contig": row[0],
                "start": row[1],
                "end": row[2],
            }
            for offset, row in enumerate(index.rows)
        ),
    )


def write_source_output(
    output: Path,
    cells: Sequence[Mapping[str, Any]],
    windows: Mapping[str, WindowIndex],
    matrices: Mapping[str, Any],
    authorities: Mapping[str, Any],
) -> None:
    from scipy import sparse

    output.mkdir(parents=True, mode=0o750)
    write_axis(output / "cell_axis.tsv", cells)
    for role in ROLES:
        sparse.save_npz(output / f"{role}.counts.npz", matrices[role], compressed=True)
        write_windows(output / f"{role}.windows.tsv", windows[role])
    write_json_exclusive(output / "input_authorities.json", authorities)
    write_json_exclusive(
        output / "contract.json",
        {
            "schema_version": "masld-bench-peakvi-source-cell-input-v1",
            "dataset_id": "gse296875",
            "cells": len(cells),
            "donors": 39,
            "lineages": list(LINEAGES),
            "roles": list(ROLES),
            "features_per_role": 16_000,
            "matrix_unit": "deduplicated_tn5_cut_site_count_per_selected_single_nucleus",
            "native_peakvi_transform": "binary_x_greater_than_zero_applied_only_in_model_adapter",
            "selection_uses_condition_or_outcome": False,
            "target_dataset_values_read": False,
            "evaluator_outcomes_read": False,
            "metrics_calculated": False,
            "champion_claim_allowed": False,
        },
    )
    freeze_tree(
        output,
        {
            "artifact_class": "peakvi_gse296875_source_cell_input",
            "dataset_id": "gse296875",
            "native_output_family": "joint_representation",
            "cells": len(cells),
            "sealed_outcomes_accessed": False,
            "status": "passed",
        },
    )


def write_target_output(
    output: Path,
    role: str,
    cells: Sequence[Mapping[str, Any]],
    windows: WindowIndex,
    matrix: Any,
    aggregate_query: Path,
    authorities: Mapping[str, Any],
) -> None:
    import numpy as np
    from scipy import sparse

    output.mkdir(parents=True, mode=0o750)
    write_axis(output / "cell_axis.tsv", cells)
    sparse.save_npz(output / "context.counts.npz", matrix, compressed=True)
    write_windows(output / "context.windows.tsv", windows)
    axis = np.load(aggregate_query / "query_unit_fragment_total.uint32.npy", mmap_mode="r", allow_pickle=False)
    if axis.shape != (12, 4, 16_000):
        raise PeakVICellMaterializationError("aggregate query geometry differs")
    donor_index = {f"Z{index + 1:02d}": index for index in range(12)}
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    selected_total = np.zeros((12, 4, 16_000), dtype=np.uint64)
    for row_index, row in enumerate(cells):
        selected_total[donor_index[row["donor_id"]], lineage_index[row["lineage_id"]]] += matrix.getrow(row_index).toarray()[0]
    if np.any(selected_total > np.asarray(axis, dtype=np.uint64)):
        raise PeakVICellMaterializationError("selected-cell counts exceed frozen aggregate query")
    write_json_exclusive(output / "input_authorities.json", authorities)
    write_json_exclusive(
        output / "contract.json",
        {
            "schema_version": "masld-bench-peakvi-target-cell-context-v1",
            "dataset_id": "gse281367",
            "context_role": role,
            "scored_target_role": "test" if role == "valid" else "valid",
            "cells": len(cells),
            "donors": 12,
            "lineages": list(LINEAGES),
            "features": 16_000,
            "matrix_unit": "deduplicated_tn5_cut_site_count_per_selected_single_nucleus",
            "native_peakvi_transform": "binary_x_greater_than_zero_applied_only_in_model_adapter",
            "selected_cell_counts_bounded_by_frozen_aggregate_query": True,
            "other_rotation_values_present": False,
            "scored_target_contigs_present": False,
            "condition_values_read": False,
            "phenotype_values_read": False,
            "evaluator_outcomes_read": False,
            "metrics_calculated": False,
            "rna_state": "structurally_missing",
            "champion_claim_allowed": False,
        },
    )
    freeze_tree(
        output,
        {
            "artifact_class": "peakvi_gse281367_target_cell_context",
            "dataset_id": "gse281367",
            "context_role": role,
            "native_output_family": "joint_representation",
            "condition_values_read": False,
            "evaluator_outcomes_accessed": False,
            "status": "passed",
        },
    )


def validate_config(root: Path, path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    config = read_json(path.resolve(strict=True))
    if (
        config.get("schema_version") != SCHEMA
        or config.get("source_dataset_id") != "gse296875"
        or config.get("target_dataset_id") != "gse281367"
        or config.get("cell_selection", {}).get("uses_condition_or_outcome") is not False
    ):
        raise PeakVICellMaterializationError("campaign identity or selection firewall differs")
    for field in (
        "condition_values_read", "phenotype_values_read", "evaluator_outcomes_read",
        "sealed_data_read", "metrics_calculated",
    ):
        if config.get("firewall", {}).get(field) is not False:
            raise PeakVICellMaterializationError(f"campaign firewall opened: {field}")
    paths = {
        "gate": verify_authority(root, config["terminal_gate"], "gse296875_observed_multiome_specialist_gate_readiness"),
        "axis": verify_authority(root, config["exchange_axis"], "gse281367_label_free_atac_exchange_axis"),
        "source_membership": verify_authority(root, config["source_membership"], "gse296875_fragment_membership"),
        "source_raw": verify_manifest_authority(
            root, config["source_raw_fragments"], "gse296875_raw_fragment_acquisition"
        ),
        "target_source": verify_authority(root, config["target_source"], "atac_transport_source_admission_audit"),
        "target_valid": verify_authority(root, config["target_queries"]["valid"], "gse281367_label_free_query_atac"),
        "target_test": verify_authority(root, config["target_queries"]["test"], "gse281367_label_free_query_atac"),
    }
    if digest(paths["axis"] / "windows.tsv") != config["exchange_axis"]["windows_sha256"]:
        raise PeakVICellMaterializationError("exchange window hash differs")
    reference = Path(config["reference"]["fai_path"]).resolve(strict=True)
    if digest(reference) != config["reference"]["fai_sha256"]:
        raise PeakVICellMaterializationError("reference FAI hash differs")
    paths["reference"] = reference
    implementation = config["implementation"]
    for name in ("builder", "test", "sbatch"):
        source = resolve_inside(root, implementation[name]["path"], f"implementation {name}", file=True)
        if digest(source) != implementation[name]["sha256"]:
            raise PeakVICellMaterializationError(f"implementation hash differs: {name}")
    return config, paths


def read_artifact_records(path: Path) -> dict[str, dict[str, Any]]:
    manifest = read_json(path / "ARTIFACTS.json")
    return {str(row["path"]): row for row in manifest["artifacts"]}


def build(
    *,
    root: Path,
    config_path: Path,
    source_output: Path,
    target_valid_output: Path,
    target_test_output: Path,
    work: Path,
    workers: int,
) -> dict[str, Any]:
    if any(path.exists() for path in (source_output, target_valid_output, target_test_output, work)):
        raise PeakVICellMaterializationError("output or work path already exists")
    if not 1 <= workers <= 16:
        raise PeakVICellMaterializationError("worker count differs")
    config, paths = validate_config(root.resolve(strict=True), config_path)
    windows = read_windows(paths["axis"] / "windows.tsv")
    reference = read_reference_lengths(paths["reference"])
    cap = int(config["cell_selection"]["maximum_per_donor_lineage"])
    seed = int(config["cell_selection"]["seed"])
    source_cells, _source_available = selected_cells(
        read_source_membership(paths["source_membership"] / "cell_membership.tsv.gz"),
        dataset_id="gse296875", cap=cap, seed=seed,
    )
    target_cells, _target_available = selected_cells(
        read_target_membership(paths["axis"] / "label_free_cell_membership.tsv.gz"),
        dataset_id="gse281367", cap=cap, seed=seed,
    )
    work.mkdir(parents=True, mode=0o750)
    source_records = read_artifact_records(paths["source_raw"])
    source_filelist: list[tuple[str, Path, dict[str, Any]]] = []
    pattern = re.compile(r"^fragments/GSM[0-9]+_(well[1-8])_atac_fragments[.]tsv[.]gz$")
    for relative, record in sorted(source_records.items()):
        match = pattern.fullmatch(relative)
        if match:
            source_filelist.append((match.group(1), paths["source_raw"] / relative, record))
    if sorted(well for well, _path, _record in source_filelist) != [f"well{i}" for i in range(1, 9)]:
        raise PeakVICellMaterializationError("source well roster differs")
    target_inventory: dict[str, dict[str, str]] = {}
    with (paths["target_source"] / "source_inventory.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            candidate = Path(row["path"])
            if len(candidate.parts) >= 3 and candidate.parts[-2:] == ("outs", "fragments.tsv.gz"):
                donor = candidate.parts[-3]
                if donor.startswith("Z"):
                    target_inventory[donor] = row
    if sorted(target_inventory) != [f"Z{i:02d}" for i in range(1, 13)]:
        raise PeakVICellMaterializationError("target fragment roster differs")
    source_by_well: dict[str, dict[str, int]] = defaultdict(dict)
    for index, row in enumerate(source_cells):
        source_by_well[row["well_id"]][row["raw_barcode"]] = index
    target_by_donor: dict[str, dict[str, int]] = defaultdict(dict)
    for index, row in enumerate(target_cells):
        target_by_donor[row["donor_id"]][row["raw_barcode"]] = index
    jobs: list[tuple[str, dict[str, Any]]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for well, path, record in source_filelist:
            shard = work / f"source-{well}"
            future = executor.submit(
                process_fragment_file,
                path=str(path), expected_size=int(record["size_bytes"]),
                expected_sha256=str(record["sha256"]), selected_lookup=source_by_well[well],
                windows=windows, reference_lengths=reference, output=str(shard),
                n_rows=len(source_cells),
            )
            futures[future] = (f"source:{well}", shard)
        for donor, record in sorted(target_inventory.items()):
            shard = work / f"target-{donor}"
            future = executor.submit(
                process_fragment_file,
                path=row_path(record), expected_size=int(record["size_bytes"]),
                expected_sha256=record["sha256"], selected_lookup=target_by_donor[donor],
                windows=windows, reference_lengths=reference, output=str(shard),
                n_rows=len(target_cells),
            )
            futures[future] = (f"target:{donor}", shard)
        for future in as_completed(futures):
            label, shard = futures[future]
            jobs.append((label, {"shard": str(shard), **future.result()}))
    source_shards = [work / f"source-well{i}" for i in range(1, 9)]
    target_shards = [work / f"target-Z{i:02d}" for i in range(1, 13)]
    source_matrices = {
        role: merge_shards(source_shards, role, len(source_cells)) for role in ROLES
    }
    target_matrices = {
        role: merge_shards(target_shards, role, len(target_cells)) for role in ROLES
    }
    for label, matrices in (("source", source_matrices), ("target", target_matrices)):
        if any((matrix.getnnz(axis=1) == 0).any() for matrix in matrices.values()):
            raise PeakVICellMaterializationError(f"{label} selected cell has zero context support")
    authorities = {
        "campaign_sha256": digest(config_path.resolve(strict=True)),
        "terminal_gate_artifacts_sha256": config["terminal_gate"]["artifacts_sha256"],
        "exchange_axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
        "source_membership_artifacts_sha256": config["source_membership"]["artifacts_sha256"],
        "source_raw_fragments_artifacts_sha256": config["source_raw_fragments"]["artifacts_sha256"],
        "target_source_artifacts_sha256": config["target_source"]["artifacts_sha256"],
        "reference_fai_sha256": config["reference"]["fai_sha256"],
        "condition_or_phenotype_authority_included": False,
        "evaluator_outcome_authority_included": False,
    }
    write_source_output(source_output, source_cells, windows, source_matrices, authorities)
    write_target_output(
        target_valid_output, "valid", target_cells, windows["valid"], target_matrices["valid"],
        paths["target_valid"], {**authorities, "aggregate_query_artifacts_sha256": config["target_queries"]["valid"]["artifacts_sha256"]},
    )
    write_target_output(
        target_test_output, "test", target_cells, windows["test"], target_matrices["test"],
        paths["target_test"], {**authorities, "aggregate_query_artifacts_sha256": config["target_queries"]["test"]["artifacts_sha256"]},
    )
    receipt = {
        "schema_version": SCHEMA,
        "status": "passed",
        "source_cells": len(source_cells),
        "target_cells": len(target_cells),
        "source_donors": 39,
        "target_donors": 12,
        "lineages": list(LINEAGES),
        "features_per_rotation": 16_000,
        "source_roles_materialized": list(ROLES),
        "target_roles_frozen_in_separate_artifacts": True,
        "target_contigs_absent_from_each_rotation_artifact": True,
        "native_output_family": "joint_representation",
        "source_only_decoder_required": True,
        "condition_values_read": False,
        "phenotype_values_read": False,
        "evaluator_outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
        "source_artifacts_sha256": digest(source_output / "ARTIFACTS.json"),
        "target_valid_artifacts_sha256": digest(target_valid_output / "ARTIFACTS.json"),
        "target_test_artifacts_sha256": digest(target_test_output / "ARTIFACTS.json"),
        "fragment_shards": {label: record for label, record in sorted(jobs)},
    }
    write_json_exclusive(work / "receipt.json", receipt)
    return receipt


def row_path(record: Mapping[str, str]) -> str:
    return str(Path(record["path"]).resolve(strict=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--source-output", required=True, type=Path)
    parser.add_argument("--target-valid-output", required=True, type=Path)
    parser.add_argument("--target-test-output", required=True, type=Path)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    arguments = parser.parse_args()
    result = build(
        root=arguments.root,
        config_path=arguments.config,
        source_output=arguments.source_output,
        target_valid_output=arguments.target_valid_output,
        target_test_output=arguments.target_test_output,
        work=arguments.work,
        workers=arguments.workers,
    )
    print(canonical_json(result))


if __name__ == "__main__":
    main()
