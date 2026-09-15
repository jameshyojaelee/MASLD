#!/usr/bin/env python3
"""Audit local GSE244832 and GSE281367 ATAC transport source assets.

This audit does not activate either dataset. It inventories read-only source
objects, checks donor and cell-type axes in the processed H5ADs, and rejects
legacy donor-pseudobulk matrices with duplicate or malformed peak columns.
Benchmark outcomes must subsequently be rebuilt from the raw fragment files on
one frozen project window set.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence


LINEAGES = ("cholangiocyte", "hepatocyte", "macrophage", "stellate")
PEAK_RE = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)):(\d+)-(\d+)$")


class ATACSourceAuditError(RuntimeError):
    """Raised when a local ATAC source does not meet the audit requirements."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_peak_names(names: Sequence[str]) -> dict[str, Any]:
    counts = Counter(names)
    duplicate_instances = sum(value - 1 for value in counts.values())
    malformed = 0
    invalid_intervals = 0
    contigs: Counter[str] = Counter()
    for name in counts:
        match = PEAK_RE.fullmatch(name)
        if match is None:
            malformed += 1
            continue
        start, end = int(match.group(2)), int(match.group(3))
        if start < 0 or end <= start:
            invalid_intervals += 1
        contigs[match.group(1)] += 1
    return {
        "feature_columns": len(names),
        "unique_feature_columns": len(counts),
        "duplicate_feature_instances": duplicate_instances,
        "duplicate_feature_ids": sum(value > 1 for value in counts.values()),
        "malformed_feature_ids": malformed,
        "invalid_intervals": invalid_intervals,
        "contigs": dict(sorted(contigs.items())),
    }


def _read_coldata(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"donor_id", "libsize"}
        if not required <= set(reader.fieldnames or ()):
            raise ATACSourceAuditError(f"coldata fields differ: {path}")
        rows = [dict(row) for row in reader]
    by_donor = {row["donor_id"]: row for row in rows}
    if len(by_donor) != len(rows) or not rows:
        raise ATACSourceAuditError(f"coldata donor IDs are empty or duplicated: {path}")
    return by_donor


def audit_pseudobulk(counts_path: Path, coldata_path: Path) -> dict[str, Any]:
    coldata = _read_coldata(coldata_path)
    with gzip.open(counts_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ATACSourceAuditError(f"empty pseudobulk matrix: {counts_path}") from error
        if len(header) < 2 or header[0] != "donor_id":
            raise ATACSourceAuditError(f"pseudobulk header differs: {counts_path}")
        features = header[1:]
        feature_audit = audit_peak_names(features)
        row_ids: list[str] = []
        row_sums: dict[str, int] = {}
        minimum = math.inf
        maximum = -math.inf
        nonzero = 0
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ATACSourceAuditError(
                    f"pseudobulk row width differs at {counts_path}:{row_number}"
                )
            donor = row[0]
            if donor in row_sums:
                raise ATACSourceAuditError(f"duplicate pseudobulk donor: {donor}")
            values: list[int] = []
            for column, value_text in enumerate(row[1:], start=2):
                try:
                    value = int(value_text)
                except ValueError as error:
                    raise ATACSourceAuditError(
                        f"non-integer count at {counts_path}:{row_number}:{column}"
                    ) from error
                if value < 0:
                    raise ATACSourceAuditError(
                        f"negative count at {counts_path}:{row_number}:{column}"
                    )
                values.append(value)
            row_ids.append(donor)
            row_sums[donor] = sum(values)
            if values:
                minimum = min(minimum, min(values))
                maximum = max(maximum, max(values))
                nonzero += sum(value != 0 for value in values)
    if set(row_ids) != set(coldata):
        raise ATACSourceAuditError(f"pseudobulk and coldata donor axes differ: {counts_path}")
    libsize_mismatches = {
        donor: {"matrix": row_sums[donor], "coldata": int(coldata[donor]["libsize"])}
        for donor in row_ids
        if row_sums[donor] != int(coldata[donor]["libsize"])
    }
    usable = (
        feature_audit["duplicate_feature_instances"] == 0
        and feature_audit["malformed_feature_ids"] == 0
        and feature_audit["invalid_intervals"] == 0
        and not libsize_mismatches
    )
    return {
        **feature_audit,
        "donor_rows": len(row_ids),
        "donor_ids": sorted(row_ids),
        "numeric_values": len(row_ids) * len(features),
        "nonzero_values": nonzero,
        "minimum_value": None if math.isinf(minimum) else minimum,
        "maximum_value": None if math.isinf(maximum) else maximum,
        "libsize_mismatches": libsize_mismatches,
        "usable_as_benchmark_outcome": usable,
    }


def _decode_h5_values(values: Any) -> list[str]:
    result: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            result.append(value.decode("utf-8"))
        else:
            result.append(str(value))
    return result


def _read_h5_obs_column(obs: Any, key: str) -> list[str]:
    import numpy as np

    if key not in obs:
        raise ATACSourceAuditError(f"H5AD obs column is missing: {key}")
    value = obs[key]
    if hasattr(value, "shape"):
        return _decode_h5_values(value[...])
    if "categories" not in value or "codes" not in value:
        raise ATACSourceAuditError(f"unsupported H5AD obs encoding: {key}")
    categories = _decode_h5_values(value["categories"][...])
    codes = np.asarray(value["codes"][...], dtype=int)
    return ["" if code < 0 else categories[code] for code in codes]


def audit_h5ad(path: Path, expected_donors: set[str]) -> tuple[dict[str, Any], list[dict[str, str]]]:
    import h5py

    with h5py.File(path, "r") as handle:
        if "obs" not in handle:
            raise ATACSourceAuditError(f"H5AD lacks obs: {path}")
        obs = handle["obs"]
        index_key = obs.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode("utf-8")
        cell_ids = _read_h5_obs_column(obs, str(index_key))
        donors = _read_h5_obs_column(obs, "donor_id")
        cell_types = _read_h5_obs_column(obs, "cell_type")
        conditions = (
            _read_h5_obs_column(obs, "condition")
            if "condition" in obs
            else [""] * len(cell_ids)
        )
        obs_columns = sorted(str(key) for key in obs.keys())
    lengths = {len(cell_ids), len(donors), len(cell_types), len(conditions)}
    if len(lengths) != 1 or not cell_ids:
        raise ATACSourceAuditError(f"H5AD obs axes differ or are empty: {path}")
    observed_donors = set(donors)
    if observed_donors != expected_donors:
        raise ATACSourceAuditError(
            f"H5AD donor axis differs for {path.name}: {sorted(observed_donors)}"
        )
    rows = [
        {
            "source_cell_id": cell_id,
            "donor_id": donor,
            "source_cell_type": cell_type,
            "condition": condition,
        }
        for cell_id, donor, cell_type, condition in zip(
            cell_ids, donors, cell_types, conditions, strict=True
        )
    ]
    report = {
        "cells": len(rows),
        "donors": len(observed_donors),
        "donor_counts": dict(sorted(Counter(donors).items())),
        "cell_type_counts": dict(sorted(Counter(cell_types).items())),
        "condition_counts": dict(sorted(Counter(conditions).items())),
        "obs_columns": obs_columns,
        "cell_ids_unique": len(set(cell_ids)) == len(cell_ids),
        "missing_cell_type": sum(not value for value in cell_types),
        "missing_donor": sum(not value for value in donors),
    }
    return report, rows


def read_gse281_metadata(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = {row["donor_id"]: row["condition"] for row in rows}
    if len(result) != 12 or Counter(result.values()) != Counter({"NORMAL": 6, "MASH": 6}):
        raise ATACSourceAuditError("GSE281367 donor metadata differs")
    return result


def read_gse244_metadata(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {row["donor_id"]: row["condition"] for row in rows}
    if len(result) != 18:
        raise ATACSourceAuditError("GSE244832 donor metadata differs")
    return result


def inventory_source(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    if path.is_symlink() or not resolved.is_file():
        raise ATACSourceAuditError(f"source must be a regular non-symlink file: {path}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse281-h5ad", type=Path, required=True)
    parser.add_argument("--gse281-meta", type=Path, required=True)
    parser.add_argument("--gse281-fragments", type=Path, required=True)
    parser.add_argument("--gse281-pseudobulk", type=Path, required=True)
    parser.add_argument("--gse244-h5ad", type=Path, required=True)
    parser.add_argument("--gse244-meta", type=Path, required=True)
    parser.add_argument("--gse244-fragments", type=Path, required=True)
    parser.add_argument("--gse244-pseudobulk", type=Path, required=True)
    parser.add_argument("--hash-workers", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists() or not 1 <= args.hash_workers <= 16:
        raise ATACSourceAuditError("invalid output or hash-worker contract")
    gse281_meta = read_gse281_metadata(args.gse281_meta.resolve(strict=True))
    gse244_meta = read_gse244_metadata(args.gse244_meta.resolve(strict=True))
    h5_reports: dict[str, Any] = {}
    h5_reports["gse281367"], gse281_cells = audit_h5ad(
        args.gse281_h5ad.resolve(strict=True), set(gse281_meta)
    )
    h5_reports["gse244832"], gse244_cells = audit_h5ad(
        args.gse244_h5ad.resolve(strict=True), set(gse244_meta)
    )

    fragment_paths = [
        args.gse281_fragments / f"Z{index:02d}" / "outs" / "fragments.tsv.gz"
        for index in range(1, 13)
    ] + [
        args.gse244_fragments / f"D{index:02d}_fragments.tsv.gz"
        for index in range(1, 19)
    ]
    source_paths = [
        args.gse281_h5ad,
        args.gse281_meta,
        args.gse244_h5ad,
        args.gse244_meta,
        *fragment_paths,
    ]
    with ThreadPoolExecutor(max_workers=args.hash_workers) as executor:
        source_inventory = list(executor.map(inventory_source, source_paths))

    pseudobulk_reports: dict[str, Any] = {}
    pseudobulk_rows: list[dict[str, object]] = []
    for dataset_id, directory, suffix in (
        ("gse281367", args.gse281_pseudobulk, "_GSE281367"),
        ("gse244832", args.gse244_pseudobulk, ""),
    ):
        pseudobulk_reports[dataset_id] = {}
        for lineage in LINEAGES:
            tag = "hep" if lineage == "hepatocyte" else lineage
            counts = directory / f"{tag}_pseudobulk_counts{suffix}.tsv.gz"
            coldata = directory / f"{tag}_pseudobulk_coldata{suffix}.tsv"
            report = audit_pseudobulk(counts.resolve(strict=True), coldata.resolve(strict=True))
            pseudobulk_reports[dataset_id][lineage] = report
            pseudobulk_rows.append(
                {
                    "dataset_id": dataset_id,
                    "lineage_id": lineage,
                    "donors": report["donor_rows"],
                    "feature_columns": report["feature_columns"],
                    "unique_feature_columns": report["unique_feature_columns"],
                    "duplicate_feature_instances": report["duplicate_feature_instances"],
                    "malformed_feature_ids": report["malformed_feature_ids"],
                    "libsize_mismatches": len(report["libsize_mismatches"]),
                    "usable_as_benchmark_outcome": str(
                        report["usable_as_benchmark_outcome"]
                    ).lower(),
                }
            )

    args.output.mkdir(mode=0o750, parents=True)
    write_tsv(
        args.output / "source_inventory.tsv",
        ("path", "size_bytes", "sha256"),
        source_inventory,
    )
    write_tsv(
        args.output / "legacy_pseudobulk_audit.tsv",
        (
            "dataset_id",
            "lineage_id",
            "donors",
            "feature_columns",
            "unique_feature_columns",
            "duplicate_feature_instances",
            "malformed_feature_ids",
            "libsize_mismatches",
            "usable_as_benchmark_outcome",
        ),
        pseudobulk_rows,
    )
    celltype_rows = []
    for dataset_id, rows in (("gse281367", gse281_cells), ("gse244832", gse244_cells)):
        counts = Counter((row["donor_id"], row["source_cell_type"]) for row in rows)
        celltype_rows.extend(
            {
                "dataset_id": dataset_id,
                "donor_id": donor,
                "source_cell_type": cell_type,
                "cells": count,
            }
            for (donor, cell_type), count in sorted(counts.items())
        )
    write_tsv(
        args.output / "donor_celltype_census.tsv",
        ("dataset_id", "donor_id", "source_cell_type", "cells"),
        celltype_rows,
    )
    summary = {
        "schema_version": "masld-bench-atac-transport-source-audit-v1",
        "cohort_family_ids": ["gse244832", "gse281367"],
        "h5ad": h5_reports,
        "legacy_pseudobulk": pseudobulk_reports,
        "source_files": len(source_inventory),
        "source_bytes": sum(int(row["size_bytes"]) for row in source_inventory),
        "source_files_full_sha256": True,
        "dataset_activated": False,
        "raw_fragments_are_preferred_source": True,
        "benchmark_outcome_rule": "rebuild_from_raw_fragments_on_one_frozen_project_window_set",
        "remaining_blockers": [
            "freeze final donor-by-cell-type raw-barcode membership against fragment sources",
            "full gzip and fragment-coordinate/reference audit during the one-pass split",
            "freeze GRCh38.p14 window crosswalk and evaluator-only outcomes",
            "freeze GSE244832 donor-level RNA context without false cell pairing or Atlas double counting",
            "bind source-specific rights, exposure, QC, TaskSpec, and outer-study role",
        ],
    }
    (args.output / "audit.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "source_files": summary["source_files"],
                "source_bytes": summary["source_bytes"],
                "gse281_cells": h5_reports["gse281367"]["cells"],
                "gse244_cells": h5_reports["gse244832"]["cells"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
