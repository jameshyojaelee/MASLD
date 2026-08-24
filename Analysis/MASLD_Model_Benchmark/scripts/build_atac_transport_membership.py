#!/usr/bin/env python3
"""Freeze donor/cell-type/raw-barcode membership for two ATAC transport cohorts."""

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
from typing import Any, Iterable, Mapping, Sequence

from scripts.audit_atac_transport_sources import (
    ATACSourceAuditError,
    _read_h5_obs_column,
    read_gse244_metadata,
    read_gse281_metadata,
    sha256_file,
)


SOURCE_AUDIT_CLASS = "atac_transport_source_admission_audit"
PRIMARY_LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
LABEL_CONTRACT: dict[str, dict[str, tuple[str, str]]] = {
    "gse281367": {
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
    },
    "gse244832": {
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
    },
}
MEMBERSHIP_FIELDS = (
    "dataset_id",
    "donor_id",
    "source_cell_id",
    "raw_barcode",
    "source_cell_type",
    "lineage_id",
    "analysis_role",
    "condition",
    "outer_fold",
)


class ATACMembershipError(RuntimeError):
    """Raised when H5AD cell identity cannot be joined to donor raw barcodes."""


class HashingWriter:
    def __init__(self, handle: io.BufferedWriter) -> None:
        self.handle = handle
        self.digest = sha256()
        self.name = str(handle.name)

    def write(self, value: bytes) -> int:
        self.digest.update(value)
        return self.handle.write(value)

    def flush(self) -> None:
        self.handle.flush()

    def tell(self) -> int:
        return self.handle.tell()

    def writable(self) -> bool:
        return True


def open_deterministic_gzip(path: Path) -> tuple[io.TextIOWrapper, HashingWriter, io.BufferedWriter]:
    raw = path.open("xb")
    hashed = HashingWriter(raw)
    compressed = gzip.GzipFile(
        filename="", mode="wb", compresslevel=6, fileobj=hashed, mtime=0
    )
    return io.TextIOWrapper(compressed, encoding="utf-8", newline=""), hashed, raw


def close_deterministic_gzip(
    text: io.TextIOWrapper, hashed: HashingWriter, raw: io.BufferedWriter
) -> str:
    text.close()
    hashed.flush()
    raw.close()
    return hashed.digest.hexdigest()


def fold_index(dataset_id: str, donor_id: str, seed: int = 20260821) -> int:
    if not dataset_id or not donor_id:
        raise ATACMembershipError("empty dataset or donor identifier")
    digest = sha256(f"{seed}\0{dataset_id}\0{donor_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % 5


def resolve_raw_barcodes(merged_ids: Sequence[str], raw_ids: Sequence[str]) -> list[str]:
    raw_set = set(raw_ids)
    if len(raw_set) != len(raw_ids):
        raise ATACMembershipError("per-donor raw barcodes are duplicated")
    resolved: list[str] = []
    used: set[str] = set()
    for merged in merged_ids:
        candidate = merged
        matches: list[str] = []
        while candidate:
            if candidate in raw_set:
                matches.append(candidate)
            suffix = re.fullmatch(r"(.+)-([0-9]+)", candidate)
            if suffix is None:
                break
            candidate = suffix.group(1)
        available = [value for value in matches if value not in used]
        if len(available) != 1:
            raise ATACMembershipError(
                f"merged cell ID has {len(available)} unused raw matches: {merged}"
            )
        raw = available[0]
        used.add(raw)
        resolved.append(raw)
    if used != raw_set:
        raise ATACMembershipError(
            f"merged/per-donor cell universes differ: used={len(used)} raw={len(raw_set)}"
        )
    return resolved


def _read_obs(path: Path, fields: Sequence[str]) -> dict[str, list[str]]:
    import h5py

    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        index_key = obs.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode("utf-8")
        result = {"source_cell_id": _read_h5_obs_column(obs, str(index_key))}
        result.update({field: _read_h5_obs_column(obs, field) for field in fields})
    lengths = {len(value) for value in result.values()}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ATACMembershipError(f"H5AD obs axes differ or are empty: {path}")
    return result


def build_dataset_membership(
    *,
    dataset_id: str,
    merged_h5ad: Path,
    per_donor_dir: Path,
    metadata: Mapping[str, str],
) -> tuple[list[dict[str, str]], list[Path]]:
    if dataset_id not in LABEL_CONTRACT:
        raise ATACMembershipError(f"unsupported dataset: {dataset_id}")
    merged = _read_obs(merged_h5ad, ("donor_id", "cell_type", "condition"))
    rows: list[dict[str, str]] = []
    per_donor_paths: list[Path] = []
    for donor in sorted(metadata):
        indexes = [
            index
            for index, observed in enumerate(merged["donor_id"])
            if observed == donor
        ]
        if not indexes:
            raise ATACMembershipError(f"merged H5AD lacks donor {dataset_id}/{donor}")
        per_donor = (per_donor_dir / f"{donor}.h5ad").resolve(strict=True)
        per_donor_paths.append(per_donor)
        raw_obs = _read_obs(per_donor, ())
        merged_ids = [merged["source_cell_id"][index] for index in indexes]
        raw_ids = raw_obs["source_cell_id"]
        if len(merged_ids) != len(raw_ids):
            raise ATACMembershipError(
                f"cell count differs for {dataset_id}/{donor}: {len(merged_ids)} != {len(raw_ids)}"
            )
        raw_barcodes = resolve_raw_barcodes(merged_ids, raw_ids)
        expected_condition = metadata[donor]
        for index, raw_barcode in zip(indexes, raw_barcodes, strict=True):
            observed_condition = merged["condition"][index]
            if observed_condition != expected_condition:
                raise ATACMembershipError(
                    f"condition differs for {dataset_id}/{donor}: {observed_condition}"
                )
            source_label = merged["cell_type"][index]
            if source_label not in LABEL_CONTRACT[dataset_id]:
                raise ATACMembershipError(
                    f"unregistered source cell type for {dataset_id}: {source_label}"
                )
            lineage, role = LABEL_CONTRACT[dataset_id][source_label]
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "donor_id": donor,
                    "source_cell_id": merged["source_cell_id"][index],
                    "raw_barcode": raw_barcode,
                    "source_cell_type": source_label,
                    "lineage_id": lineage,
                    "analysis_role": role,
                    "condition": expected_condition,
                    "outer_fold": str(fold_index(dataset_id, donor)),
                }
            )
    if len({(row["donor_id"], row["raw_barcode"]) for row in rows}) != len(rows):
        raise ATACMembershipError(f"donor-local raw barcode duplication in {dataset_id}")
    return rows, per_donor_paths


def verify_source_audit(path: Path, expected_sha256: str) -> dict[str, Any]:
    artifacts = path / "ARTIFACTS.json"
    if sha256_file(artifacts) != expected_sha256:
        raise ATACMembershipError("source-audit ARTIFACTS SHA-256 differs")
    manifest = json.loads(artifacts.read_text(encoding="utf-8"))
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("artifact_class") != SOURCE_AUDIT_CLASS
        or metadata.get("cohort_family_ids") != ["gse244832", "gse281367"]
        or metadata.get("activates_dataset") is not False
    ):
        raise ATACMembershipError("source-audit metadata differs")
    return manifest


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--gse281-merged", type=Path, required=True)
    parser.add_argument("--gse281-per-donor", type=Path, required=True)
    parser.add_argument("--gse281-meta", type=Path, required=True)
    parser.add_argument("--gse244-merged", type=Path, required=True)
    parser.add_argument("--gse244-per-donor", type=Path, required=True)
    parser.add_argument("--gse244-meta", type=Path, required=True)
    parser.add_argument("--hash-workers", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists() or not 1 <= args.hash_workers <= 16:
        raise ATACMembershipError("invalid output or hash-worker contract")
    source_audit = args.source_audit.resolve(strict=True)
    verify_source_audit(source_audit, args.source_artifacts_sha256)
    gse281_meta = read_gse281_metadata(args.gse281_meta.resolve(strict=True))
    gse244_meta = read_gse244_metadata(args.gse244_meta.resolve(strict=True))
    rows281, sources281 = build_dataset_membership(
        dataset_id="gse281367",
        merged_h5ad=args.gse281_merged.resolve(strict=True),
        per_donor_dir=args.gse281_per_donor.resolve(strict=True),
        metadata=gse281_meta,
    )
    rows244, sources244 = build_dataset_membership(
        dataset_id="gse244832",
        merged_h5ad=args.gse244_merged.resolve(strict=True),
        per_donor_dir=args.gse244_per_donor.resolve(strict=True),
        metadata=gse244_meta,
    )
    rows = sorted(
        [*rows281, *rows244],
        key=lambda row: (row["dataset_id"], row["donor_id"], row["raw_barcode"]),
    )
    args.output.mkdir(mode=0o750, parents=True)
    membership, hashed, raw = open_deterministic_gzip(args.output / "cell_membership.tsv.gz")
    writer = csv.DictWriter(membership, fieldnames=MEMBERSHIP_FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    membership_sha256 = close_deterministic_gzip(membership, hashed, raw)

    census_rows = []
    for dataset_id in ("gse244832", "gse281367"):
        dataset_rows = [row for row in rows if row["dataset_id"] == dataset_id]
        for donor in sorted({row["donor_id"] for row in dataset_rows}):
            donor_rows = [row for row in dataset_rows if row["donor_id"] == donor]
            for lineage in PRIMARY_LINEAGES:
                selected = [
                    row
                    for row in donor_rows
                    if row["lineage_id"] == lineage and row["analysis_role"] == "primary"
                ]
                census_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "donor_id": donor,
                        "condition": donor_rows[0]["condition"],
                        "outer_fold": donor_rows[0]["outer_fold"],
                        "lineage_id": lineage,
                        "cells": len(selected),
                        "observed": str(bool(selected)).lower(),
                    }
                )
    write_tsv(
        args.output / "donor_lineage_census.tsv",
        (
            "dataset_id",
            "donor_id",
            "condition",
            "outer_fold",
            "lineage_id",
            "cells",
            "observed",
        ),
        census_rows,
    )
    with ThreadPoolExecutor(max_workers=args.hash_workers) as executor:
        hashes = list(executor.map(sha256_file, [*sources244, *sources281]))
    inventory = [
        {
            "dataset_id": "gse244832" if path in sources244 else "gse281367",
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "sha256": digest,
        }
        for path, digest in zip([*sources244, *sources281], hashes, strict=True)
    ]
    write_tsv(
        args.output / "per_donor_h5ad_inventory.tsv",
        ("dataset_id", "path", "size_bytes", "sha256"),
        inventory,
    )
    summary = {
        "schema_version": "masld-bench-atac-transport-membership-v1",
        "source_audit_artifacts_sha256": args.source_artifacts_sha256,
        "pairing_by_dataset": {
            "gse244832": "same_sample_different_aliquot",
            "gse281367": "same_study_unpaired",
        },
        "biological_unit": "donor",
        "cells": Counter(row["dataset_id"] for row in rows),
        "donors": {
            dataset_id: len({row["donor_id"] for row in rows if row["dataset_id"] == dataset_id})
            for dataset_id in ("gse244832", "gse281367")
        },
        "primary_lineages": list(PRIMARY_LINEAGES),
        "membership_sha256": membership_sha256,
        "per_donor_h5ads_full_sha256": True,
        "fragment_membership_verified": False,
        "dataset_activated": False,
        "remaining_gate": "one-pass raw-fragment membership, gzip, coordinate, reference, and fixed-window outcome build",
    }
    (args.output / "membership.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True, default=dict))


if __name__ == "__main__":
    main()
