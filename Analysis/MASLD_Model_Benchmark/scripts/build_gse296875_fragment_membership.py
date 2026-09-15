#!/usr/bin/env python3
"""Freeze the complete GSE296875 well/barcode/donor/lineage join."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file


MAPPING_SHA256 = "5d3d946607f55b140ffa4a44369cf35f6e9e91f0f17276afc0e3e9f008554781"
LABELS_SHA256 = "45c077e9beb8607402bd2ed7bce0586596430d401f290ca86978f91e230ea67d"
FILELIST_SHA256 = "983a0764965f832d90ef1e008dcc6bf553450d35c508f4afebf3a1bab2b9768e"
EXPECTED_CELLS = 68_398
EXPECTED_DONORS = 39
EXPECTED_WELLS = tuple(f"well{index}" for index in range(1, 9))
EXPECTED_LABEL_CENSUS = {
    "B cells": (1_030, 38),
    "Cholangiocytes": (1_886, 39),
    "Hepatocytes": (46_286, 39),
    "Kupffer": (4_335, 39),
    "LSEC": (7_379, 39),
    "Mesenchymal": (4_781, 39),
    "NK-T": (2_701, 39),
}
LABEL_CONTRACT = {
    "B cells": ("b_cell", "secondary"),
    "Cholangiocytes": ("cholangiocyte", "primary"),
    "Hepatocytes": ("hepatocyte", "primary"),
    "Kupffer": ("macrophage", "primary"),
    "LSEC": ("endothelial_cell", "secondary"),
    "Mesenchymal": ("fibroblast", "primary"),
    "NK-T": ("t_cell", "primary"),
}
PRIMARY_LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
SECONDARY_LINEAGES = ("b_cell", "endothelial_cell")
SPLIT_SEED = 20260821
OUTER_FOLDS = 5


class MembershipError(ValueError):
    """Raised when the frozen source join does not meet its requirements."""


def fold_index(donor_id: str, *, seed: int = SPLIT_SEED) -> int:
    if not donor_id:
        raise MembershipError("donor identifier is empty")
    digest = sha256(f"{seed}\0{donor_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % OUTER_FOLDS


def read_exact_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise MembershipError(f"{path.name} columns differ from the contract")
        return [dict(row) for row in reader]


def parse_fragment_roster(filelist: Path) -> dict[str, dict[str, object]]:
    if sha256_file(filelist) != FILELIST_SHA256:
        raise MembershipError("GSE296875 filelist SHA-256 differs")
    with filelist.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "#Archive/File",
            "Name",
            "Time",
            "Size",
            "Type",
        ):
            raise MembershipError("GSE296875 filelist columns differ")
        rows = list(reader)
    roster: dict[str, dict[str, object]] = {}
    pattern = re.compile(
        r"^(GSM[0-9]+)_(well[1-8])_atac_fragments\.tsv\.gz$"
    )
    for row in rows:
        match = pattern.fullmatch(row["Name"])
        if row["#Archive/File"] != "File" or match is None:
            continue
        well = match.group(2)
        if well in roster:
            raise MembershipError(f"duplicate fragment file for {well}")
        roster[well] = {
            "well_id": well,
            "gsm": match.group(1),
            "fragment_filename": row["Name"],
            "fragment_size_bytes": int(row["Size"]),
        }
    if tuple(sorted(roster, key=lambda value: int(value[4:]))) != EXPECTED_WELLS:
        raise MembershipError("fragment file roster does not cover eight wells")
    if sum(int(row["fragment_size_bytes"]) for row in roster.values()) != 70_290_848_297:
        raise MembershipError("fragment byte total differs")
    return roster


def join_sources(mapping: Path, labels: Path) -> list[dict[str, str]]:
    if sha256_file(mapping) != MAPPING_SHA256:
        raise MembershipError("mapping metadata SHA-256 differs")
    if sha256_file(labels) != LABELS_SHA256:
        raise MembershipError("author labels SHA-256 differs")
    mapping_rows = read_exact_tsv(
        mapping, ("cell_id", "donor_id", "well_id", "raw_barcode")
    )
    label_rows = read_exact_tsv(labels, ("cell_id", "author_label"))
    if len(mapping_rows) != EXPECTED_CELLS or len(label_rows) != EXPECTED_CELLS:
        raise MembershipError("source tables do not contain 68,398 nuclei")
    by_cell = {row["cell_id"]: row for row in mapping_rows}
    if len(by_cell) != EXPECTED_CELLS:
        raise MembershipError("mapping cell identifiers are not unique")
    labels_by_cell = {row["cell_id"]: row["author_label"] for row in label_rows}
    if len(labels_by_cell) != EXPECTED_CELLS or set(labels_by_cell) != set(by_cell):
        raise MembershipError("mapping and author-label cell universes differ")

    joined: list[dict[str, str]] = []
    seen_well_barcodes: set[tuple[str, str]] = set()
    for cell_id, row in by_cell.items():
        well = row["well_id"]
        prefixed_barcode = row["raw_barcode"]
        prefix = f"{well}_"
        if well not in EXPECTED_WELLS or not prefixed_barcode.startswith(prefix):
            raise MembershipError("well-namespaced barcode encoding differs")
        raw_barcode = prefixed_barcode[len(prefix) :]
        if re.fullmatch(r"[ACGT]+-[0-9]+", raw_barcode) is None:
            raise MembershipError("raw barcode format differs")
        if cell_id != prefixed_barcode:
            raise MembershipError("cell ID and well-namespaced barcode differ")
        key = (well, raw_barcode)
        if key in seen_well_barcodes:
            raise MembershipError("well-local barcode is duplicated")
        seen_well_barcodes.add(key)
        source_label = labels_by_cell[cell_id]
        if source_label not in LABEL_CONTRACT:
            raise MembershipError(f"unexpected author label: {source_label}")
        lineage_id, analysis_role = LABEL_CONTRACT[source_label]
        joined.append(
            {
                "well_id": well,
                "raw_barcode": raw_barcode,
                "cell_id": cell_id,
                "donor_id": row["donor_id"],
                "source_label": source_label,
                "lineage_id": lineage_id,
                "analysis_role": analysis_role,
                "outer_fold": str(fold_index(row["donor_id"])),
            }
        )
    return joined


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with gzip.open(path, "xt", encoding="utf-8", newline="", compresslevel=6) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build(
    *,
    acquisition: Path,
    acquisition_artifacts_sha256: str,
    mapping: Path,
    labels: Path,
    output: Path,
    source_paths: Sequence[Path] = (),
) -> str:
    acquisition = acquisition.resolve(strict=True)
    if acquisition.is_symlink() or output.exists():
        raise MembershipError("unsafe acquisition or existing output")
    artifacts = acquisition / "ARTIFACTS.json"
    complete = acquisition / "COMPLETE"
    if not artifacts.is_file() or not complete.is_file():
        raise MembershipError("fragment acquisition is not complete")
    if sha256_file(artifacts) != acquisition_artifacts_sha256:
        raise MembershipError("fragment acquisition ARTIFACTS SHA-256 differs")
    manifest = json.loads(artifacts.read_text(encoding="utf-8"))
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("artifact_class") != "gse296875_raw_fragment_acquisition"
        or metadata.get("dataset_id") != "gse296875"
        or metadata.get("fragment_files") != 8
        or metadata.get("fragment_bytes") != 70_290_848_297
    ):
        raise MembershipError("fragment acquisition metadata differs")
    roster = parse_fragment_roster(acquisition / "filelist.txt")
    for well, row in roster.items():
        fragment = acquisition / "fragments" / str(row["fragment_filename"])
        index = Path(f"{fragment}.tbi")
        if (
            fragment.is_symlink()
            or not fragment.is_file()
            or fragment.stat().st_size != row["fragment_size_bytes"]
            or not index.is_file()
        ):
            raise MembershipError(f"fragment acquisition is incomplete for {well}")

    joined = join_sources(mapping.resolve(strict=True), labels.resolve(strict=True))
    donor_ids = sorted({row["donor_id"] for row in joined}, key=lambda value: int(value))
    if len(donor_ids) != EXPECTED_DONORS:
        raise MembershipError("donor census differs")

    label_census: dict[str, tuple[int, int]] = {}
    for label in LABEL_CONTRACT:
        rows = [row for row in joined if row["source_label"] == label]
        label_census[label] = (len(rows), len({row["donor_id"] for row in rows}))
    if label_census != EXPECTED_LABEL_CENSUS:
        raise MembershipError("author-label nucleus or donor census differs")

    output.mkdir(mode=0o750)
    barcode_dir = output / "barcodes"
    barcode_dir.mkdir(mode=0o750)
    fields = (
        "well_id",
        "raw_barcode",
        "cell_id",
        "donor_id",
        "source_label",
        "lineage_id",
        "analysis_role",
        "outer_fold",
    )
    joined.sort(key=lambda row: (int(row["well_id"][4:]), row["raw_barcode"]))
    write_tsv_gz(output / "cell_membership.tsv.gz", fields, joined)
    well_rows = []
    for well in EXPECTED_WELLS:
        rows = [row for row in joined if row["well_id"] == well]
        write_tsv_gz(barcode_dir / f"{well}.tsv.gz", fields, rows)
        well_rows.append(
            {
                **roster[well],
                "nuclei": len(rows),
                "donors": len({row["donor_id"] for row in rows}),
            }
        )
    write_tsv(
        output / "well_census.tsv",
        (
            "well_id",
            "gsm",
            "fragment_filename",
            "fragment_size_bytes",
            "nuclei",
            "donors",
        ),
        well_rows,
    )

    donor_lineage_rows = []
    for donor in donor_ids:
        for lineage in (*PRIMARY_LINEAGES, *SECONDARY_LINEAGES):
            rows = [
                row
                for row in joined
                if row["donor_id"] == donor and row["lineage_id"] == lineage
            ]
            donor_lineage_rows.append(
                {
                    "donor_id": donor,
                    "outer_fold": fold_index(donor),
                    "lineage_id": lineage,
                    "analysis_role": (
                        "primary" if lineage in PRIMARY_LINEAGES else "secondary"
                    ),
                    "nuclei": len(rows),
                    "wells": ",".join(sorted({row["well_id"] for row in rows})),
                    "observed": str(bool(rows)).lower(),
                }
            )
    write_tsv(
        output / "donor_lineage_census.tsv",
        (
            "donor_id",
            "outer_fold",
            "lineage_id",
            "analysis_role",
            "nuclei",
            "wells",
            "observed",
        ),
        donor_lineage_rows,
    )
    write_tsv(
        output / "donor_folds.tsv",
        ("donor_id", "outer_fold"),
        ({"donor_id": donor, "outer_fold": fold_index(donor)} for donor in donor_ids),
    )
    contract = {
        "schema_version": "masld-bench-gse296875-fragment-membership-v1",
        "dataset_id": "gse296875",
        "pairing": "same_nucleus",
        "biological_unit": "donor",
        "technical_batch": "well",
        "cells": EXPECTED_CELLS,
        "donors": EXPECTED_DONORS,
        "wells": list(EXPECTED_WELLS),
        "primary_lineages": list(PRIMARY_LINEAGES),
        "secondary_lineages": list(SECONDARY_LINEAGES),
        "primary_endpoint_unchanged_from_smoke": True,
        "secondary_lineages_are_not_primary_gate_members": True,
        "split_id": "donor_outer",
        "split_seed": SPLIT_SEED,
        "outer_folds": OUTER_FOLDS,
        "missing_barcode_encoding": "absent_from_filtered_author_cell_roster",
        "acquisition_artifacts_sha256": acquisition_artifacts_sha256,
        "mapping_sha256": MAPPING_SHA256,
        "author_labels_sha256": LABELS_SHA256,
        "filelist_sha256": FILELIST_SHA256,
    }
    (output / "contract.json").write_text(
        json.dumps(contract, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    if source_paths:
        source_rows = []
        for source_path in source_paths:
            resolved = source_path.resolve(strict=True)
            if not resolved.is_file() or resolved.is_symlink():
                raise MembershipError("source provenance path is unsafe")
            source_rows.append(f"{sha256_file(resolved)}  {resolved}\n")
        (output / "source.sha256").write_text(
            "".join(source_rows), encoding="utf-8"
        )
    return freeze_tree(
        output,
        {
            "artifact_class": "gse296875_fragment_membership",
            "dataset_id": "gse296875",
            "cells": EXPECTED_CELLS,
            "donors": EXPECTED_DONORS,
            "primary_lineages": len(PRIMARY_LINEAGES),
            "secondary_lineages": len(SECONDARY_LINEAGES),
            "acquisition_artifacts_sha256": acquisition_artifacts_sha256,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--acquisition-artifacts-sha256", required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-path", action="append", type=Path, default=[])
    arguments = parser.parse_args()
    digest = build(
        acquisition=arguments.acquisition,
        acquisition_artifacts_sha256=arguments.acquisition_artifacts_sha256,
        mapping=arguments.mapping,
        labels=arguments.labels,
        output=arguments.output,
        source_paths=arguments.source_path,
    )
    print(arguments.output.resolve(strict=True))
    print(digest)


if __name__ == "__main__":
    main()
