#!/usr/bin/env python3
"""Freeze outcome-blind GSE256398 TranscriptFormer activation memberships."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import h5py

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


EXPECTED_SOURCE_ARTIFACTS = (
    "f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20"
)
EXPECTED_SOURCE_EXACT_ARTIFACTS = (
    "3499eed58c6da15e644f5114da33b2874e3252a86fad6d6937960424e781ecbc"
)
EXPECTED_PROSPECTIVE_ARTIFACTS = (
    "77b7b62b30c49fdf306039a3c31930040a5a140d23a64dd24875370cc11e91ad"
)
EXPECTED_SOURCE_EXACT_ROWS = 165_372
EXPECTED_PROSPECTIVE_ROWS = 172_997
EXPECTED_DONORS = 26
H5_PATTERN = re.compile(
    r"^(GSM809\d{4})_(S\d+)_CB_raw_feature_bc_matrix_filtered\.h5$"
)
FIELDS = ("row_id", "gsm", "source_sample_id", "barcode")


class MembershipLockError(RuntimeError):
    """Raised when a frozen input or activation membership differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def donor_sort_key(value: str) -> tuple[int, str]:
    if not re.fullmatch(r"S\d+", value):
        raise MembershipLockError(f"unexpected source_sample_id: {value}")
    return int(value[1:]), value


def load_membership(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise MembershipLockError(f"membership schema differs: {path}")
        rows = list(reader)
    if not rows or len({row["row_id"] for row in rows}) != len(rows):
        raise MembershipLockError(f"membership is empty or duplicates row IDs: {path}")
    for row in rows:
        expected = f"gse256398:{row['source_sample_id']}:{row['barcode']}"
        if row["row_id"] != expected:
            raise MembershipLockError(f"row identity differs: {row['row_id']}")
    return rows


def canonical(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (donor_sort_key(row["source_sample_id"]), row["barcode"]),
    )


def write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def source_barcode_roster(source: Path) -> dict[str, tuple[str, set[str]]]:
    roster: dict[str, tuple[str, set[str]]] = {}
    for path in sorted((source / "h5").glob("*.h5")):
        match = H5_PATTERN.fullmatch(path.name)
        if match is None:
            raise MembershipLockError(f"source H5 filename differs: {path.name}")
        gsm, donor = match.groups()
        if donor in roster:
            raise MembershipLockError(f"duplicate source H5 donor: {donor}")
        with h5py.File(path, "r") as handle:
            values = handle["matrix/barcodes"][:]
        barcodes = {
            value.decode("utf-8") if isinstance(value, bytes) else str(value)
            for value in values
        }
        if len(barcodes) != len(values):
            raise MembershipLockError(f"duplicate source H5 barcode: {donor}")
        roster[donor] = (gsm, barcodes)
    return roster


def validate_against_source(
    rows: list[dict[str, str]], roster: dict[str, tuple[str, set[str]]], label: str
) -> None:
    donors = {row["source_sample_id"] for row in rows}
    if donors != set(roster):
        raise MembershipLockError(f"{label} donor roster differs from source H5 roster")
    for row in rows:
        gsm, barcodes = roster[row["source_sample_id"]]
        if row["gsm"] != gsm or row["barcode"] not in barcodes:
            raise MembershipLockError(f"{label} row is absent from source H5: {row['row_id']}")


def counts_by_donor(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        donor = row["source_sample_id"]
        counts[donor] = counts.get(donor, 0) + 1
    return counts


def build(
    *, source: Path, source_exact_root: Path, prospective_root: Path, output: Path
) -> dict[str, Any]:
    if output.exists():
        raise MembershipLockError("membership lock output already exists")
    expected = (
        (source, EXPECTED_SOURCE_ARTIFACTS),
        (source_exact_root, EXPECTED_SOURCE_EXACT_ARTIFACTS),
        (prospective_root, EXPECTED_PROSPECTIVE_ARTIFACTS),
    )
    for root, digest in expected:
        verify_frozen_tree(root)
        if sha256_file(root / "ARTIFACTS.json") != digest:
            raise MembershipLockError(f"frozen input differs: {root}")

    source_exact = load_membership(source_exact_root / "retained_membership.tsv")
    prospective = load_membership(
        prospective_root / "expected_rate_capped_membership.tsv"
    )
    if len(source_exact) != EXPECTED_SOURCE_EXACT_ROWS:
        raise MembershipLockError("source-exact row count differs")
    if len(prospective) != EXPECTED_PROSPECTIVE_ROWS:
        raise MembershipLockError("prospective-capped row count differs")

    source_by_id = {row["row_id"]: row for row in source_exact}
    prospective_by_id = {row["row_id"]: row for row in prospective}
    if not set(source_by_id).issubset(prospective_by_id):
        raise MembershipLockError("source-exact membership is not a prospective-capped subset")
    for row_id, row in source_by_id.items():
        if row != prospective_by_id[row_id]:
            raise MembershipLockError(f"membership identity fields differ: {row_id}")

    roster = source_barcode_roster(source)
    if len(roster) != EXPECTED_DONORS:
        raise MembershipLockError("source H5 donor count differs")
    validate_against_source(source_exact, roster, "source_exact")
    validate_against_source(prospective, roster, "prospective_capped")

    source_exact = canonical(source_exact)
    prospective = canonical(prospective)
    source_ids = set(source_by_id)
    union = [
        {
            **row,
            "in_source_exact": "true" if row["row_id"] in source_ids else "false",
            "in_prospective_capped": "true",
        }
        for row in prospective
    ]
    source_without_s35 = [row for row in source_exact if row["source_sample_id"] != "S35"]
    prospective_without_s35 = [
        row for row in prospective if row["source_sample_id"] != "S35"
    ]

    output.mkdir(parents=True, exist_ok=False)
    files = {
        "source_exact": output / "source_exact_membership.tsv",
        "prospective_capped": output / "prospective_capped_membership.tsv",
        "union": output / "prospective_union_membership.tsv",
        "source_exact_leave_s35_out": output / "source_exact_leave_s35_out.tsv",
        "prospective_capped_leave_s35_out": output / "prospective_capped_leave_s35_out.tsv",
    }
    write_tsv(files["source_exact"], FIELDS, source_exact)
    write_tsv(files["prospective_capped"], FIELDS, prospective)
    write_tsv(
        files["union"],
        FIELDS + ("in_source_exact", "in_prospective_capped"),
        union,
    )
    write_tsv(files["source_exact_leave_s35_out"], FIELDS, source_without_s35)
    write_tsv(
        files["prospective_capped_leave_s35_out"],
        FIELDS,
        prospective_without_s35,
    )

    exact_counts = counts_by_donor(source_exact)
    prospective_counts = counts_by_donor(prospective)
    donor_rows = [
        {
            "source_sample_id": donor,
            "source_exact_rows": exact_counts[donor],
            "prospective_capped_rows": prospective_counts[donor],
            "prospective_only_rows": prospective_counts[donor] - exact_counts[donor],
        }
        for donor in sorted(roster, key=donor_sort_key)
    ]
    write_tsv(
        output / "donor_counts.tsv",
        (
            "source_sample_id",
            "source_exact_rows",
            "prospective_capped_rows",
            "prospective_only_rows",
        ),
        donor_rows,
    )

    receipt = {
        "schema_version": "masld-bench-gse256398-membership-lock-v1",
        "status": "pass_outcome_blind_membership_lock",
        "dataset_id": "gse256398",
        "arms": {
            "source_exact": len(source_exact),
            "prospective_capped": len(prospective),
            "source_exact_leave_s35_out": len(source_without_s35),
            "prospective_capped_leave_s35_out": len(prospective_without_s35),
        },
        "donors": len(roster),
        "source_exact_is_prospective_subset": True,
        "prospective_only_rows": len(prospective) - len(source_exact),
        "membership_files": {
            label: {
                "path": path.name,
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for label, path in files.items()
        },
        "donor_counts_sha256": sha256_file(output / "donor_counts.tsv"),
        "phenotype_metadata_read": False,
        "barcode_cell_labels_read": False,
        "model_predictions_read": False,
        "sealed_outcomes_read": False,
        "supervised_accuracy_claim_allowed": False,
        "next_step": "materialize_prospective_union_raw_counts_and_transcriptformer_tokens",
    }
    write_json_exclusive(output / "membership_lock_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse256398_outcome_blind_transcriptformer_membership_lock",
            "dataset_id": "gse256398",
            "donors": len(roster),
            "source_exact_rows": len(source_exact),
            "prospective_capped_rows": len(prospective),
            "phenotype_metadata_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-exact-root", type=Path, required=True)
    parser.add_argument("--prospective-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = build(
        source=args.source,
        source_exact_root=args.source_exact_root,
        prospective_root=args.prospective_root,
        output=args.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
