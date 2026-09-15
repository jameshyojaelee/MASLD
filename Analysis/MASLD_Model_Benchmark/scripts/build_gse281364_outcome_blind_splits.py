#!/usr/bin/env python3
"""Freeze outcome-blind merged-locus folds for the GSE281364 MPRA substrate."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Iterable, Mapping


EXPECTED_ELEMENTS_SHA256 = "990d12227dad16781ead0b87abb465d48969170b59298b7dd09da61af39d0dd8"
EXPECTED_QUALIFIED = 4_359
EXPECTED_GROUPS = 1_033
HEX64 = re.compile(r"^[0-9a-f]{64}$")
GROUP_ID = re.compile(r"^outer_[0-9a-f]{20}$")
BASES = frozenset("ACGT")
OUTPUT_FIELDS = (
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "reference_match_start0",
    "reference_match_end0",
    "variant_pos0",
    "variant_pos1",
    "genomic_ref",
    "genomic_alt",
    "oligo_difference_index0",
    "ref_sequence_107bp",
    "alt_sequence_107bp",
    "ref_sequence_sha256",
    "alt_sequence_sha256",
    "canonical_ref_sequence_sha256",
    "canonical_alt_sequence_sha256",
)


class SplitContractError(ValueError):
    """Raised when the outcome-blind split substrate differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SplitContractError("elements table has no header")
        return [dict(row) for row in reader]


def write_tsv(path: Path, fields: Iterable[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_source(elements: Path, artifacts: Path, expected_artifacts_sha256: str) -> None:
    if digest(artifacts) != expected_artifacts_sha256:
        raise SplitContractError("source artifact manifest changed")
    manifest = json.loads(artifacts.read_text(encoding="utf-8"))
    records = {row["path"]: row for row in manifest.get("artifacts", [])}
    record = records.get("validated/elements.tsv")
    if (
        manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or manifest.get("metadata", {}).get("qualified_paired_snv_elements") != EXPECTED_QUALIFIED
        or manifest.get("metadata", {}).get("outcome_role") != "exposed_development_MPRA_only"
        or record is None
        or record.get("sha256") != EXPECTED_ELEMENTS_SHA256
        or digest(elements) != EXPECTED_ELEMENTS_SHA256
    ):
        raise SplitContractError("validated element authority differs")


def qualified_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = [row for row in rows if row.get("pair_state") == "paired_snv"]
    if len(selected) != EXPECTED_QUALIFIED:
        raise SplitContractError("qualified paired-SNV count differs")
    identifiers: set[str] = set()
    sequence_owners: dict[str, str] = {}
    by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        identifier = row["element_id"]
        group = row["outer_locus_sequence_group_id"]
        reference, alternative = row["ref_sequence_107bp"], row["alt_sequence_107bp"]
        difference = int(row["oligo_difference_index0"])
        if (
            identifier in identifiers
            or not GROUP_ID.fullmatch(group)
            or row["sequence_state"] != "observed_reconstructed_from_raw_plasmid_reads"
            or row["sequence_qc_reason"] != "not_applicable"
            or len(reference) != 107
            or len(alternative) != 107
            or set(reference) - BASES
            or set(alternative) - BASES
            or not 0 <= difference < 107
            or sum(left != right for left, right in zip(reference, alternative)) != 1
            or reference[difference] != row["genomic_ref"]
            or alternative[difference] != row["genomic_alt"]
            or row["genomic_ref"] == row["genomic_alt"]
        ):
            raise SplitContractError("qualified element contract differs")
        identifiers.add(identifier)
        for field in (
            "ref_sequence_sha256",
            "alt_sequence_sha256",
            "canonical_ref_sequence_sha256",
            "canonical_alt_sequence_sha256",
        ):
            if not HEX64.fullmatch(row[field]):
                raise SplitContractError("sequence hash differs")
        for field in ("canonical_ref_sequence_sha256", "canonical_alt_sequence_sha256"):
            owner = sequence_owners.setdefault(row[field], group)
            if owner != group:
                raise SplitContractError("exact or reverse-complement sequence identity spans groups")
        by_contig[row["contig"]].append(row)
    for rows_on_contig in by_contig.values():
        ordered = sorted(rows_on_contig, key=lambda row: int(row["variant_pos0"]))
        for left, right in zip(ordered, ordered[1:]):
            if (
                int(right["variant_pos0"]) - int(left["variant_pos0"]) < 6_000
                and left["outer_locus_sequence_group_id"]
                != right["outer_locus_sequence_group_id"]
            ):
                raise SplitContractError("nearby loci span outer groups")
    return selected


def assign_groups(rows: list[dict[str, str]], folds: int, seed: str) -> dict[str, int]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["outer_locus_sequence_group_id"]].append(row)
    if len(grouped) != EXPECTED_GROUPS or folds != 5:
        raise SplitContractError("outer group or fold count differs")
    summaries: list[tuple[str, str, int]] = []
    for group, members in grouped.items():
        contigs = {row["contig"] for row in members}
        if len(contigs) != 1:
            raise SplitContractError("outer group spans contigs")
        summaries.append((group, next(iter(contigs)), len(members)))
    summaries.sort(
        key=lambda item: (
            -item[2],
            sha256(f"{seed}\0{item[0]}".encode()).hexdigest(),
        )
    )
    totals = [0] * folds
    contig_totals: dict[str, list[int]] = defaultdict(lambda: [0] * folds)
    assignments: dict[str, int] = {}
    for group, contig, size in summaries:
        fold = min(
            range(folds),
            key=lambda candidate: (
                contig_totals[contig][candidate],
                totals[candidate],
                sha256(f"{seed}\0{group}\0{candidate}".encode()).hexdigest(),
            ),
        )
        assignments[group] = fold
        totals[fold] += size
        contig_totals[contig][fold] += size
    if set(assignments.values()) != set(range(folds)):
        raise SplitContractError("an outer fold is empty")
    return assignments


def build(
    elements: Path,
    artifacts: Path,
    expected_artifacts_sha256: str,
    output: Path,
    seed: str,
) -> dict[str, object]:
    if output.exists():
        raise SplitContractError("output exists")
    validate_source(elements, artifacts, expected_artifacts_sha256)
    rows = qualified_rows(read_tsv(elements))
    assignments = assign_groups(rows, 5, seed)
    output.mkdir(parents=True)
    element_rows = [
        {
            field: assignments[row["outer_locus_sequence_group_id"]]
            if field == "outer_fold"
            else row[field]
            for field in OUTPUT_FIELDS
        }
        for row in sorted(rows, key=lambda row: row["element_id"])
    ]
    write_tsv(output / "elements.tsv", OUTPUT_FIELDS, element_rows)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["outer_locus_sequence_group_id"]].append(row)
    group_rows = []
    for group in sorted(grouped):
        members = grouped[group]
        group_rows.append(
            {
                "outer_locus_sequence_group_id": group,
                "outer_fold": assignments[group],
                "contig": members[0]["contig"],
                "elements": len(members),
                "minimum_variant_pos0": min(int(row["variant_pos0"]) for row in members),
                "maximum_variant_pos0": max(int(row["variant_pos0"]) for row in members),
            }
        )
    write_tsv(
        output / "groups.tsv",
        (
            "outer_locus_sequence_group_id",
            "outer_fold",
            "contig",
            "elements",
            "minimum_variant_pos0",
            "maximum_variant_pos0",
        ),
        group_rows,
    )
    fold_elements = {
        str(fold): sum(int(row["outer_fold"]) == fold for row in element_rows)
        for fold in range(5)
    }
    fold_groups = {
        str(fold): sum(int(row["outer_fold"]) == fold for row in group_rows)
        for fold in range(5)
    }
    assignment_hash = sha256(
        "\n".join(f"{group}\t{assignments[group]}" for group in sorted(assignments)).encode()
    ).hexdigest()
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-outcome-blind-splits-v1",
        "status": "pass_outcome_blind_split_contract",
        "dataset_id": "gse281364",
        "qualified_paired_snv_elements": len(rows),
        "outer_locus_sequence_groups": len(assignments),
        "outer_folds": 5,
        "fold_elements": fold_elements,
        "fold_groups": fold_groups,
        "split_seed": seed,
        "assignment_sha256": assignment_hash,
        "grouping_contract": "variant_distance_lt_6000_or_exact_or_reverse_complement_sequence_identity",
        "replication_unit": "experimental_replicate_not_donor",
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "champion_eligible": False,
        "source_artifacts_sha256": expected_artifacts_sha256,
        "elements_sha256": EXPECTED_ELEMENTS_SHA256,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--elements", type=Path, required=True)
    parser.add_argument("--source-artifacts", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", default="gse281364-mpra-outer-v1")
    arguments = parser.parse_args()
    print(
        json.dumps(
            build(
                arguments.elements,
                arguments.source_artifacts,
                arguments.source_artifacts_sha256,
                arguments.output,
                arguments.seed,
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
