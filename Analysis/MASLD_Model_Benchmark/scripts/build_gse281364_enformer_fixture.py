#!/usr/bin/env python3
"""Build outcome-blind 196,608-bp Enformer windows for the MPRA locus screen."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from scripts.alphagenome_sei_build_fixture import IndexedFasta


WINDOW = 196_608
CENTER = WINDOW // 2
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "variant_pos1",
    "input_start0",
    "input_end0",
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reference_reverse_complement_sha256",
    "alternative_reverse_complement_sha256",
)


class EnformerFixtureError(ValueError):
    """Raised when an outcome-blind Enformer fixture invariant differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def allele_window(
    reference: IndexedFasta, row: dict[str, str]
) -> tuple[str, str, int, int]:
    position = int(row["variant_pos0"])
    start, end = position - CENTER, position - CENTER + WINDOW
    sequence = reference.fetch(row["contig"], start, end)
    if (
        len(sequence) != WINDOW
        or set(sequence) - set("ACGTN")
        or sequence[CENTER] != row["genomic_ref"]
        or row["genomic_ref"] == row["genomic_alt"]
        or row["genomic_alt"] not in "ACGT"
    ):
        raise EnformerFixtureError("reference window or allele differs")
    alternative = sequence[:CENTER] + row["genomic_alt"] + sequence[CENTER + 1 :]
    if sum(left != right for left, right in zip(sequence, alternative)) != 1:
        raise EnformerFixtureError("alternative window differs at more than one base")
    return sequence, alternative, start, end


def build(
    *,
    split_root: Path,
    split_sha256: str,
    sei_fixture_root: Path,
    sei_fixture_sha256: str,
    fasta_root: Path,
    fasta_sha256: str,
    output: Path,
) -> dict[str, object]:
    if output.exists():
        raise EnformerFixtureError("output exists")
    authorities = (
        (split_root, split_sha256),
        (sei_fixture_root, sei_fixture_sha256),
        (fasta_root, fasta_sha256),
    )
    if any(digest(root / "ARTIFACTS.json") != expected for root, expected in authorities):
        raise EnformerFixtureError("upstream artifact changed")
    split_receipt = json.loads((split_root / "split/receipt.json").read_text())
    sei_receipt = json.loads((sei_fixture_root / "fixture/receipt.json").read_text())
    fasta_artifacts = json.loads((fasta_root / "ARTIFACTS.json").read_text())
    if (
        split_receipt.get("status") != "pass_outcome_blind_split_contract"
        or split_receipt.get("outcomes_read") is not False
        or sei_receipt.get("status") != "pass_outcome_blind_sei_fixture"
        or sei_receipt.get("outcomes_read") is not False
        or fasta_artifacts.get("metadata", {}).get("build") != "GRCh38.p14"
        or fasta_artifacts.get("metadata", {}).get("indexed") is not True
    ):
        raise EnformerFixtureError("upstream outcome firewall or reference differs")

    split_rows = {row["element_id"]: row for row in read_tsv(split_root / "split/elements.tsv")}
    sei_rows = read_tsv(sei_fixture_root / "fixture/manifest.tsv")
    if (
        len(split_rows) != 4_359
        or len(sei_rows) != 1_033
        or len({row["outer_locus_sequence_group_id"] for row in sei_rows}) != 1_033
    ):
        raise EnformerFixtureError("source census differs")
    selected: list[dict[str, str]] = []
    for sei_row in sei_rows:
        row = split_rows.get(sei_row["element_id"])
        if row is None or any(
            row[field] != sei_row[field]
            for field in ("outer_locus_sequence_group_id", "outer_fold", "contig", "variant_pos0")
        ):
            raise EnformerFixtureError("Sei common-locus selection differs")
        selected.append(row)

    reference = IndexedFasta(
        fasta_root / "GRCh38.p14.sequence_model.fa",
        fasta_root / "GRCh38.p14.sequence_model.fa.fai",
    )
    manifest_rows: list[dict[str, object]] = []
    for row in selected:
        ref_sequence, alt_sequence, start, end = allele_window(reference, row)
        ref_rc, alt_rc = reverse_complement(ref_sequence), reverse_complement(alt_sequence)
        manifest_rows.append(
            {
                "fixture_id": "enformer_" + sha256(row["element_id"].encode()).hexdigest()[:24],
                "element_id": row["element_id"],
                "outer_locus_sequence_group_id": row["outer_locus_sequence_group_id"],
                "outer_fold": row["outer_fold"],
                "contig": row["contig"],
                "variant_pos0": row["variant_pos0"],
                "variant_pos1": row["variant_pos1"],
                "input_start0": start,
                "input_end0": end,
                "variant_index0": CENTER,
                "ref": row["genomic_ref"],
                "alt": row["genomic_alt"],
                "reference_sequence_sha256": sha256(ref_sequence.encode()).hexdigest(),
                "alternative_sequence_sha256": sha256(alt_sequence.encode()).hexdigest(),
                "reference_reverse_complement_sha256": sha256(ref_rc.encode()).hexdigest(),
                "alternative_reverse_complement_sha256": sha256(alt_rc.encode()).hexdigest(),
            }
        )

    output.mkdir(parents=True)
    with (output / "manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-enformer-fixture-v1",
        "status": "pass_outcome_blind_enformer_fixture",
        "dataset_id": "gse281364",
        "model_id": "enformer_crested_restricted_port",
        "elements": 1_033,
        "outer_locus_sequence_groups": 1_033,
        "outer_folds": 5,
        "input_length_bp": WINDOW,
        "variant_index0": CENTER,
        "selection": "exact_same_elements_as_frozen_sei_one_per_locus_group_fixture",
        "split_artifacts_sha256": split_sha256,
        "sei_fixture_artifacts_sha256": sei_fixture_sha256,
        "fasta_artifacts_sha256": fasta_sha256,
        "reference_build": "GRCh38.p14",
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "champion_eligible": False,
        "restriction": "internal_restricted_comparator_only",
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--split-sha256", required=True)
    parser.add_argument("--sei-fixture-root", type=Path, required=True)
    parser.add_argument("--sei-fixture-sha256", required=True)
    parser.add_argument("--fasta-root", type=Path, required=True)
    parser.add_argument("--fasta-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
