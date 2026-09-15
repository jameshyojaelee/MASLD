#!/usr/bin/env python3
"""Independently validate Enformer MPRA windows, alleles, grouping, and separation."""

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


class FixtureValidationError(ValueError):
    """Raised when the frozen Enformer fixture fails independent validation."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def validate(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists():
        raise FixtureValidationError("output exists")
    expected_artifacts = (
        (arguments.reconstruction_root, arguments.reconstruction_sha256),
        (arguments.split_root, arguments.split_sha256),
        (arguments.sei_fixture_root, arguments.sei_fixture_sha256),
        (arguments.fasta_root, arguments.fasta_sha256),
        (arguments.fixture_root, arguments.fixture_sha256),
    )
    if any(digest(root / "ARTIFACTS.json") != expected for root, expected in expected_artifacts):
        raise FixtureValidationError("upstream artifact hash differs")
    fixture_files = sorted(
        path.relative_to(arguments.fixture_root).as_posix()
        for path in arguments.fixture_root.rglob("*")
        if path.is_file() and path.name not in {"ARTIFACTS.json", "COMPLETE"}
    )
    if any(
        token in path.lower()
        for path in fixture_files
        for token in ("outcome", "reporter", "activity", "label")
    ):
        raise FixtureValidationError("fixture contains a forbidden outcome-named file")

    reconstruction = {
        row["element_id"]: row
        for row in read_tsv(arguments.reconstruction_root / "validated/elements.tsv")
        if row["pair_state"] == "paired_snv"
    }
    split = {row["element_id"]: row for row in read_tsv(arguments.split_root / "split/elements.tsv")}
    sei = {row["element_id"]: row for row in read_tsv(arguments.sei_fixture_root / "fixture/manifest.tsv")}
    fixture = read_tsv(arguments.fixture_root / "fixture/manifest.tsv")
    receipts = (
        json.loads((arguments.reconstruction_root / "validated/reconstruction_receipt.json").read_text()),
        json.loads((arguments.split_root / "split/receipt.json").read_text()),
        json.loads((arguments.sei_fixture_root / "fixture/receipt.json").read_text()),
        json.loads((arguments.fixture_root / "fixture/receipt.json").read_text()),
    )
    reconstruction_receipt, split_receipt, sei_receipt, fixture_receipt = receipts
    if (
        len(reconstruction) != 4_359
        or len(split) != 4_359
        or len(sei) != 1_033
        or len(fixture) != 1_033
        or reconstruction_receipt.get("sealed_outcomes_loaded") is not False
        or any(receipt.get("sealed_outcomes_read") is not False for receipt in receipts[1:])
        or any(receipt.get("outcomes_read") is not False for receipt in receipts[1:])
    ):
        raise FixtureValidationError("census or outcome firewall differs")

    reference = IndexedFasta(
        arguments.fasta_root / "GRCh38.p14.sequence_model.fa",
        arguments.fasta_root / "GRCh38.p14.sequence_model.fa.fai",
    )
    fold_counts = {str(fold): 0 for fold in range(5)}
    observed_groups: set[str] = set()
    for row in fixture:
        element = row["element_id"]
        original, split_row, sei_row = reconstruction.get(element), split.get(element), sei.get(element)
        if original is None or split_row is None or sei_row is None:
            raise FixtureValidationError("fixture element is absent from a source authority")
        exact_fields = (
            "outer_locus_sequence_group_id", "outer_fold", "contig", "variant_pos0",
            "variant_pos1", "genomic_ref", "genomic_alt",
        )
        if any(split_row[field] != sei_row[field] for field in exact_fields[:4]):
            raise FixtureValidationError("Sei shared-element grouping differs")
        if any(
            row[target] != split_row[source]
            for target, source in (
                ("outer_locus_sequence_group_id", "outer_locus_sequence_group_id"),
                ("outer_fold", "outer_fold"),
                ("contig", "contig"),
                ("variant_pos0", "variant_pos0"),
                ("variant_pos1", "variant_pos1"),
                ("ref", "genomic_ref"),
                ("alt", "genomic_alt"),
            )
        ):
            raise FixtureValidationError("Enformer-to-split allele or grouping join differs")
        if (
            original["reference_match_orientation"] != "forward"
            or original["sequence_state"] != "observed_reconstructed_from_raw_plasmid_reads"
            or original["sequence_qc_reason"] != "not_applicable"
        ):
            raise FixtureValidationError("construct orientation or QC state differs")
        start, end = int(row["input_start0"]), int(row["input_end0"])
        position = int(row["variant_pos0"])
        if end - start != WINDOW or position - start != CENTER or int(row["variant_index0"]) != CENTER:
            raise FixtureValidationError("window anchoring differs")
        sequence = reference.fetch(row["contig"], start, end)
        if sequence[CENTER] != row["ref"] or sha256(sequence.encode()).hexdigest() != row["reference_sequence_sha256"]:
            raise FixtureValidationError("reference sequence or REF allele differs")
        alternative = sequence[:CENTER] + row["alt"] + sequence[CENTER + 1 :]
        if (
            sum(left != right for left, right in zip(sequence, alternative)) != 1
            or sha256(alternative.encode()).hexdigest() != row["alternative_sequence_sha256"]
        ):
            raise FixtureValidationError("ALT substitution differs")
        construct_start = int(original["reference_match_start0"]) - start
        construct_end = int(original["reference_match_end0"]) - start
        if (
            sequence[construct_start:construct_end] != original["ref_sequence_107bp"]
            or alternative[construct_start:construct_end] != original["alt_sequence_107bp"]
            or construct_start + int(original["oligo_difference_index0"]) != CENTER
        ):
            raise FixtureValidationError("genomic window and reconstructed 107-bp construct differ")
        ref_rc, alt_rc = reverse_complement(sequence), reverse_complement(alternative)
        if (
            reverse_complement(ref_rc) != sequence
            or reverse_complement(alt_rc) != alternative
            or sha256(ref_rc.encode()).hexdigest() != row["reference_reverse_complement_sha256"]
            or sha256(alt_rc.encode()).hexdigest() != row["alternative_reverse_complement_sha256"]
        ):
            raise FixtureValidationError("reverse-complement identity or grouping differs")
        fold_counts[row["outer_fold"]] += 1
        observed_groups.add(row["outer_locus_sequence_group_id"])

    if len(observed_groups) != 1_033 or set(fold_counts) != {"0", "1", "2", "3", "4"}:
        raise FixtureValidationError("outer group or fold coverage differs")
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-enformer-fixture-validation-v1",
        "status": "pass_independent_full_census",
        "dataset_id": "gse281364",
        "model_id": "enformer_crested_restricted_port",
        "elements": 1_033,
        "outer_locus_sequence_groups": 1_033,
        "outer_fold_counts": fold_counts,
        "forward_construct_orientations": 1_033,
        "ref_alt_single_substitutions": 1_033,
        "reverse_complement_involutions": 2_066,
        "construct_107bp_exact_matches": 2_066,
        "fixture_file_inventory": fixture_files,
        "fixture_artifacts_sha256": arguments.fixture_sha256,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "champion_eligible": False,
    }
    arguments.output.mkdir(parents=True)
    (arguments.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("reconstruction", "split", "sei_fixture", "fasta", "fixture"):
        parser.add_argument(f"--{name.replace('_', '-')}-root", type=Path, required=True)
        parser.add_argument(f"--{name.replace('_', '-')}-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(validate(parser.parse_args()), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
