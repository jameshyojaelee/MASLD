#!/usr/bin/env python3
"""Build outcome-blind 230-bp GSE281364 inputs for MPRALegNet HepG2."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Iterable

from scripts.alphagenome_sei_build_fixture import IndexedFasta


WINDOW = 230
FORWARD_VARIANT_INDEX0 = 115
REVERSE_VARIANT_INDEX0 = WINDOW - 1 - FORWARD_VARIANT_INDEX0
EXPECTED_ELEMENTS = 4_359
EXPECTED_GROUPS = 1_033
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MANIFEST_FIELDS = (
    "fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "variant_pos1",
    "input_start0",
    "input_end0",
    "forward_variant_index0",
    "reverse_complement_variant_index0",
    "ref",
    "alt",
    "cell_context",
    "input_lane",
    "source_assay_insert_length_bp",
    "source_oligo_difference_index0",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reference_reverse_complement_sha256",
    "alternative_reverse_complement_sha256",
)


class MPRALegNetFixtureError(ValueError):
    """Raised when the outcome-blind MPRALegNet input requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def allele_window(reference: IndexedFasta, row: dict[str, str]) -> tuple[str, str, int, int]:
    position = int(row["variant_pos0"])
    start, end = position - FORWARD_VARIANT_INDEX0, position - FORWARD_VARIANT_INDEX0 + WINDOW
    sequence = reference.fetch(row["contig"], start, end)
    if (
        len(sequence) != WINDOW
        or set(sequence) - set("ACGT")
        or sequence[FORWARD_VARIANT_INDEX0] != row["genomic_ref"]
        or row["genomic_ref"] == row["genomic_alt"]
        or row["genomic_alt"] not in "ACGT"
        or len(row["ref_sequence_107bp"]) != 107
        or len(row["alt_sequence_107bp"]) != 107
        or int(row["oligo_difference_index0"]) != 62
    ):
        raise MPRALegNetFixtureError("reference, allele, or source oligo contract differs")
    alternative = (
        sequence[:FORWARD_VARIANT_INDEX0]
        + row["genomic_alt"]
        + sequence[FORWARD_VARIANT_INDEX0 + 1 :]
    )
    if sum(left != right for left, right in zip(sequence, alternative, strict=True)) != 1:
        raise MPRALegNetFixtureError("alternative window differs at more than one base")
    return sequence, alternative, start, end


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as handle:
                for name, sequence in records:
                    handle.write(f">{name}\n{sequence}\n")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    groups = {row["outer_locus_sequence_group_id"] for row in rows}
    if (
        len(rows) != EXPECTED_ELEMENTS
        or len(groups) != EXPECTED_GROUPS
        or {int(row["outer_fold"]) for row in rows} != set(range(5))
        or len({row["element_id"] for row in rows}) != EXPECTED_ELEMENTS
    ):
        raise MPRALegNetFixtureError("split element, group, or fold census differs")
    return rows


def build(
    *,
    split_root: Path,
    split_sha256: str,
    reference_root: Path,
    reference_sha256: str,
    output: Path,
) -> dict[str, object]:
    if output.exists():
        raise MPRALegNetFixtureError("output exists")
    if digest(split_root / "ARTIFACTS.json") != split_sha256:
        raise MPRALegNetFixtureError("split artifact changed")
    if digest(reference_root / "ARTIFACTS.json") != reference_sha256:
        raise MPRALegNetFixtureError("reference artifact changed")
    split_receipt = json.loads((split_root / "split/receipt.json").read_text())
    reference_manifest = json.loads((reference_root / "ARTIFACTS.json").read_text())
    if (
        split_receipt.get("schema_version")
        != "masld-bench-gse281364-outcome-blind-splits-v1"
        or split_receipt.get("status") != "pass_outcome_blind_split_contract"
        or split_receipt.get("outcomes_read") is not False
        or split_receipt.get("reporter_counts_read") is not False
        or reference_manifest.get("metadata", {}).get("build") != "GRCh38.p14"
        or reference_manifest.get("metadata", {}).get("indexed") is not True
    ):
        raise MPRALegNetFixtureError("split or reference receipt differs")
    rows = read_rows(split_root / "split/elements.tsv")
    reference = IndexedFasta(
        reference_root / "GRCh38.p14.sequence_model.fa",
        reference_root / "GRCh38.p14.sequence_model.fa.fai",
    )
    output.mkdir(parents=True)
    manifest_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []
    for row in rows:
        ref_sequence, alt_sequence, start, end = allele_window(reference, row)
        ref_rc = reverse_complement(ref_sequence)
        alt_rc = reverse_complement(alt_sequence)
        fixture_id = "mpralegnet_" + sha256(row["element_id"].encode()).hexdigest()[:24]
        fasta_records.extend(
            (
                (f"{fixture_id}|REF", ref_sequence),
                (f"{fixture_id}|ALT", alt_sequence),
                (f"{fixture_id}|REF_RC", ref_rc),
                (f"{fixture_id}|ALT_RC", alt_rc),
            )
        )
        manifest_rows.append(
            {
                "fixture_id": fixture_id,
                "element_id": row["element_id"],
                "outer_locus_sequence_group_id": row["outer_locus_sequence_group_id"],
                "outer_fold": row["outer_fold"],
                "contig": row["contig"],
                "variant_pos0": row["variant_pos0"],
                "variant_pos1": row["variant_pos1"],
                "input_start0": start,
                "input_end0": end,
                "forward_variant_index0": FORWARD_VARIANT_INDEX0,
                "reverse_complement_variant_index0": REVERSE_VARIANT_INDEX0,
                "ref": row["genomic_ref"],
                "alt": row["genomic_alt"],
                "cell_context": "HepG2",
                "input_lane": "GRCh38p14_230bp_genomic_window_transfer",
                "source_assay_insert_length_bp": 107,
                "source_oligo_difference_index0": row["oligo_difference_index0"],
                "reference_sequence_sha256": sha256(ref_sequence.encode()).hexdigest(),
                "alternative_sequence_sha256": sha256(alt_sequence.encode()).hexdigest(),
                "reference_reverse_complement_sha256": sha256(ref_rc.encode()).hexdigest(),
                "alternative_reverse_complement_sha256": sha256(alt_rc.encode()).hexdigest(),
            }
        )
    write_fasta(output / "mpralegnet.alleles.fa.gz", fasta_records)
    with (output / "manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    fold_counts = {
        str(fold): sum(int(row["outer_fold"]) == fold for row in manifest_rows)
        for fold in range(5)
    }
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-mpralegnet-fixture-v1",
        "status": "pass_outcome_blind_mpralegnet_fixture",
        "dataset_id": "gse281364",
        "model_id": "mpralegnet_hepg2_test1_val2",
        "elements": len(manifest_rows),
        "outer_locus_sequence_groups": len(
            {row["outer_locus_sequence_group_id"] for row in manifest_rows}
        ),
        "sequences": len(fasta_records),
        "input_length_bp": WINDOW,
        "forward_variant_index0": FORWARD_VARIANT_INDEX0,
        "reverse_complement_variant_index0": REVERSE_VARIANT_INDEX0,
        "fold_elements": fold_counts,
        "cell_context": "HepG2",
        "native_output": "uncalibrated_lentiMPRA_reporter_expression_score",
        "input_lane": "GRCh38p14_230bp_genomic_window_transfer",
        "assay_insert_input_status": "unsupported_107bp_not_padded_or_truncated",
        "strand_policy": "score_forward_and_reverse_complement_then_average",
        "allele_effect_sign": "ALT_minus_REF",
        "split_artifacts_sha256": split_sha256,
        "reference_artifacts_sha256": reference_sha256,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "standalone_champion_eligible": False,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--split-sha256", required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
