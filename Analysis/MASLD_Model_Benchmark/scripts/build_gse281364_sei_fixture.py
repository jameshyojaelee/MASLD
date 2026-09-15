#!/usr/bin/env python3
"""Build outcome-blind 4,096-bp REF/ALT/RC inputs for the Sei MPRA screen."""

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


WINDOW = 4_096
CENTER = WINDOW // 2
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
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reference_reverse_complement_sha256",
    "alternative_reverse_complement_sha256",
)


class SeiFixtureError(ValueError):
    """Raised when the Sei outcome-blind fixture requirement differs."""


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
    start, end = position - CENTER, position - CENTER + WINDOW
    sequence = reference.fetch(row["contig"], start, end)
    if (
        len(sequence) != WINDOW
        or set(sequence) - set("ACGT")
        or sequence[CENTER] != row["genomic_ref"]
        or row["genomic_ref"] == row["genomic_alt"]
        or row["genomic_alt"] not in "ACGT"
    ):
        raise SeiFixtureError("reference window or allele differs")
    alternative = sequence[:CENTER] + row["genomic_alt"] + sequence[CENTER + 1 :]
    if sum(left != right for left, right in zip(sequence, alternative)) != 1:
        raise SeiFixtureError("alternative window differs at more than one base")
    return sequence, alternative, start, end


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as handle:
                for name, sequence in records:
                    handle.write(f">{name}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if len(rows) != 4_359:
        raise SeiFixtureError("split element count differs")
    return rows


def select_one_per_group(rows: list[dict[str, str]], seed: str) -> list[dict[str, str]]:
    selected: dict[str, tuple[str, dict[str, str]]] = {}
    for row in rows:
        group = row["outer_locus_sequence_group_id"]
        key = sha256(f"{seed}\0{row['element_id']}".encode()).hexdigest()
        if group not in selected or key < selected[group][0]:
            selected[group] = (key, row)
    if len(selected) != 1_033:
        raise SeiFixtureError("smoke group count differs")
    return [selected[group][1] for group in sorted(selected)]


def build(
    *,
    split_root: Path,
    split_sha256: str,
    reference_root: Path,
    reference_sha256: str,
    output: Path,
    seed: str,
) -> dict[str, object]:
    if output.exists():
        raise SeiFixtureError("output exists")
    if digest(split_root / "ARTIFACTS.json") != split_sha256:
        raise SeiFixtureError("split artifact changed")
    if digest(reference_root / "ARTIFACTS.json") != reference_sha256:
        raise SeiFixtureError("reference artifact changed")
    split_receipt = json.loads((split_root / "split/receipt.json").read_text())
    reference_manifest = json.loads((reference_root / "ARTIFACTS.json").read_text())
    if (
        split_receipt.get("status") != "pass_outcome_blind_split_contract"
        or split_receipt.get("outcomes_read") is not False
        or reference_manifest.get("metadata", {}).get("build") != "GRCh38.p14"
        or reference_manifest.get("metadata", {}).get("indexed") is not True
    ):
        raise SeiFixtureError("split or reference receipt differs")
    rows = select_one_per_group(read_rows(split_root / "split/elements.tsv"), seed)
    reference = IndexedFasta(
        reference_root / "GRCh38.p14.sequence_model.fa",
        reference_root / "GRCh38.p14.sequence_model.fa.fai",
    )
    output.mkdir(parents=True)
    manifest_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []
    for row in rows:
        ref_sequence, alt_sequence, start, end = allele_window(reference, row)
        fixture_id = "sei_" + sha256(row["element_id"].encode()).hexdigest()[:24]
        ref_rc, alt_rc = reverse_complement(ref_sequence), reverse_complement(alt_sequence)
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
                "variant_index0": CENTER,
                "ref": row["genomic_ref"],
                "alt": row["genomic_alt"],
                "reference_sequence_sha256": sha256(ref_sequence.encode()).hexdigest(),
                "alternative_sequence_sha256": sha256(alt_sequence.encode()).hexdigest(),
                "reference_reverse_complement_sha256": sha256(ref_rc.encode()).hexdigest(),
                "alternative_reverse_complement_sha256": sha256(alt_rc.encode()).hexdigest(),
            }
        )
    write_fasta(output / "sei.alleles.fa.gz", fasta_records)
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
        "schema_version": "masld-bench-gse281364-sei-fixture-v1",
        "status": "pass_outcome_blind_sei_fixture",
        "dataset_id": "gse281364",
        "model_id": "sei",
        "elements": len(manifest_rows),
        "outer_locus_sequence_groups": len(manifest_rows),
        "sequences": len(fasta_records),
        "input_length_bp": WINDOW,
        "variant_index0": CENTER,
        "fold_elements": fold_counts,
        "selection": "one_deterministic_element_per_outer_locus_sequence_group",
        "selection_seed": seed,
        "split_artifacts_sha256": split_sha256,
        "reference_artifacts_sha256": reference_sha256,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "champion_eligible": False,
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
    parser.add_argument("--seed", default="gse281364-sei-smoke-v1")
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
