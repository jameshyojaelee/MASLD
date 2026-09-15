#!/usr/bin/env python3
"""Build outcome-blind 131,072-bp Caduceus windows for fixed MPRA loci."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from scripts.alphagenome_sei_build_fixture import IndexedFasta


WINDOW = 131_072
CENTER = WINDOW // 2
POOL = 1_536
POOL_START = CENTER - POOL // 2
POOL_END = POOL_START + POOL
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
FIELDS = (
    "fixture_id",
    "source_fixture_id",
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "variant_pos1",
    "input_start0",
    "input_end0",
    "input_length_bp",
    "forward_variant_index0",
    "reverse_complement_variant_index0",
    "pool_start0",
    "pool_end0",
    "pool_length_bp",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reverse_complement_reference_sha256",
    "reverse_complement_alternative_sha256",
)


class CaduceusFixtureError(ValueError):
    """Raised when a long-context fixture invariant differs."""


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
    ref = row.get("ref", row.get("genomic_ref", ""))
    alt = row.get("alt", row.get("genomic_alt", ""))
    if (
        len(sequence) != WINDOW
        or set(sequence) - set("ACGTN")
        or sequence[CENTER] != ref
        or ref == alt
        or alt not in "ACGT"
    ):
        raise CaduceusFixtureError("reference window or allele differs")
    alternative = sequence[:CENTER] + alt + sequence[CENTER + 1 :]
    if sum(left != right for left, right in zip(sequence, alternative)) != 1:
        raise CaduceusFixtureError("alternative window differs at more than one base")
    return sequence, alternative, start, end


def build(
    *,
    enformer_fixture_root: Path,
    enformer_fixture_sha256: str,
    fasta_root: Path,
    fasta_sha256: str,
    output: Path,
) -> dict[str, object]:
    if output.exists() or any(
        path.is_symlink() for path in (enformer_fixture_root, fasta_root)
    ):
        raise CaduceusFixtureError("fixture request differs")
    if (
        digest(enformer_fixture_root / "ARTIFACTS.json") != enformer_fixture_sha256
        or digest(fasta_root / "ARTIFACTS.json") != fasta_sha256
    ):
        raise CaduceusFixtureError("upstream artifact changed")
    enformer_receipt = json.loads(
        (enformer_fixture_root / "fixture/receipt.json").read_text(encoding="utf-8")
    )
    fasta_artifacts = json.loads(
        (fasta_root / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    if (
        enformer_receipt.get("status") != "pass_outcome_blind_enformer_fixture"
        or enformer_receipt.get("elements") != 1_033
        or enformer_receipt.get("outer_locus_sequence_groups") != 1_033
        or any(
            enformer_receipt.get(field) is not False
            for field in ("outcomes_read", "reporter_counts_read", "sealed_outcomes_read")
        )
        or fasta_artifacts.get("metadata", {}).get("build") != "GRCh38.p14"
        or fasta_artifacts.get("metadata", {}).get("indexed") is not True
    ):
        raise CaduceusFixtureError("upstream outcome firewall or reference differs")
    source_rows = read_tsv(enformer_fixture_root / "fixture/manifest.tsv")
    if (
        len(source_rows) != 1_033
        or len({row["element_id"] for row in source_rows}) != 1_033
        or len({row["outer_locus_sequence_group_id"] for row in source_rows}) != 1_033
        or {int(row["outer_fold"]) for row in source_rows} != set(range(5))
    ):
        raise CaduceusFixtureError("source locus census differs")

    reference = IndexedFasta(
        fasta_root / "GRCh38.p14.sequence_model.fa",
        fasta_root / "GRCh38.p14.sequence_model.fa.fai",
    )
    manifest_rows: list[dict[str, object]] = []
    for source in source_rows:
        ref_sequence, alt_sequence, start, end = allele_window(reference, source)
        ref_rc = reverse_complement(ref_sequence)
        alt_rc = reverse_complement(alt_sequence)
        if (
            ref_rc[CENTER - 1] != reverse_complement(source["ref"])
            or alt_rc[CENTER - 1] != reverse_complement(source["alt"])
        ):
            raise CaduceusFixtureError("reverse-complement variant geometry differs")
        manifest_rows.append(
            {
                "fixture_id": "caduceus_" + sha256(source["element_id"].encode()).hexdigest()[:24],
                "source_fixture_id": source["fixture_id"],
                "element_id": source["element_id"],
                "outer_locus_sequence_group_id": source["outer_locus_sequence_group_id"],
                "outer_fold": source["outer_fold"],
                "contig": source["contig"],
                "variant_pos0": source["variant_pos0"],
                "variant_pos1": source["variant_pos1"],
                "input_start0": start,
                "input_end0": end,
                "input_length_bp": WINDOW,
                "forward_variant_index0": CENTER,
                "reverse_complement_variant_index0": CENTER - 1,
                "pool_start0": POOL_START,
                "pool_end0": POOL_END,
                "pool_length_bp": POOL,
                "ref": source["ref"],
                "alt": source["alt"],
                "reference_sequence_sha256": sha256(ref_sequence.encode()).hexdigest(),
                "alternative_sequence_sha256": sha256(alt_sequence.encode()).hexdigest(),
                "reverse_complement_reference_sha256": sha256(ref_rc.encode()).hexdigest(),
                "reverse_complement_alternative_sha256": sha256(alt_rc.encode()).hexdigest(),
            }
        )

    output.mkdir(parents=True)
    with (output / "sequence_manifest.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-caduceus-131k-fixture-v1",
        "status": "pass_outcome_blind_caduceus_131k_fixture",
        "dataset_id": "gse281364",
        "model_id": "caduceus",
        "elements": len(manifest_rows),
        "outer_locus_sequence_groups": len(manifest_rows),
        "outer_folds": 5,
        "input_length_bp": WINDOW,
        "forward_variant_index0": CENTER,
        "reverse_complement_variant_index0": CENTER - 1,
        "pool_start0": POOL_START,
        "pool_end0": POOL_END,
        "pool_length_bp": POOL,
        "selection": "exact_same_elements_and_locus_groups_as_frozen_enformer_and_common_DNA_LM_fixtures",
        "reference_build": "GRCh38.p14",
        "enformer_fixture_artifacts_sha256": enformer_fixture_sha256,
        "fasta_artifacts_sha256": fasta_sha256,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "checkpoint_bytes_loaded": False,
        "model_forward_executed": False,
        "downstream_head_fit": False,
        "champion_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enformer-fixture-root", type=Path, required=True)
    parser.add_argument("--enformer-fixture-sha256", required=True)
    parser.add_argument("--fasta-root", type=Path, required=True)
    parser.add_argument("--fasta-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args())), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
