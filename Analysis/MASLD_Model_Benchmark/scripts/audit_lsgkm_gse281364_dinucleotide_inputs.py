#!/usr/bin/env python3
"""Independently audit emitted LS-GKM exact-dinucleotide-null FASTAs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

from scripts.freeze_lsgkm_gse281364_dinucleotide_null_readiness import (
    SEEDS,
    canonical_sequence,
    dinucleotide_counts,
    file_sha256,
    monomer_counts,
)


class DinucleotideInputAuditError(RuntimeError):
    """Raised when emitted materialization output files fail an invariant."""


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise DinucleotideInputAuditError(f"JSON object differs: {path}")
    return value


def read_tsv_gzip(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if not rows:
        raise DinucleotideInputAuditError(f"TSV is empty: {path}")
    return rows


def read_fasta_gzip(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence_parts: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence_parts)))
                identifier = line[1:]
                sequence_parts = []
            elif identifier is None or not line:
                raise DinucleotideInputAuditError(f"FASTA structure differs: {path}")
            else:
                sequence_parts.append(line)
    if identifier is not None:
        records.append((identifier, "".join(sequence_parts)))
    if not records or len({identifier for identifier, _sequence in records}) != len(records):
        raise DinucleotideInputAuditError(f"FASTA roster differs: {path}")
    return records


def audit_pair(positive: str, negative: str) -> None:
    if (
        len(positive) != 300
        or len(negative) != 300
        or positive[0] != negative[0]
        or positive[-1] != negative[-1]
        or monomer_counts(positive) != monomer_counts(negative)
        or dinucleotide_counts(positive) != dinucleotide_counts(negative)
        or canonical_sequence(positive) == canonical_sequence(negative)
    ):
        raise DinucleotideInputAuditError("emitted positive/negative invariant differs")


def build_audit(materialization: Path) -> dict[str, Any]:
    contract_path = materialization / "contract/materialization_contract.json"
    contract = load_json(contract_path)
    if (
        contract.get("status")
        != "inputs_materialized_pending_independent_emitted_fasta_audit"
        or contract.get("design_id")
        != "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1"
        or contract.get("split_count") != 5
        or contract.get("shared_fit_count") != 25
        or contract.get("outcomes_read") is not False
        or contract.get("prediction_values_read") is not False
        or contract.get("production_fits_executed") != 0
    ):
        raise DinucleotideInputAuditError("materialization contract differs")

    scoring_manifest_path = materialization / "contract/scoring_elements.tsv.gz"
    scoring_fasta_path = materialization / "scoring/alleles.fa.gz"
    scoring_contract = contract["scoring_universe"]
    if (
        file_sha256(scoring_manifest_path)
        != scoring_contract["scoring_element_manifest_sha256"]
        or file_sha256(scoring_fasta_path)
        != scoring_contract["scoring_allele_fasta_sha256"]
    ):
        raise DinucleotideInputAuditError("scoring materialization hash differs")
    scoring_rows = read_tsv_gzip(scoring_manifest_path)
    scoring_records = read_fasta_gzip(scoring_fasta_path)
    if (
        len(scoring_rows) != 1033
        or len(scoring_records) != 2066
        or len({row["long_range_block_id"] for row in scoring_rows}) != 239
    ):
        raise DinucleotideInputAuditError("scoring materialization denominator differs")
    scoring_map = dict(scoring_records)
    scoring_canonical: set[str] = set()
    scoring_split_by_element: dict[str, str] = {}
    for row in scoring_rows:
        element = row["element_id"]
        reference = scoring_map.get(f"{element}|REF")
        alternative = scoring_map.get(f"{element}|ALT")
        if reference is None or alternative is None:
            raise DinucleotideInputAuditError("scoring FASTA REF/ALT pair is missing")
        index = int(row["variant_index0"])
        if (
            len(reference) != 300
            or len(alternative) != 300
            or index != 150
            or reference[index] != row["ref"]
            or alternative[index] != row["alt"]
            or reference[:index] != alternative[:index]
            or reference[index + 1 :] != alternative[index + 1 :]
            or hashlib.sha256(reference.encode("ascii")).hexdigest()
            != row["reference_sequence_sha256"]
            or hashlib.sha256(alternative.encode("ascii")).hexdigest()
            != row["alternative_sequence_sha256"]
            or row["allele_effect_sign"] != "ALT_minus_REF"
        ):
            raise DinucleotideInputAuditError("emitted scoring REF/ALT invariant differs")
        scoring_canonical.update(
            (canonical_sequence(reference), canonical_sequence(alternative))
        )
        scoring_split_by_element[element] = row["scoring_split_id"]

    summaries = contract["split_summaries"]
    if len(summaries) != 5 or {summary["split_id"] for summary in summaries} != {
        f"donor{fold}_genomic{fold}" for fold in range(5)
    }:
        raise DinucleotideInputAuditError("split summary roster differs")
    held_test_by_split = {
        summary["split_id"]: set(summary["held_test_contigs"]) for summary in summaries
    }
    scoring_contig_by_element = {row["element_id"]: row["contig"] for row in scoring_rows}
    if any(
        scoring_contig_by_element[element]
        not in held_test_by_split[scoring_split_by_element[element]]
        for element in scoring_split_by_element
    ):
        raise DinucleotideInputAuditError("scoring element is not held by its fit split")

    total_pairs = 0
    for summary in summaries:
        split_id = summary["split_id"]
        positive_path = materialization / "inputs" / f"{split_id}.positive.fa.gz"
        if file_sha256(positive_path) != summary["positive_fasta_sha256"]:
            raise DinucleotideInputAuditError("positive FASTA hash differs")
        positives = read_fasta_gzip(positive_path)
        if len(positives) != 10000 or summary["selected_positive_windows"] != 10000:
            raise DinucleotideInputAuditError("positive FASTA count differs")
        positive_ids = [identifier for identifier, _sequence in positives]
        positive_sequences = [sequence for _identifier, sequence in positives]
        positive_canonical = {canonical_sequence(sequence) for sequence in positive_sequences}
        if len(positive_canonical) != 10000 or positive_canonical & scoring_canonical:
            raise DinucleotideInputAuditError("positive canonical collision differs")

        negative_by_seed: dict[int, list[str]] = {}
        for seed in SEEDS:
            seed_summary = summary["seeds"][str(seed)]
            negative_path = (
                materialization / "inputs" / f"{split_id}.seed{seed}.negative.fa.gz"
            )
            pair_path = (
                materialization / "inputs" / f"{split_id}.seed{seed}.pairs.tsv.gz"
            )
            if (
                file_sha256(negative_path) != seed_summary["negative_fasta_sha256"]
                or file_sha256(pair_path) != seed_summary["pair_manifest_sha256"]
            ):
                raise DinucleotideInputAuditError("seed input hash differs")
            negatives = read_fasta_gzip(negative_path)
            pairs = read_tsv_gzip(pair_path)
            if len(negatives) != 10000 or len(pairs) != 10000:
                raise DinucleotideInputAuditError("seed pair count differs")
            negative_sequences = [sequence for _identifier, sequence in negatives]
            negative_canonical: set[str] = set()
            for index, ((negative_id, negative), pair) in enumerate(zip(negatives, pairs)):
                positive_id, positive = positives[index]
                audit_pair(positive, negative)
                canonical = canonical_sequence(negative)
                if (
                    pair["split_id"] != split_id
                    or int(pair["seed"]) != seed
                    or int(pair["pair_index"]) != index
                    or pair["positive_id"] != positive_id
                    or pair["negative_id"] != negative_id
                    or pair["contig"] in held_test_by_split[split_id]
                    or hashlib.sha256(positive.encode("ascii")).hexdigest()
                    != pair["positive_sequence_sha256"]
                    or hashlib.sha256(negative.encode("ascii")).hexdigest()
                    != pair["negative_sequence_sha256"]
                    or canonical in negative_canonical
                    or canonical in positive_canonical
                    or canonical in scoring_canonical
                ):
                    raise DinucleotideInputAuditError("emitted pair manifest differs")
                negative_canonical.add(canonical)
            negative_by_seed[seed] = negative_sequences
            total_pairs += len(pairs)
        for index in range(10000):
            if len(
                {
                    canonical_sequence(negative_by_seed[seed][index])
                    for seed in SEEDS
                }
            ) != 5:
                raise DinucleotideInputAuditError("cross-seed shuffle diversity differs")

    if total_pairs != 250000:
        raise DinucleotideInputAuditError("full five-seed rectangle differs")
    return {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-input-audit-v1",
        "status": "pass_independent_emitted_fasta_invariant_audit",
        "design_id": contract["design_id"],
        "materialization_contract_sha256": file_sha256(contract_path),
        "split_count": 5,
        "fit_count": 25,
        "positive_FASTA_count": 5,
        "negative_FASTA_count": 25,
        "pair_manifest_count": 25,
        "positive_windows_per_split": 10000,
        "negative_windows_per_split_seed": 10000,
        "audited_pair_rows": total_pairs,
        "scoring_elements": 1033,
        "scoring_long_range_blocks": 239,
        "scoring_REF_ALT_pairs": 1033,
        "sequence_length_bp": 300,
        "monomer_counts_rechecked_from_emitted_FASTA": True,
        "dinucleotide_counts_rechecked_from_emitted_FASTA": True,
        "first_last_base_rechecked_from_emitted_FASTA": True,
        "canonical_collisions_rechecked_from_emitted_FASTA": True,
        "five_seed_diversity_rechecked_from_emitted_FASTA": True,
        "REF_ALT_pairing_rechecked_from_emitted_FASTA": True,
        "held_test_contig_exclusion_rechecked_from_emitted_manifests": True,
        "production_training_authorized": False,
        "production_fits_executed": 0,
        "production_predictions_generated": 0,
        "benchmark_metrics_calculated": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "prediction_values_read": False,
        "sealed_assets_read": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--materialization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise DinucleotideInputAuditError("refusing to overwrite audit receipt")
    audit = build_audit(arguments.materialization.resolve(strict=True))
    arguments.output.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
