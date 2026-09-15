#!/usr/bin/env python3
"""Independently audit frozen-shape LS-GKM input activation outputs."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
STEM = re.compile(r"^(donor([0-4])_genomic\2)\.seed(1103|2909|4721|6673|8111)$")
MATCH_FIELDS = (
    "seed",
    "pair_id",
    "split_id",
    "genomic_test_fold",
    "positive_id",
    "negative_id",
    "positive_contig",
    "negative_contig",
    "positive_start0",
    "positive_end0",
    "negative_start0",
    "negative_end0",
    "positive_atac_decile",
    "negative_atac_decile",
    "positive_tss_distance_bin",
    "negative_tss_distance_bin",
    "positive_gc_fraction",
    "negative_gc_fraction",
    "positive_mappability",
    "negative_mappability",
    "positive_log1p_atac",
    "negative_log1p_atac",
    "positive_tss_distance",
    "negative_tss_distance",
    "positive_repeat_fraction",
    "negative_repeat_fraction",
    "standardized_covariate_distance",
    "positive_sequence_sha256",
    "negative_sequence_sha256",
)


class ActivationAuditError(RuntimeError):
    """Raised when a materialized input differs from the activation requirement."""


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ActivationAuditError(f"JSON object differs: {path}")
    return value


def hash_id(*parts: object) -> str:
    return hashlib.sha256("\0".join(map(str, parts)).encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical(sequence: str) -> str:
    return min(sequence, reverse_complement(sequence))


def read_fasta(path: Path, suffix: str) -> dict[str, str]:
    records: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            if not header.startswith(">") or not sequence or handle.closed:
                raise ActivationAuditError(f"FASTA structure differs: {path}")
            identifier = header[1:].rstrip("\n")
            expected_suffix = f"|{suffix}"
            if not identifier.endswith(expected_suffix):
                raise ActivationAuditError(f"FASTA role differs: {path}")
            pair_id = identifier[: -len(expected_suffix)]
            if pair_id in records or len(sequence) != 300 or set(sequence) - set("ACGT"):
                raise ActivationAuditError(f"FASTA sequence differs: {path}")
            records[pair_id] = sequence
    if not records:
        raise ActivationAuditError(f"FASTA is empty: {path}")
    return records


def audit_scoring_map(path: Path) -> dict[str, Any]:
    counts: Counter[tuple[str, str]] = Counter()
    elements: set[str] = set()
    rows = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seed = int(row["seed"])
            context = row["assay_context_id"]
            split_match = re.fullmatch(r"donor([0-4])_genomic([0-4])", row["scoring_split_id"])
            if (
                seed not in SEEDS
                or context not in CONTEXTS
                or row["study_id"] != "gse281364"
                or row["model_state_id"] != "hepatocyte"
                or row["context_specific_prediction"] != "false"
                or split_match is None
                or split_match.group(1) != split_match.group(2)
                or int(split_match.group(2)) != int(row["genomic_test_fold"])
            ):
                raise ActivationAuditError("scoring-row topology differs")
            counts[(row["element_id"], context)] += 1
            elements.add(row["element_id"])
            rows += 1
    if rows != 10_330 or len(elements) != 1_033 or set(counts.values()) != {5}:
        raise ActivationAuditError("scoring-row denominator differs")
    return {"rows": rows, "elements": len(elements), "contexts": list(CONTEXTS)}


def audit_fit(stem: Path, contract_summary: dict[str, Any]) -> dict[str, Any]:
    match = STEM.fullmatch(stem.name)
    if match is None:
        raise ActivationAuditError(f"fit stem differs: {stem.name}")
    split_id, fold_text, seed_text = match.group(1), match.group(2), match.group(3)
    fold, seed = int(fold_text), int(seed_text)
    match_path = Path(f"{stem}.matches.tsv.gz")
    positive_path = Path(f"{stem}.positive.fa.gz")
    negative_path = Path(f"{stem}.negative.fa.gz")
    positive = read_fasta(positive_path, "positive")
    negative = read_fasta(negative_path, "negative")
    rows = 0
    negative_ids: set[str] = set()
    positive_sequences: set[str] = set()
    negative_sequences: set[str] = set()
    with gzip.open(match_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MATCH_FIELDS:
            raise ActivationAuditError(f"match fields differ: {stem.name}")
        for row in reader:
            pair_id = row["pair_id"]
            positive_sequence = positive.get(pair_id)
            negative_sequence = negative.get(pair_id)
            if positive_sequence is None or negative_sequence is None:
                raise ActivationAuditError(f"FASTA pair join differs: {stem.name}")
            if (
                int(row["seed"]) != seed
                or row["split_id"] != split_id
                or int(row["genomic_test_fold"]) != fold
                or row["positive_contig"] != row["negative_contig"]
                or row["positive_atac_decile"] != row["negative_atac_decile"]
                or row["positive_tss_distance_bin"] != row["negative_tss_distance_bin"]
                or int(row["positive_end0"]) - int(row["positive_start0"]) != 300
                or int(row["negative_end0"]) - int(row["negative_start0"]) != 300
                or pair_id
                != hash_id(
                    "lsgkm-pair-v1",
                    seed,
                    split_id,
                    row["positive_id"],
                    row["negative_id"],
                )
                or hashlib.sha256(positive_sequence.encode("ascii")).hexdigest()
                != row["positive_sequence_sha256"]
                or hashlib.sha256(negative_sequence.encode("ascii")).hexdigest()
                != row["negative_sequence_sha256"]
            ):
                raise ActivationAuditError(f"match identity differs: {stem.name}")
            gc_delta = abs(float(row["positive_gc_fraction"]) - float(row["negative_gc_fraction"]))
            map_delta = abs(float(row["positive_mappability"]) - float(row["negative_mappability"]))
            repeat_delta = abs(float(row["positive_repeat_fraction"]) - float(row["negative_repeat_fraction"]))
            distance = gc_delta / 0.02 + map_delta / 0.05 + repeat_delta / 0.05
            if (
                not 0 <= int(row["positive_atac_decile"]) <= 9
                or not 0 <= int(row["positive_tss_distance_bin"]) <= 4
                or float(row["positive_mappability"]) < 0.80
                or float(row["negative_mappability"]) < 0.80
                or gc_delta > 0.02 + 1e-9
                or map_delta > 0.05 + 1e-9
                or repeat_delta > 0.05 + 1e-9
                or not math.isclose(
                    distance,
                    float(row["standardized_covariate_distance"]),
                    rel_tol=0.0,
                    abs_tol=2e-7,
                )
                or row["negative_id"] in negative_ids
            ):
                raise ActivationAuditError(f"match tolerance or replacement differs: {stem.name}")
            negative_ids.add(row["negative_id"])
            positive_sequences.add(canonical(positive_sequence))
            negative_sequences.add(canonical(negative_sequence))
            rows += 1
    seed_summary = contract_summary["seeds"][str(seed)]
    expected = int(seed_summary["matched_pairs"])
    if (
        rows != expected
        or file_sha256(match_path) != seed_summary["match_manifest_sha256"]
        or file_sha256(positive_path) != seed_summary["positive_fasta_sha256"]
        or file_sha256(negative_path) != seed_summary["negative_fasta_sha256"]
        or set(positive) != set(negative)
        or len(positive) != rows
        or positive_sequences.intersection(negative_sequences)
        or len(positive_sequences) != rows
        or len(negative_sequences) != rows
    ):
        raise ActivationAuditError(f"fit denominator or sequence leakage differs: {stem.name}")
    return {"split_id": split_id, "seed": seed, "matched_pairs": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    root = arguments.activation.resolve(strict=True)
    if arguments.output.exists():
        raise ActivationAuditError("audit output exists")
    contract = load_json(root / "contract/production_activation_contract.json")
    if (
        contract.get("status") != "inputs_frozen_production_fit_not_executed"
        or contract.get("split_design", {}).get("shared_fits") != 25
        or contract.get("active_state_roster") != ["hepatocyte"]
        or contract.get("context_specific_prediction") is not False
        or contract.get("outcomes_read") is not False
        or contract.get("reporter_counts_read") is not False
        or contract.get("sealed_assets_read") is not False
        or contract.get("production_training_executed") is not False
        or contract.get("production_predictions_generated") is not False
    ):
        raise ActivationAuditError("activation contract firewall differs")
    scoring = audit_scoring_map(root / "contract/scoring_row_map.tsv.gz")
    summaries = {summary["split_id"]: summary for summary in contract["input_summary"]}
    expected_stems = {
        f"donor{fold}_genomic{fold}.seed{seed}"
        for fold in range(5)
        for seed in SEEDS
    }
    observed_stems = {
        path.name.removesuffix(".matches.tsv.gz")
        for path in (root / "inputs").glob("*.matches.tsv.gz")
    }
    if observed_stems != expected_stems or set(summaries) != {
        f"donor{fold}_genomic{fold}" for fold in range(5)
    }:
        raise ActivationAuditError("fit roster differs")
    fits = [
        audit_fit(root / "inputs" / stem, summaries[stem.split(".seed", 1)[0]])
        for stem in sorted(observed_stems)
    ]
    receipt = {
        "schema_version": "masld-bench-lsgkm-gse281364-input-activation-audit-v1",
        "status": "pass_outcome_blind_input_activation",
        "scoring": scoring,
        "fits": fits,
        "fit_count": len(fits),
        "matched_pairs_minimum": min(fit["matched_pairs"] for fit in fits),
        "matched_pairs_maximum": max(fit["matched_pairs"] for fit in fits),
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_assets_read": False,
        "production_training_executed": False,
        "production_predictions_generated": False,
    }
    arguments.output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
