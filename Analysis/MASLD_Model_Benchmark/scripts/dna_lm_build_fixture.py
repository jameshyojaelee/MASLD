#!/usr/bin/env python3
"""Derive identical outcome-free DNA-LM REF/ALT/RC fixtures."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Iterable, Mapping


COMPLEMENT = str.maketrans("ACGT", "TGCA")
MANIFEST_FIELDS = (
    "fixture_id",
    "model_id",
    "context_id",
    "genomic_fold",
    "contig",
    "anchor0",
    "variant_pos_1based",
    "ref",
    "alt",
    "input_length_bp",
    "forward_variant_index0",
    "reverse_complement_variant_index0",
    "pool_start0",
    "pool_end0",
    "pool_length_bp",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reverse_complement_reference_sha256",
    "reverse_complement_alternative_sha256",
    "source_reference_sequence_sha256",
    "source_alternative_sequence_sha256",
    "allele_effect_sign",
    "native_task_contract",
)


class DnaLmFixtureError(ValueError):
    """Raised when a DNA-LM fixture violates its frozen contract."""


def _digest_bytes(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def _digest_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def _read_fasta_gz(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            value = line.rstrip("\n")
            if value.startswith(">"):
                name = value[1:]
                if not name or name in records:
                    raise DnaLmFixtureError("upstream FASTA header differs")
                records[name] = []
            elif name is None:
                raise DnaLmFixtureError("upstream FASTA sequence precedes header")
            else:
                records[name].append(value)
    return {name: "".join(parts) for name, parts in records.items()}


def _write_fasta_gz(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as handle:
                for name, sequence in records:
                    handle.write(f">{name}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    required = {
        "fixture_id",
        "model_id",
        "genomic_fold",
        "contig",
        "anchor0",
        "variant_pos_1based",
        "ref",
        "alt",
        "input_length",
        "reference_sequence_sha256",
        "alternative_sequence_sha256",
    }
    if not rows or not required.issubset(rows[0]):
        raise DnaLmFixtureError("upstream sequence manifest differs")
    selected = [row for row in rows if row["model_id"] == "alphagenome_all_folds"]
    if len(selected) != 3 or {int(row["genomic_fold"]) for row in selected} != {2, 3, 4}:
        raise DnaLmFixtureError("upstream anchor census differs")
    return sorted(selected, key=lambda row: int(row["genomic_fold"]))


def _slice(sequence: str, length: int) -> str:
    if len(sequence) < length or (len(sequence) - length) % 2:
        raise DnaLmFixtureError("context cannot be centered exactly")
    start = (len(sequence) - length) // 2
    return sequence[start : start + length]


def _validate_alleles(ref: str, alt: str, rc_ref: str, rc_alt: str) -> int:
    if (
        len({len(ref), len(alt), len(rc_ref), len(rc_alt)}) != 1
        or set(ref + alt + rc_ref + rc_alt) - set("ACGT")
        or rc_ref != _reverse_complement(ref)
        or rc_alt != _reverse_complement(alt)
    ):
        raise DnaLmFixtureError("allele or reverse-complement sequence differs")
    differences = [index for index, pair in enumerate(zip(ref, alt)) if pair[0] != pair[1]]
    if differences != [len(ref) // 2]:
        raise DnaLmFixtureError("forward variant is not the single right-center change")
    rc_differences = [
        index for index, pair in enumerate(zip(rc_ref, rc_alt)) if pair[0] != pair[1]
    ]
    if rc_differences != [len(ref) // 2 - 1]:
        raise DnaLmFixtureError("reverse-complement variant is not realigned")
    return differences[0]


def build_fixture(config_path: Path, upstream: Path, output: Path) -> dict[str, object]:
    if output.exists() or config_path.is_symlink() or upstream.is_symlink():
        raise DnaLmFixtureError("fixture request differs")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "masld-bench-dna-lm-fixture-v1":
        raise DnaLmFixtureError("fixture config schema differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or any(firewall.values()):
        raise DnaLmFixtureError("fixture firewall differs")
    fixture_root = upstream / "fixture"
    manifest_path = fixture_root / "sequence_manifest.tsv"
    fasta_path = fixture_root / "alphagenome_all_folds.alleles.fa.gz"
    for path in (manifest_path, fasta_path):
        if path.is_symlink() or not path.is_file():
            raise DnaLmFixtureError("upstream fixture input differs")
    source_rows = _read_manifest(manifest_path)
    source_sequences = _read_fasta_gz(fasta_path)
    output.mkdir(mode=0o750)

    contexts: dict[str, dict[str, object]] = {
        config["common_lane"]["context_id"]: {
            "input_length_bp": config["common_lane"]["input_length_bp"],
            "models": config["common_lane"]["models"],
        }
    }
    for context_id, record in config["native_contexts"].items():
        contexts[context_id] = {
            "input_length_bp": record["input_length_bp"],
            "models": [record["model_id"]],
        }

    fasta_records: dict[str, list[tuple[str, str]]] = {
        context_id: [] for context_id in contexts
    }
    manifest_rows: list[dict[str, object]] = []
    for source in source_rows:
        fixture_id = source["fixture_id"]
        prefix = f"{fixture_id}|alphagenome_all_folds|"
        sequences = {
            allele: source_sequences[prefix + allele]
            for allele in ("REF", "ALT", "REF_RC", "ALT_RC")
        }
        if (
            len(sequences["REF"]) != config["upstream_input_length_bp"]
            or _digest_bytes(sequences["REF"]) != source["reference_sequence_sha256"]
            or _digest_bytes(sequences["ALT"]) != source["alternative_sequence_sha256"]
        ):
            raise DnaLmFixtureError("upstream sequence identity differs")
        for context_id, context in contexts.items():
            length = int(context["input_length_bp"])
            sliced = {allele: _slice(sequence, length) for allele, sequence in sequences.items()}
            variant_index = _validate_alleles(
                sliced["REF"], sliced["ALT"], sliced["REF_RC"], sliced["ALT_RC"]
            )
            for allele, sequence in sliced.items():
                fasta_records[context_id].append(
                    (f"{fixture_id}|{context_id}|{allele}", sequence)
                )
            pool_length = int(config["common_lane"]["pool_length_bp"])
            pool_start = variant_index - pool_length // 2
            pool_end = pool_start + pool_length
            for model_id in context["models"]:
                manifest_rows.append(
                    {
                        "fixture_id": fixture_id,
                        "model_id": model_id,
                        "context_id": context_id,
                        "genomic_fold": source["genomic_fold"],
                        "contig": source["contig"],
                        "anchor0": source["anchor0"],
                        "variant_pos_1based": source["variant_pos_1based"],
                        "ref": source["ref"],
                        "alt": source["alt"],
                        "input_length_bp": length,
                        "forward_variant_index0": variant_index,
                        "reverse_complement_variant_index0": length // 2 - 1,
                        "pool_start0": pool_start,
                        "pool_end0": pool_end,
                        "pool_length_bp": pool_length,
                        "reference_sequence_sha256": _digest_bytes(sliced["REF"]),
                        "alternative_sequence_sha256": _digest_bytes(sliced["ALT"]),
                        "reverse_complement_reference_sha256": _digest_bytes(sliced["REF_RC"]),
                        "reverse_complement_alternative_sha256": _digest_bytes(sliced["ALT_RC"]),
                        "source_reference_sequence_sha256": source[
                            "reference_sequence_sha256"
                        ],
                        "source_alternative_sequence_sha256": source[
                            "alternative_sequence_sha256"
                        ],
                        "allele_effect_sign": "ALT_minus_REF",
                        "native_task_contract": config["native_task_contracts"][model_id],
                    }
                )

    for context_id, records in fasta_records.items():
        _write_fasta_gz(output / f"{context_id}.alleles.fa.gz", records)
    with (output / "sequence_manifest.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    firewall_receipt = {
        "schema_version": "masld-bench-dna-lm-outcome-firewall-v1",
        **firewall,
        "allowed_inputs": [
            "frozen outcome-free AlphaGenome/Sei training-fold sequence fixture",
            "frozen DNA-LM fixture configuration",
        ],
        "fixture_genomic_roles": ["train"],
        "eligible_genomic_folds": [2, 3, 4],
    }
    (output / "outcome_firewall.json").write_text(
        json.dumps(firewall_receipt, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    summary = {
        "schema_version": "masld-bench-dna-lm-sequence-fixture-v1",
        "status": "pass",
        "reference_build": config["reference_build"],
        "upstream_manifest_sha256": _digest_file(manifest_path),
        "upstream_fasta_sha256": _digest_file(fasta_path),
        "anchors": len(source_rows),
        "models": list(config["common_lane"]["models"]),
        "contexts": {
            context_id: {
                "input_length_bp": context["input_length_bp"],
                "models": context["models"],
                "fasta_records": len(fasta_records[context_id]),
            }
            for context_id, context in contexts.items()
        },
        "sequence_manifest_rows": len(manifest_rows),
        "common_pool_length_bp": config["common_lane"]["pool_length_bp"],
        "common_head_inputs": config["common_lane"]["head_inputs"],
        "common_heads": config["common_lane"]["heads"],
        "embeddings_biologically_interpretable_without_head": False,
        **firewall,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build_fixture(args.config, args.upstream, args.output)


if __name__ == "__main__":
    main()
