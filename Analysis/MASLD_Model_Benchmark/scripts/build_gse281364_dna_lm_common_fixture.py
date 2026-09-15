#!/usr/bin/env python3
"""Derive shared outcome-blind DNA-LM inputs from the frozen GSE281364 Sei fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Iterable, Mapping


EXPECTED_ELEMENTS = 1_033
EXPECTED_SEQUENCES = EXPECTED_ELEMENTS * 4
EXPECTED_FOLDS = 5
EXPECTED_MODELS = ("dnabert2", "nucleotide_transformer", "hyenadna")
EXPECTED_ALLELES = ("REF", "ALT", "REF_RC", "ALT_RC")
BASES = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MANIFEST_FIELDS = (
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
    "fasta_record_prefix",
)
MODEL_FIELDS = (
    "model_id",
    "context_id",
    "input_length_bp",
    "input_adapter",
    "pool_start0",
    "pool_end0",
    "head_feature_blocks",
    "downstream_heads",
    "head_training",
)


class CommonFixtureError(ValueError):
    """Raised when the shared GSE281364 DNA-LM fixture requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def sequence_digest(sequence: str) -> str:
    return sha256(sequence.encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CommonFixtureError("TSV has no header")
        return [dict(row) for row in reader]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                if not name or name in records:
                    raise CommonFixtureError("upstream FASTA header differs")
                records[name] = []
            elif name is None:
                raise CommonFixtureError("upstream FASTA sequence precedes header")
            else:
                records[name].append(line)
    return {record_id: "".join(parts) for record_id, parts in records.items()}


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="ascii", newline="\n") as handle:
                for name, sequence in records:
                    handle.write(f">{name}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80] + "\n")


def write_tsv(
    path: Path, fields: Iterable[str], rows: Iterable[Mapping[str, object]]
) -> None:
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


def validate_pair(
    row: Mapping[str, str],
    sequences: Mapping[str, str],
    *,
    input_length: int,
    forward_index: int,
    reverse_index: int,
) -> None:
    if set(sequences) != set(EXPECTED_ALLELES):
        raise CommonFixtureError("allele sequence set differs")
    ref, alt = sequences["REF"], sequences["ALT"]
    ref_rc, alt_rc = sequences["REF_RC"], sequences["ALT_RC"]
    if (
        {len(sequence) for sequence in sequences.values()} != {input_length}
        or any(set(sequence) - BASES for sequence in sequences.values())
        or ref_rc != reverse_complement(ref)
        or alt_rc != reverse_complement(alt)
    ):
        raise CommonFixtureError("allele or reverse-complement sequence differs")
    differences = [index for index, pair in enumerate(zip(ref, alt)) if pair[0] != pair[1]]
    reverse_differences = [
        index for index, pair in enumerate(zip(ref_rc, alt_rc)) if pair[0] != pair[1]
    ]
    if (
        differences != [forward_index]
        or reverse_differences != [reverse_index]
        or ref[forward_index] != row["ref"]
        or alt[forward_index] != row["alt"]
        or ref_rc[reverse_index] != reverse_complement(row["ref"])
        or alt_rc[reverse_index] != reverse_complement(row["alt"])
    ):
        raise CommonFixtureError("REF/ALT variant geometry differs")
    expected_hashes = {
        "REF": row["reference_sequence_sha256"],
        "ALT": row["alternative_sequence_sha256"],
        "REF_RC": row["reference_reverse_complement_sha256"],
        "ALT_RC": row["alternative_reverse_complement_sha256"],
    }
    if any(sequence_digest(sequences[key]) != value for key, value in expected_hashes.items()):
        raise CommonFixtureError("upstream sequence digest differs")


def validate_config(config: Mapping[str, object]) -> dict[str, object]:
    if (
        config.get("schema_version")
        != "masld-bench-gse281364-dna-lm-common-fixture-v1"
        or config.get("dataset_id") != "gse281364"
    ):
        raise CommonFixtureError("fixture configuration differs")
    lane = config.get("common_lane")
    firewall = config.get("firewall")
    adapters = config.get("model_input_adapters")
    if (
        not isinstance(lane, dict)
        or lane.get("context_id") != "common_4096"
        or lane.get("input_length_bp") != 4_096
        or lane.get("forward_variant_index0") != 2_048
        or lane.get("reverse_complement_variant_index0") != 2_047
        or lane.get("pool_start0") != 1_280
        or lane.get("pool_end0") != 2_816
        or lane.get("pool_length_bp") != 1_536
        or tuple(lane.get("models", [])) != EXPECTED_MODELS
        or tuple(lane.get("alleles", [])) != EXPECTED_ALLELES
        or not isinstance(adapters, dict)
        or tuple(adapters) != EXPECTED_MODELS
        or not isinstance(firewall, dict)
        or any(firewall.values())
    ):
        raise CommonFixtureError("common lane or outcome firewall differs")
    return lane


def build(
    *,
    config_path: Path,
    split_root: Path,
    sei_root: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists() or any(path.is_symlink() for path in (config_path, split_root, sei_root)):
        raise CommonFixtureError("fixture request differs")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    lane = validate_config(config)
    split_sha256 = str(config["split_artifacts_sha256"])
    sei_sha256 = str(config["sei_fixture_artifacts_sha256"])
    if (
        digest(split_root / "ARTIFACTS.json") != split_sha256
        or digest(sei_root / "ARTIFACTS.json") != sei_sha256
    ):
        raise CommonFixtureError("frozen upstream artifact changed")
    split_receipt = json.loads((split_root / "split/receipt.json").read_text())
    sei_receipt = json.loads((sei_root / "fixture/receipt.json").read_text())
    if (
        split_receipt.get("status") != "pass_outcome_blind_split_contract"
        or split_receipt.get("outer_locus_sequence_groups") != EXPECTED_ELEMENTS
        or split_receipt.get("outer_folds") != EXPECTED_FOLDS
        or any(
            split_receipt.get(field) is not False
            for field in ("outcomes_read", "reporter_counts_read", "sealed_outcomes_read")
        )
        or sei_receipt.get("status") != "pass_outcome_blind_sei_fixture"
        or sei_receipt.get("elements") != EXPECTED_ELEMENTS
        or sei_receipt.get("sequences") != EXPECTED_SEQUENCES
        or sei_receipt.get("split_artifacts_sha256") != split_sha256
        or any(
            sei_receipt.get(field) is not False
            for field in ("outcomes_read", "reporter_counts_read", "sealed_outcomes_read")
        )
    ):
        raise CommonFixtureError("outcome-blind upstream receipt differs")

    group_rows = read_tsv(split_root / "split/groups.tsv")
    group_folds = {
        row["outer_locus_sequence_group_id"]: row["outer_fold"] for row in group_rows
    }
    if len(group_rows) != EXPECTED_ELEMENTS or len(group_folds) != EXPECTED_ELEMENTS:
        raise CommonFixtureError("split group census differs")
    source_rows = read_tsv(sei_root / "fixture/manifest.tsv")
    source_fasta = read_fasta(sei_root / "fixture/sei.alleles.fa.gz")
    if len(source_rows) != EXPECTED_ELEMENTS or len(source_fasta) != EXPECTED_SEQUENCES:
        raise CommonFixtureError("Sei fixture census differs")

    output.mkdir(parents=True)
    manifest_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []
    seen_groups: set[str] = set()
    input_length = int(lane["input_length_bp"])
    forward_index = int(lane["forward_variant_index0"])
    reverse_index = int(lane["reverse_complement_variant_index0"])
    context_id = str(lane["context_id"])
    for source in source_rows:
        group = source["outer_locus_sequence_group_id"]
        if (
            group in seen_groups
            or group_folds.get(group) != source["outer_fold"]
            or int(source["variant_index0"]) != forward_index
            or int(source["input_end0"]) - int(source["input_start0"]) != input_length
        ):
            raise CommonFixtureError("locus-group fold or geometry differs")
        seen_groups.add(group)
        source_id = source["fixture_id"]
        sequences = {
            allele: source_fasta[f"{source_id}|{allele}"] for allele in EXPECTED_ALLELES
        }
        validate_pair(
            source,
            sequences,
            input_length=input_length,
            forward_index=forward_index,
            reverse_index=reverse_index,
        )
        fixture_id = "dna_lm_" + sha256(source["element_id"].encode()).hexdigest()[:24]
        prefix = f"{fixture_id}|{context_id}"
        fasta_records.extend(
            (f"{prefix}|{allele}", sequences[allele]) for allele in EXPECTED_ALLELES
        )
        manifest_rows.append(
            {
                "fixture_id": fixture_id,
                "source_fixture_id": source_id,
                "element_id": source["element_id"],
                "outer_locus_sequence_group_id": group,
                "outer_fold": source["outer_fold"],
                "contig": source["contig"],
                "variant_pos0": source["variant_pos0"],
                "variant_pos1": source["variant_pos1"],
                "input_start0": source["input_start0"],
                "input_end0": source["input_end0"],
                "input_length_bp": input_length,
                "forward_variant_index0": forward_index,
                "reverse_complement_variant_index0": reverse_index,
                "pool_start0": lane["pool_start0"],
                "pool_end0": lane["pool_end0"],
                "pool_length_bp": lane["pool_length_bp"],
                "ref": source["ref"],
                "alt": source["alt"],
                "reference_sequence_sha256": source["reference_sequence_sha256"],
                "alternative_sequence_sha256": source["alternative_sequence_sha256"],
                "reverse_complement_reference_sha256": source[
                    "reference_reverse_complement_sha256"
                ],
                "reverse_complement_alternative_sha256": source[
                    "alternative_reverse_complement_sha256"
                ],
                "fasta_record_prefix": prefix,
            }
        )
    if seen_groups != set(group_folds):
        raise CommonFixtureError("Sei fixture does not cover every frozen locus group")

    write_fasta(output / f"{context_id}.alleles.fa.gz", fasta_records)
    write_tsv(output / "sequence_manifest.tsv", MANIFEST_FIELDS, manifest_rows)
    model_rows = [
        {
            "model_id": model,
            "context_id": context_id,
            "input_length_bp": input_length,
            "input_adapter": config["model_input_adapters"][model],
            "pool_start0": lane["pool_start0"],
            "pool_end0": lane["pool_end0"],
            "head_feature_blocks": ",".join(lane["head_feature_blocks"]),
            "downstream_heads": ",".join(lane["identical_downstream_heads"]),
            "head_training": lane["head_training"],
        }
        for model in EXPECTED_MODELS
    ]
    write_tsv(output / "model_contracts.tsv", MODEL_FIELDS, model_rows)
    firewall = {
        "schema_version": "masld-bench-gse281364-dna-lm-outcome-firewall-v1",
        **config["firewall"],
        "allowed_inputs": [
            "frozen GSE281364 outcome-blind locus-group split receipt and groups",
            "frozen GSE281364 outcome-blind Sei sequence fixture",
            "frozen common fixture configuration",
        ],
    }
    (output / "outcome_firewall.json").write_text(
        json.dumps(firewall, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    fold_elements = {
        str(fold): sum(int(row["outer_fold"]) == fold for row in manifest_rows)
        for fold in range(EXPECTED_FOLDS)
    }
    summary: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-dna-lm-common-sequence-fixture-v1",
        "status": "pass_outcome_blind_common_fixture",
        "dataset_id": "gse281364",
        "context_id": context_id,
        "models": list(EXPECTED_MODELS),
        "elements": len(manifest_rows),
        "outer_locus_sequence_groups": len(seen_groups),
        "outer_folds": EXPECTED_FOLDS,
        "fold_elements": fold_elements,
        "sequences": len(fasta_records),
        "input_length_bp": input_length,
        "forward_variant_index0": forward_index,
        "reverse_complement_variant_index0": reverse_index,
        "pool_start0": lane["pool_start0"],
        "pool_end0": lane["pool_end0"],
        "pool_length_bp": lane["pool_length_bp"],
        "head_feature_blocks": lane["head_feature_blocks"],
        "identical_downstream_heads": lane["identical_downstream_heads"],
        "split_artifacts_sha256": split_sha256,
        "sei_fixture_artifacts_sha256": sei_sha256,
        "upstream_manifest_sha256": digest(sei_root / "fixture/manifest.tsv"),
        "upstream_fasta_sha256": digest(sei_root / "fixture/sei.alleles.fa.gz"),
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "checkpoint_bytes_loaded": False,
        "model_forward_executed": False,
        "downstream_head_fit": False,
        "champion_eligible": False,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--sei-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(build(
        config_path=arguments.config,
        split_root=arguments.split_root,
        sei_root=arguments.sei_root,
        output=arguments.output,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
