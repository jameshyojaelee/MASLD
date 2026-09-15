#!/usr/bin/env python3
"""Freeze the outcome-blind 230-bp fixture for task-native sequence controls."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-gse281364-task-native-sequence-controls-v1"
SEEDS = (1103, 2909, 4721, 6673, 8111)
FOLDS = tuple(f"fold-{index}" for index in range(5))
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class SequenceFixtureError(RuntimeError):
    """Raised when an outcome or sequence authority differs."""


def file_sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return sha256(sequence.encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    if set(sequence) - set("ACGTN"):
        raise SequenceFixtureError("sequence alphabet differs")
    return sequence.translate(COMPLEMENT)[::-1]


def load_config(path: Path) -> dict[str, Any]:
    value = tomllib.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SequenceFixtureError("config must be a table")
    return value


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SequenceFixtureError("JSON authority must be an object")
    return value


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return tuple(reader.fieldnames or ()), [dict(row) for row in reader]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        for line in handle:
            value = line.strip()
            if value.startswith(">"):
                if name is not None:
                    if name in records:
                        raise SequenceFixtureError("duplicate FASTA name")
                    records[name] = "".join(pieces)
                name, pieces = value[1:], []
            elif name is None:
                raise SequenceFixtureError("FASTA sequence precedes name")
            else:
                pieces.append(value.upper())
    if name is not None:
        if name in records:
            raise SequenceFixtureError("duplicate FASTA name")
        records[name] = "".join(pieces)
    return records


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != "prespecified_outcome_blind_task_and_fixture_freeze"
        or config.get("dataset_id") != "gse281364"
        or config.get("outcome_access_authorized") is not False
        or config.get("reporter_count_access_authorized") is not False
        or config.get("metric_calculation_authorized") is not False
        or config.get("model_training_authorized") is not False
        or config.get("prediction_generation_authorized") is not False
        or config.get("selection_authorized") is not False
        or config.get("promotion_authorized") is not False
    ):
        raise SequenceFixtureError("outcome and action firewall differs")
    row = config.get("row_contract", {})
    if (
        row.get("elements") != 1033
        or row.get("source_locus_groups") != 1033
        or row.get("long_range_blocks") != 239
        or row.get("outer_folds") != 5
        or tuple(row.get("fixed_seeds", ())) != SEEDS
        or row.get("seeds_are_genuinely_distinct_fits") is not True
        or row.get("schema_seed_repeats_permitted") is not False
    ):
        raise SequenceFixtureError("row contract differs")
    sequence = config.get("sequence_contract", {})
    if (
        sequence.get("input_length_bp") != 230
        or tuple(sequence.get("channel_order", ())) != ("A", "G", "C", "T")
        or sequence.get("forward_variant_index0") != 115
        or sequence.get("reverse_complement_variant_index0") != 114
        or sequence.get("allele_sign") != "ALT_minus_REF"
        or sequence.get("shared_ref_alt_encoder") is not True
        or sequence.get("observed_query_atac_allowed") is not False
        or sequence.get("local_atac_weight_reuse_allowed") is not False
    ):
        raise SequenceFixtureError("sequence contract differs")


def validate_trees(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    trees: dict[str, Path] = {}
    for label, authority in config["immutable_trees"].items():
        relative = Path(authority["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise SequenceFixtureError("unsafe tree path")
        tree = (root / relative).resolve(strict=True)
        tree.relative_to(root)
        if file_sha256(tree / "ARTIFACTS.json") != authority["artifacts_sha256"]:
            raise SequenceFixtureError(f"tree manifest differs: {label}")
        verify_frozen_tree(tree)
        trees[label] = tree
    selected = trees["reconciliation"] / config["immutable_trees"]["reconciliation"][
        "selected_member"
    ]
    if file_sha256(selected) != config["immutable_trees"]["reconciliation"][
        "selected_member_sha256"
    ]:
        raise SequenceFixtureError("selected-element member differs")
    return trees


def task_spec(config: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-outcome-blind-training-task-spec-v1",
        "status": "frozen_training_not_authorized_by_this_artifact",
        "dataset_id": "gse281364",
        "task_id": config["task_id"],
        "models": [
            config["architectures"]["sequence_cnn_control"]["implementation_id"],
            config["architectures"]["sequence_transformer_control"][
                "implementation_id"
            ],
        ],
        "row_contract": config["row_contract"],
        "sequence_contract": config["sequence_contract"],
        "target_contract": config["target_contract"],
        "split_contract": config["split_contract"],
        "optimization": config["optimization"],
        "architectures": config["architectures"],
        "artifact_contract": config["artifact_contract"],
        "claim_contract": config["claim_contract"],
        "outcomes_read": False,
        "metrics_calculated": False,
        "training_executed": False,
        "predictions_generated": False,
    }


def license_disposition(config: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-project-control-license-disposition-v1",
        "status": "internal_execution_allowed_release_blocked",
        "models": ["sequence_cnn_control", "sequence_transformer_control"],
        **config["license_disposition"],
        "checkpoint_exposure": "none_upstream_random_initialization",
        "future_project_fit_exposure": (
            "continual_seen_for_outer_training_loci_and_targets; "
            "target_label_unexposed_for_held_outer_blocks_only"
        ),
        "no_local_atac_weights": True,
        "no_pretrained_weights": True,
    }


def build(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise SequenceFixtureError("fixture output exists")
    config = load_config(arguments.config)
    validate_config(config)
    root = arguments.project_root.resolve(strict=True)
    trees = validate_trees(root, config)
    selected_fields, selected = read_tsv(
        trees["reconciliation"]
        / config["immutable_trees"]["reconciliation"]["selected_member"]
    )
    source_fields, source_rows = read_tsv(
        trees["sequence_source"] / "fixture/manifest.tsv"
    )
    if (
        len(selected) != 1033
        or len({row["element_id"] for row in selected}) != 1033
        or len(source_rows) != 4359
        or "alt_minus_ref_score" not in selected_fields
        or "reference_sequence_sha256" not in source_fields
    ):
        raise SequenceFixtureError("source row census or schema differs")
    source = {row["element_id"]: row for row in source_rows}
    records = read_fasta(
        trees["sequence_source"] / "fixture/mpralegnet.alleles.fa.gz"
    )
    if len(records) != 4 * 4359:
        raise SequenceFixtureError("source FASTA census differs")

    manifest: list[dict[str, Any]] = []
    output_sequences: dict[str, str] = {}
    identity_blocks: dict[str, set[str]] = {}
    for selected_row in sorted(selected, key=lambda row: row["element_id"]):
        element = selected_row["element_id"]
        source_row = source.get(element)
        if source_row is None:
            raise SequenceFixtureError("selected element lacks sequence source")
        fixture = source_row["fixture_id"]
        values = {
            allele: records.get(f"{fixture}|{allele}", "")
            for allele in ("REF", "ALT", "REF_RC", "ALT_RC")
        }
        ref, alt, ref_rc, alt_rc = (values[key] for key in values)
        fold = selected_row["outer_fold"]
        if (
            fold not in FOLDS
            or any(len(value) != 230 for value in values.values())
            or set("".join(values.values())) - set("ACGT")
            or ref[115] != source_row["ref"]
            or alt[115] != source_row["alt"]
            or [index for index, pair in enumerate(zip(ref, alt, strict=True)) if pair[0] != pair[1]]
            != [115]
            or ref_rc != reverse_complement(ref)
            or alt_rc != reverse_complement(alt)
            or ref_rc[114] != reverse_complement(source_row["ref"])
            or alt_rc[114] != reverse_complement(source_row["alt"])
            or sequence_sha256(ref) != source_row["reference_sequence_sha256"]
            or sequence_sha256(alt) != source_row["alternative_sequence_sha256"]
            or source_row["outer_locus_sequence_group_id"]
            != selected_row["source_locus_group_id"]
        ):
            raise SequenceFixtureError("allele, orientation, or grouping differs")
        for allele, sequence in values.items():
            name = f"{element}|{allele}"
            output_sequences[name] = sequence
            identity = min(sequence, reverse_complement(sequence))
            identity_blocks.setdefault(identity, set()).add(
                selected_row["long_range_block_id"]
            )
        manifest.append(
            {
                "element_id": element,
                "source_locus_group_id": selected_row["source_locus_group_id"],
                "long_range_block_id": selected_row["long_range_block_id"],
                "outer_fold": fold,
                "contig": source_row["contig"],
                "variant_pos0": source_row["variant_pos0"],
                "ref": source_row["ref"],
                "alt": source_row["alt"],
                "input_length_bp": 230,
                "forward_variant_index0": 115,
                "reverse_complement_variant_index0": 114,
                "reference_sequence_sha256": sequence_sha256(ref),
                "alternative_sequence_sha256": sequence_sha256(alt),
                "reference_reverse_complement_sha256": sequence_sha256(ref_rc),
                "alternative_reverse_complement_sha256": sequence_sha256(alt_rc),
                "input_lane": "GRCh38p14_230bp_genomic_window_transfer",
                "outcomes_present": "false",
            }
        )
    if any(len(blocks) != 1 for blocks in identity_blocks.values()):
        raise SequenceFixtureError("exact or reverse-complement identity crosses blocks")
    fold_elements = [sum(row["outer_fold"] == fold for row in manifest) for fold in FOLDS]
    fold_blocks = [
        len({row["long_range_block_id"] for row in manifest if row["outer_fold"] == fold})
        for fold in FOLDS
    ]
    if (
        fold_elements != config["row_contract"]["fold_element_counts"]
        or fold_blocks != config["row_contract"]["fold_block_counts"]
        or len({row["source_locus_group_id"] for row in manifest}) != 1033
        or len({row["long_range_block_id"] for row in manifest}) != 239
        or len(output_sequences) != 4132
    ):
        raise SequenceFixtureError("frozen universe census differs")

    arguments.output.mkdir(mode=0o750)
    fixture_dir = arguments.output / "fixture"
    contract_dir = arguments.output / "contract"
    fixture_dir.mkdir(mode=0o750)
    contract_dir.mkdir(mode=0o750)
    write_tsv(fixture_dir / "manifest.tsv", tuple(manifest[0]), manifest)
    fasta_buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=fasta_buffer, mtime=0) as zipped:
        for name in sorted(output_sequences):
            zipped.write(f">{name}\n{output_sequences[name]}\n".encode("ascii"))
    (fixture_dir / "sequence_controls.alleles.fa.gz").write_bytes(fasta_buffer.getvalue())
    split_rows = []
    for outer_index, outer_fold in enumerate(FOLDS):
        inner_valid = FOLDS[(outer_index + 1) % 5]
        inner_train = [fold for fold in FOLDS if fold not in {outer_fold, inner_valid}]
        refit = [fold for fold in FOLDS if fold != outer_fold]
        split_rows.append(
            {
                "outer_test_fold": outer_fold,
                "inner_validation_fold": inner_valid,
                "inner_training_folds": ";".join(inner_train),
                "final_refit_folds": ";".join(refit),
                "test_elements": fold_elements[outer_index],
                "test_long_range_blocks": fold_blocks[outer_index],
                "held_test_outcomes_visible_to_adapter": "false",
            }
        )
    write_tsv(fixture_dir / "split_plan.tsv", tuple(split_rows[0]), split_rows)
    smoke = []
    for fold in FOLDS:
        candidates = [row for row in manifest if row["outer_fold"] == fold]
        smoke.extend(
            sorted(candidates, key=lambda row: sha256(row["element_id"].encode()).hexdigest())[:2]
        )
    write_tsv(fixture_dir / "smoke_manifest.tsv", tuple(manifest[0]), smoke)
    (contract_dir / "task_spec.json").write_text(
        json.dumps(task_spec(config), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (contract_dir / "code_license_disposition.json").write_text(
        json.dumps(license_disposition(config), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_outcome_blind_task_fixture_freeze",
        "dataset_id": "gse281364",
        "models": ["sequence_cnn_control", "sequence_transformer_control"],
        "elements": 1033,
        "source_locus_groups": 1033,
        "long_range_blocks": 239,
        "outer_folds": 5,
        "fold_element_counts": fold_elements,
        "fold_block_counts": fold_blocks,
        "allele_sequences": 4132,
        "smoke_elements": 10,
        "input_length_bp": 230,
        "shared_ref_alt_encoder": True,
        "reverse_complement_fixture": True,
        "identity_cross_block_violations": 0,
        "genuine_training_seeds": list(SEEDS),
        "schema_seed_repeats_permitted": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "metrics_calculated": False,
        "model_training_executed": False,
        "predictions_generated": False,
        "local_atac_weights_read": False,
    }
    (fixture_dir / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
