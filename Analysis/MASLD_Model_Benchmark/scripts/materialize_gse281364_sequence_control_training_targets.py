#!/usr/bin/env python3
"""Materialize one outer-fold training-only target view for sequence controls."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-gse281364-task-native-sequence-production-v1"
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
ALLELES = ("ref", "alt")
FOLDS = tuple(f"fold-{index}" for index in range(5))


class TargetMaterializationError(RuntimeError):
    """Raised when target inclusion or held-fold isolation differs."""


def file_sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    value = tomllib.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TargetMaterializationError("campaign config must be a table")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_campaign(config: Mapping[str, Any]) -> None:
    execution = config.get("execution", {})
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status")
        != "authorized_prediction_only_development_campaign"
        or config.get("dataset_id") != "gse281364"
        or config.get("outcome_role") != "exposed_development_MPRA_only"
        or config.get("model_training_authorized") is not True
        or config.get("prediction_generation_authorized") is not True
        or config.get("adapter_metric_calculation_authorized") is not False
        or config.get("adapter_scoring_authorized") is not False
        or config.get("promotion_authorized") is not False
        or tuple(execution.get("outer_folds", ())) != FOLDS
        or execution.get("held_test_outcomes_visible_to_training_adapter")
        is not False
        or execution.get("target_materialization_process_separate_from_training_adapter")
        is not True
        or execution.get("private_target_views_retained_in_frozen_output")
        is not False
    ):
        raise TargetMaterializationError("campaign authorization differs")


def resolve_tree(root: Path, authority: Mapping[str, Any], label: str) -> Path:
    relative = Path(authority["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise TargetMaterializationError(f"unsafe {label} tree path")
    tree = (root / relative).resolve(strict=True)
    tree.relative_to(root)
    if file_sha256(tree / "ARTIFACTS.json") != authority["artifacts_sha256"]:
        raise TargetMaterializationError(f"{label} tree hash differs")
    verify_frozen_tree(tree)
    return tree


def aggregate_training_targets(
    outcome_path: Path,
    manifest: Sequence[Mapping[str, str]],
    outer_fold: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    fold_by_element = {row["element_id"]: row["outer_fold"] for row in manifest}
    if (
        outer_fold not in FOLDS
        or len(fold_by_element) != 1033
        or set(fold_by_element.values()) != set(FOLDS)
    ):
        raise TargetMaterializationError("manifest or outer fold differs")
    training = {
        element for element, fold in fold_by_element.items() if fold != outer_fold
    }
    held = set(fold_by_element) - training
    values: dict[tuple[str, str, str], dict[int, float]] = defaultdict(dict)
    held_rows_seen = 0
    with gzip.open(
        outcome_path, "rt", encoding="utf-8", errors="strict", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "element_id",
            "allele",
            "context_id",
            "experimental_replicate",
            "sample_id",
            "DNA",
            "RNA",
            "assay_state",
            "missing_reason",
            "pairing",
            "biological_unit",
            "donor_id",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise TargetMaterializationError("outcome schema differs")
        for row in reader:
            element = row["element_id"]
            if element not in fold_by_element or row["context_id"] not in CONTEXTS:
                continue
            if element in held:
                held_rows_seen += 1
                continue
            allele = row["allele"]
            if allele not in ALLELES:
                raise TargetMaterializationError("selected element allele differs")
            try:
                replicate = int(row["experimental_replicate"])
                dna = float(row["DNA"])
                rna = float(row["RNA"])
            except ValueError as error:
                raise TargetMaterializationError(
                    "selected training count is not numeric"
                ) from error
            expected_sample = f"{row['context_id']}_r{replicate}"
            if (
                replicate not in {1, 2, 3, 4}
                or row["sample_id"] != expected_sample
                or row["assay_state"] != "observed"
                or row["missing_reason"] != "not_applicable"
                or row["pairing"] != "same_sample_different_aliquot"
                or row["biological_unit"] != "experimental_replicate"
                or row["donor_id"] != "not_applicable"
                or not math.isfinite(dna)
                or not math.isfinite(rna)
                or dna < 0
                or rna < 0
            ):
                raise TargetMaterializationError(
                    "selected training replicate contract differs"
                )
            key = (element, allele, row["context_id"])
            if replicate in values[key]:
                raise TargetMaterializationError("duplicate training replicate")
            values[key][replicate] = math.log2((rna + 0.5) / (dna + 0.5))
    rows: list[dict[str, Any]] = []
    for element in sorted(training):
        for context in CONTEXTS:
            allele_means: dict[str, float] = {}
            for allele in ALLELES:
                replicates = values.get((element, allele, context), {})
                if set(replicates) != {1, 2, 3, 4}:
                    raise TargetMaterializationError(
                        "training target lacks four complete replicates"
                    )
                allele_means[allele] = sum(replicates.values()) / 4.0
            rows.append(
                {
                    "element_id": element,
                    "context_id": context,
                    "reference_activity": format(allele_means["ref"], ".17g"),
                    "alternative_activity": format(allele_means["alt"], ".17g"),
                }
            )
    expected_held_rows = len(held) * len(CONTEXTS) * len(ALLELES) * 4
    if (
        len(rows) != len(training) * len(CONTEXTS)
        or len(values) != len(training) * len(CONTEXTS) * len(ALLELES)
        or held_rows_seen != expected_held_rows
        or any(row["element_id"] in held for row in rows)
    ):
        raise TargetMaterializationError("training-only target census differs")
    return rows, {
        "outer_test_fold": outer_fold,
        "training_elements": len(training),
        "held_test_elements": len(held),
        "training_target_rows": len(rows),
        "training_replicate_rows_used": len(values) * 4,
        "held_test_replicate_rows_skipped_before_count_parsing": held_rows_seen,
        "held_test_target_rows_emitted": 0,
    }


def materialize(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise TargetMaterializationError("target-view output exists")
    config = load_config(arguments.config)
    validate_campaign(config)
    root = arguments.project_root.resolve(strict=True)
    task_tree = resolve_tree(root, config["frozen_task"], "frozen task")
    outcome_tree = resolve_tree(root, config["outcome_authority"], "outcome")
    task_member = task_tree / config["frozen_task"]["task_spec_member"]
    manifest_member = task_tree / config["frozen_task"]["fixture_manifest_member"]
    outcome_member = outcome_tree / config["outcome_authority"]["member"]
    for path, expected, label in (
        (task_member, config["frozen_task"]["task_spec_sha256"], "task spec"),
        (
            manifest_member,
            config["frozen_task"]["fixture_manifest_sha256"],
            "fixture manifest",
        ),
        (
            outcome_member,
            config["outcome_authority"]["member_sha256"],
            "outcome member",
        ),
    ):
        if file_sha256(path) != expected:
            raise TargetMaterializationError(f"{label} hash differs")
    task = json.loads(task_member.read_text(encoding="utf-8"))
    if (
        task.get("status") != "frozen_training_not_authorized_by_this_artifact"
        or task.get("outcomes_read")
        or task.get("metrics_calculated")
        or task.get("training_executed")
        or task.get("predictions_generated")
    ):
        raise TargetMaterializationError("frozen TaskSpec differs")
    rows, census = aggregate_training_targets(
        outcome_member, read_tsv(manifest_member), arguments.outer_fold
    )
    arguments.output.mkdir(mode=0o700)
    target_path = arguments.output / "training_targets.tsv"
    write_tsv(
        target_path,
        (
            "element_id",
            "context_id",
            "reference_activity",
            "alternative_activity",
        ),
        rows,
    )
    receipt = {
        "schema_version": "masld-bench-training-only-target-view-v1",
        "status": "pass_outer_training_targets_only",
        "dataset_id": "gse281364",
        **census,
        "contexts": list(CONTEXTS),
        "replicates_per_allele_context": 4,
        "activity_transform": "mean_log2_RNA_plus_0.5_over_DNA_plus_0.5",
        "target_view_sha256": file_sha256(target_path),
        "held_test_outcomes_visible_to_training_adapter": False,
        "benchmark_metrics_calculated": False,
        "raw_target_view_intended_for_frozen_output": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outer-fold", choices=FOLDS, required=True)
    parser.add_argument("--output", type=Path, required=True)
    materialize(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
