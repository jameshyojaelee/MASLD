#!/usr/bin/env python3
"""Inventory released EpiBERT checkpoint shapes without constructing a model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence


class EpiBERTShapeAuditError(ValueError):
    """Raised when released checkpoint shape inclusion differs."""


def audit_one(tf: Any, prefix: Path) -> dict[str, Any]:
    if not prefix.with_suffix(".index").is_file():
        raise EpiBERTShapeAuditError(f"checkpoint index is absent: {prefix}")
    variables = tf.train.list_variables(str(prefix))
    if not variables:
        raise EpiBERTShapeAuditError("checkpoint variable inventory is empty")
    shape_counts: dict[str, int] = {}
    attention_candidates: list[dict[str, Any]] = []
    for name, shape in variables:
        shape_key = "x".join(str(value) for value in shape) if shape else "scalar"
        shape_counts[shape_key] = shape_counts.get(shape_key, 0) + 1
        if len(shape) == 3 and shape[0] == 1024 and shape[1] in {4, 8}:
            attention_candidates.append({"name": name, "shape": list(shape)})
    head_widths = sorted(
        {(item["shape"][1], item["shape"][2]) for item in attention_candidates}
    )
    if not attention_candidates or not head_widths:
        raise EpiBERTShapeAuditError("attention kernel candidates are absent")
    return {
        "checkpoint_prefix": str(prefix),
        "variable_count": len(variables),
        "shape_counts": dict(sorted(shape_counts.items())),
        "attention_kernel_candidates": attention_candidates,
        "inferred_head_width_pairs": [list(value) for value in head_widths],
    }


def run(checkpoints: Sequence[tuple[str, Path]], output: Path) -> dict[str, Any]:
    import tensorflow as tf

    if output.exists() or len(checkpoints) != 3:
        raise EpiBERTShapeAuditError("output or checkpoint roster differs")
    names = [name for name, _path in checkpoints]
    if names != ["pretrained_model1", "pretrained_model2", "fine_tuned_rampage"]:
        raise EpiBERTShapeAuditError("checkpoint names differ")
    records = {name: audit_one(tf, path) for name, path in checkpoints}
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-epibert-checkpoint-shape-audit-v1",
        "status": "pass",
        "tensorflow": tf.__version__,
        "models": records,
        "published_pretraining_command_num_heads": 4,
        "published_finetuning_command_num_heads": 8,
        "constructor_changed_by_audit": False,
        "project_data_read": False,
        "outcomes_read": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pretrained-model1", type=Path, required=True)
    parser.add_argument("--pretrained-model2", type=Path, required=True)
    parser.add_argument("--fine-tuned-rampage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run(
        (
            ("pretrained_model1", arguments.pretrained_model1),
            ("pretrained_model2", arguments.pretrained_model2),
            ("fine_tuned_rampage", arguments.fine_tuned_rampage),
        ),
        arguments.output,
    )


if __name__ == "__main__":
    main()
