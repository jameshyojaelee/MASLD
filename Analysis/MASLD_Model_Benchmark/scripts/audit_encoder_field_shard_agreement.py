#!/usr/bin/env python3
"""Compare this experiment's shards with the already-frozen ones they overlap.

Geneformer v2 316M and TranscriptFormer TF-Sapiens were already fit through the
common-head producer on the same embedding bundles, the same split, the same
three screen seeds and the same five outer folds.  If this experiment's head is
genuinely the same head, its shards for those two blocks must agree with the
frozen ones to float noise.  Sixty shard tables therefore act as an independent
reproduction control on the whole loop -- not a digest comparison, an actual
re-derivation of every predicted probability.

Disagreement here means the head differs, and every delta in the panel that
compares a newly fit block against those two is then not like-for-like.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for entry in (str(ROOT), str(ROOT / "src")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from scripts.fit_predict_common_cell_heads_study_50000 import (  # noqa: E402
    HEADS,
    OUTER_FOLDS,
    ROSTER,
    SCREEN_SEEDS,
)

FROZEN_BUNDLES = {
    "geneformer_v2_316m": "executions/model-training-21100295/shards",
    "transcriptformer_tf_sapiens": "executions/model-training-21082720/shards",
}
TOLERANCE = 1.0e-9


def read_table(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    values = np.asarray(
        [[float(row[f"probability::{label}"]) for label in ROSTER] for row in rows],
        dtype=np.float64,
    )
    return [row["row_id"] for row in rows], values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    records: list[dict[str, Any]] = []
    for block_id, frozen_relative in FROZEN_BUNDLES.items():
        frozen_root = ROOT / frozen_relative
        if not frozen_root.is_dir():
            records.append(
                {
                    "block_id": block_id,
                    "status": "FROZEN_BUNDLE_ABSENT",
                    "frozen_root": frozen_relative,
                }
            )
            continue
        differences: list[float] = []
        compared = 0
        misaligned: list[str] = []
        for head_id in HEADS:
            for seed in SCREEN_SEEDS:
                for fold in OUTER_FOLDS:
                    frozen_path = (
                        frozen_root / f"{head_id}__seed{seed}__fold{fold}" / "predictions.tsv"
                    )
                    produced_path = (
                        arguments.shards
                        / f"{block_id}__{head_id}__seed{seed}__fold{fold}"
                        / "predictions.tsv"
                    )
                    if not frozen_path.is_file() or not produced_path.is_file():
                        continue
                    frozen_ids, frozen = read_table(frozen_path)
                    produced_ids, produced = read_table(produced_path)
                    if frozen_ids != produced_ids:
                        misaligned.append(frozen_path.parent.name)
                        continue
                    differences.append(float(np.max(np.abs(frozen - produced))))
                    compared += 1
        max_difference = max(differences) if differences else None
        records.append(
            {
                "block_id": block_id,
                "frozen_root": frozen_relative,
                "shards_compared": compared,
                "row_misaligned_shards": misaligned,
                "max_absolute_probability_difference": max_difference,
                "tolerance": TOLERANCE,
                "status": (
                    "verified_positively"
                    if compared == len(HEADS) * len(SCREEN_SEEDS) * len(OUTER_FOLDS)
                    and not misaligned
                    and max_difference is not None
                    and max_difference < TOLERANCE
                    else "REPRODUCTION_FAILED"
                ),
            }
        )
    summary = {
        "schema_version": "masld-bench-encoder-field-shard-agreement-v1",
        "check": "independent_reproduction_of_already_frozen_common_head_shards",
        "blocks": records,
        "all_reproduced": all(
            record["status"] == "verified_positively" for record in records
        ),
    }
    arguments.output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["all_reproduced"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
