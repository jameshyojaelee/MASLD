"""Tests for the development-only conditional-model evidence producer."""

from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import unittest

from masld_bench.stacking import (
    STACK_SEEDS,
    StackingError,
    _derive_stack,
    _read_component_rows,
)
from masld_bench.tournament import VARIANT_TASK


def _hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _variant_rows() -> list[dict[str, str]]:
    rows = []
    for seed in STACK_SEEDS:
        for fold in range(4):
            for offset in range(2):
                value = float(fold * 2 + offset + 1)
                rows.append(
                    {
                        "seed": str(seed),
                        "row_hash": _hash(f"row-{fold}-{offset}"),
                        "unit_hash": _hash(f"unit-{fold}-{offset}"),
                        "block_hash": _hash(f"block-{fold}"),
                        "stratum": "all",
                        "outer_fold": f"fold-{fold}",
                        "study_id": "gse281364",
                        "observed": repr(value),
                        "sequence_prediction": repr(value),
                        "context_prediction": repr(-value),
                    }
                )
    return rows


class StackingTests(unittest.TestCase):
    def test_outer_fold_fit_is_nonnegative_and_uses_other_blocks(self) -> None:
        derived = _derive_stack(_variant_rows(), VARIANT_TASK)
        self.assertEqual(len(derived["weight_rows"]), 4 * len(STACK_SEEDS))
        for row in derived["weight_rows"]:
            self.assertAlmostEqual(float(row["sequence_weight"]), 1.0)
            self.assertAlmostEqual(
                float(row["sequence_weight"]) + float(row["context_weight"]),
                1.0,
            )
            self.assertNotEqual(
                row["fit_row_set_sha256"], row["held_row_set_sha256"]
            )
        self.assertTrue(
            all(
                float(row["stack_prediction"])
                == float(row["sequence_prediction"])
                for row in derived["stack_rows"]
            )
        )

    def test_component_table_requires_all_five_fixed_seeds(self) -> None:
        rows = [row for row in _variant_rows() if row["seed"] != str(STACK_SEEDS[-1])]
        fields = tuple(rows[0])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "components.tsv"
            lines = ["\t".join(fields)]
            lines.extend("\t".join(row[field] for field in fields) for row in rows)
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(StackingError, "exactly the five fixed seeds"):
                _read_component_rows(path, VARIANT_TASK)


if __name__ == "__main__":
    unittest.main()
