from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "aggregate_midas_fold_predictions.py"
SPEC = importlib.util.spec_from_file_location("midas_aggregate_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
aggregate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(aggregate)


class MIDASAggregateTests(unittest.TestCase):
    def test_parser_requires_exact_fold_census(self) -> None:
        roots = aggregate.parse_fold_roots(
            [f"{fold}=/tmp/fold_{fold}" for fold in range(5)]
        )
        self.assertEqual(set(roots), set(range(5)))
        with self.assertRaises(aggregate.MIDASAggregationError):
            aggregate.parse_fold_roots(["0=/tmp/fold_0"])

    def test_aggregation_never_reads_outcomes(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn('"prediction_frozen_before_outcomes": True', source)
        self.assertIn('"outcomes_read": False', source)
        self.assertNotIn("source_h5", source)


if __name__ == "__main__":
    unittest.main()
