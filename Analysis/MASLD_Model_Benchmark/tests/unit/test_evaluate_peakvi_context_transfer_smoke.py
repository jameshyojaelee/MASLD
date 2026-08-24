from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_peakvi_context_transfer_smoke.py"
SPEC = importlib.util.spec_from_file_location("peakvi_context_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class PeakVIContextEvaluatorTests(unittest.TestCase):
    def test_candidate_identity_is_peakvi_context(self) -> None:
        self.assertEqual(evaluator.base.CANDIDATE_MODEL_ID, "peakvi_training_context")
        self.assertEqual(evaluator.base.PAIRING_TOPOLOGY, "same_nucleus_training")

    def test_candidate_and_training_baselines_are_exact(self) -> None:
        self.assertEqual(
            evaluator.base.MODEL_IDS,
            (
                "peakvi_training_context",
                "training_lineage_mean",
                "training_global_mean",
            ),
        )


if __name__ == "__main__":
    unittest.main()
