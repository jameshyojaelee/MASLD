from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_scbutterfly_context_transfer_smoke.py"
SPEC = importlib.util.spec_from_file_location("scbutterfly_context_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class ScButterflyContextEvaluatorTests(unittest.TestCase):
    def test_candidate_identity_is_safe_scbutterfly(self) -> None:
        self.assertEqual(evaluator.base.CANDIDATE_MODEL_ID, "scbutterfly_b_safe")
        self.assertEqual(evaluator.base.PAIRING_TOPOLOGY, "same_nucleus_training")

    def test_candidate_and_training_baselines_are_exact(self) -> None:
        self.assertEqual(
            evaluator.base.MODEL_IDS,
            (
                "scbutterfly_b_safe",
                "training_lineage_mean",
                "training_global_mean",
            ),
        )


if __name__ == "__main__":
    unittest.main()
