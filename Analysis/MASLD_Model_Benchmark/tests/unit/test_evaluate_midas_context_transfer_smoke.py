from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_midas_context_transfer_smoke.py"
SPEC = importlib.util.spec_from_file_location("midas_context_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class MIDASContextEvaluatorTests(unittest.TestCase):
    def test_candidate_and_baseline_roster_is_exact(self) -> None:
        self.assertEqual(evaluator.base.CANDIDATE_MODEL_ID, "midas_inductive")
        self.assertEqual(evaluator.base.PAIRING_TOPOLOGY, "same_nucleus_training")
        self.assertEqual(
            evaluator.base.MODEL_IDS,
            (
                "midas_inductive",
                "training_lineage_mean",
                "training_global_mean",
            ),
        )


if __name__ == "__main__":
    unittest.main()
