from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_scglue_retrieval_smoke.py"
SPEC = importlib.util.spec_from_file_location("scglue_retrieval_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class SCGLUERetrievalEvaluatorTests(unittest.TestCase):
    def test_perfect_directional_retrieval_has_rank_one(self) -> None:
        ids = ["a", "b", "c"]
        embeddings = np.eye(3)
        ranks = evaluator.directional_ranks(
            ids, ids, embeddings, embeddings, {value: value for value in ids}
        )
        self.assertEqual(ranks, {"a": 1, "b": 1, "c": 1})

    def test_normalization_rejects_zero_embedding(self) -> None:
        with self.assertRaisesRegex(
            evaluator.SCGLUERetrievalEvaluationError, "embedding rows"
        ):
            evaluator.normalize_rows(np.zeros((2, 3)))

    def test_model_roster_is_exact(self) -> None:
        self.assertEqual(evaluator.MODELS, ("scglue", "paired_scglue", "linear_cca"))


if __name__ == "__main__":
    unittest.main()
