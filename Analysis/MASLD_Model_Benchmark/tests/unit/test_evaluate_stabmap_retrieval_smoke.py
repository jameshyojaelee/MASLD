from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_stabmap_retrieval_smoke.py"
SPEC = importlib.util.spec_from_file_location("stabmap_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class StabMapRetrievalEvaluatorTests(unittest.TestCase):
    def test_model_roster_is_candidate_and_linear_baseline(self) -> None:
        self.assertEqual(evaluator.MODELS, ("stabmap", "linear_cca"))

    def test_fold_root_parser_requires_exact_five_fold_census(self) -> None:
        roots = evaluator.parse_fold_roots([f"{fold}=/tmp/fold{fold}" for fold in range(5)])
        self.assertEqual(set(roots), set(range(5)))
        with self.assertRaises(evaluator.StabMapRetrievalEvaluationError):
            evaluator.parse_fold_roots(["0=/tmp/fold0"])

    def test_embedding_loader_requires_exact_16d_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "embedding.tsv"
            header = ["query_id"] + [f"reference_rna_PC{i}" for i in range(1, 17)]
            path.write_text(
                "\t".join(header) + "\n" + "id_a\t" + "\t".join(["1"] * 16) + "\n",
                encoding="utf-8",
            )
            identifiers, values = evaluator.load_stabmap_embedding(path)
            self.assertEqual(identifiers, ["id_a"])
            self.assertEqual(values.shape, (1, 16))


if __name__ == "__main__":
    unittest.main()
