from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_seurat_wnn_bridge_retrieval_smoke.py"
SPEC = importlib.util.spec_from_file_location("seurat_wnn_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class SeuratWNNBridgeEvaluatorTests(unittest.TestCase):
    def test_load_wnn_embedding_preserves_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "embedding.tsv"
            fields = ["query_id"] + [f"bridge-lap-{index}" for index in range(1, 17)]
            path.write_text(
                "\t".join(fields) + "\n" + "x\t" + "\t".join(["1"] * 16) + "\n",
                encoding="utf-8",
            )
            identifiers, values = evaluator.load_wnn_embedding(path)
            self.assertEqual(identifiers, ["x"])
            self.assertEqual(values.shape, (1, 16))

    def test_model_roster_is_exact(self) -> None:
        self.assertEqual(evaluator.MODELS, ("seurat_wnn_bridge", "linear_cca"))


if __name__ == "__main__":
    unittest.main()
