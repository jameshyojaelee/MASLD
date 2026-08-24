from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "midas_fit_predict_fold.py"
SPEC = importlib.util.spec_from_file_location("midas_fit_predict_fold_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
fit_predict = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fit_predict)


class MIDASFitPredictFoldTests(unittest.TestCase):
    def test_chromosome_chunks_require_contiguous_blocks(self) -> None:
        self.assertEqual(
            fit_predict.chromosome_chunks(
                [
                    {"chromosome": "chr1"},
                    {"chromosome": "chr1"},
                    {"chromosome": "chr2"},
                ]
            ),
            [2, 1],
        )
        with self.assertRaises(fit_predict.MIDASFitError):
            fit_predict.chromosome_chunks(
                [
                    {"chromosome": "chr1"},
                    {"chromosome": "chr2"},
                    {"chromosome": "chr1"},
                ]
            )

    def test_query_is_excluded_from_training_contract(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn('model.combs != [["rna", "atac"]]', source)
        self.assertIn('"query_rna_used_for_training": False', source)
        self.assertIn('"held_atac_input_exposed": False', source)

    def test_only_tensor_state_is_loaded(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn("weights_only=True", source)
        self.assertNotIn("load_from_checkpoint", source)
        self.assertNotIn("MIDAS.load(", source)


if __name__ == "__main__":
    unittest.main()
