from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/preflight_geneformer_v2_316m_common.py"
SPEC = importlib.util.spec_from_file_location("geneformer_v2_316m_preflight", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

TOKENS = (
    ROOT
    / "executions/geneformer-sweep/21015001/"
    "geneformer_v2_316m__two_layer_mlp/prepare/rank_value_tokens.jsonl"
)
ROW_ORDER = (
    ROOT / "executions/embeddings/geneformer_v2_316m/embedding_row_order.txt"
)
EMBEDDINGS = ROOT / "executions/embeddings/geneformer_v2_316m/embeddings.npy"


class GeneformerV2316MPreflightTest(unittest.TestCase):
    def test_prior_smoke_raw_embedding_contract(self) -> None:
        observed = MODULE.validate_prior_smoke(TOKENS, ROW_ORDER, EMBEDDINGS)
        self.assertEqual(observed["rows"], 1000)
        self.assertEqual(observed["embedding_shape"], [1000, 1152])
        self.assertEqual(observed["embedding_dtype"], "float32")
        self.assertFalse(observed["labels_read_by_extraction"])
        self.assertEqual(observed["row_order"], sorted(observed["row_order"]))
        self.assertEqual(set(observed["tokens"]), set(observed["row_order"]))

    def test_parser_has_no_label_or_outcome_surface(self) -> None:
        option_strings = {
            option
            for action in MODULE.parser()._actions
            for option in action.option_strings
        }
        self.assertNotIn("--labels", option_strings)
        self.assertNotIn("--outcomes", option_strings)
        self.assertNotIn("--cell-labels", option_strings)
        self.assertIn("--prior-embeddings", option_strings)
        self.assertIn("--output-embeddings", option_strings)

    def test_frozen_raw_embedding_identity(self) -> None:
        self.assertEqual(MODULE.MODEL_ID, "geneformer_v2_316m")
        self.assertEqual(MODULE.HIDDEN_SIZE, 1152)
        self.assertEqual(MODULE.GPU_BATCH_SIZE, 32)
        self.assertEqual(
            MODULE.EMBEDDING_POLICY,
            "last_hidden_state_mean_over_non_special_nonpadding_tokens",
        )
        self.assertEqual(
            MODULE.CHECKPOINT_SHA256,
            "965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3",
        )
        self.assertEqual(
            MODULE.PRIOR_EMBEDDINGS_SHA256,
            "361ac7ec520bd6dcc2a9628af0c45f8b8dd7c12b4cbb895998e4b66c851150aa",
        )


if __name__ == "__main__":
    unittest.main()
