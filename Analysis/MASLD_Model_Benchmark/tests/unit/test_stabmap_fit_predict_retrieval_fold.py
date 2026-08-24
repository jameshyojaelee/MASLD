from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "stabmap_fit_predict_retrieval_fold.R"


class StabMapFoldSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = SCRIPT.read_text(encoding="utf-8")

    def test_only_training_bridge_fits_atac_idf(self) -> None:
        self.assertIn(
            "fit_atac_tfidf(\n    training_atac[, bridge_indices, drop = FALSE], query_atac_counts",
            self.source,
        )
        self.assertIn("detection <- Matrix::rowSums(bridge_binary)", self.source)

    def test_hidden_identity_map_is_not_an_argument(self) -> None:
        self.assertNotIn('file.path(input, "pair_map', self.source)
        self.assertNotIn("pair_map.tsv", self.source)
        self.assertIn('"hidden_pair_map_read",', self.source)
        self.assertIn('"false", "false",', self.source)

    def test_held_modalities_are_separate_terminal_queries(self) -> None:
        self.assertIn("query_rna = query_rna", self.source)
        self.assertIn("query_atac = query_atac", self.source)
        self.assertIn('reference_list = c("reference_rna")', self.source)

    def test_reference_only_rebase_and_query_order_check_are_explicit(self) -> None:
        self.assertIn("center <- colMeans(embedding[reference_ids, , drop = FALSE])", self.source)
        self.assertIn('fail("StabMap query-order invariance differs")', self.source)

    def test_feature_ranking_uses_training_matrices_only(self) -> None:
        self.assertIn("select_training_features(reference_rna, 1000L)", self.source)
        self.assertIn("select_training_features(bridge_atac, 1000L)", self.source)
        self.assertNotIn("select_training_features(query_", self.source)


if __name__ == "__main__":
    unittest.main()
