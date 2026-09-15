from __future__ import annotations

import unittest

import numpy as np

from scripts.audit_lsgkm_gse281364_development_evaluation import fisher_macro, spearman


class EvaluationAuditTest(unittest.TestCase):
    def test_spearman_handles_ties(self) -> None:
        left = np.asarray([1.0, 1.0, 2.0, 3.0])
        right = np.asarray([4.0, 4.0, 3.0, 1.0])
        self.assertAlmostEqual(spearman(left, right), -1.0)

    def test_fisher_macro_matches_single_value(self) -> None:
        self.assertAlmostEqual(fisher_macro([0.25]), 0.25)


if __name__ == "__main__":
    unittest.main()
