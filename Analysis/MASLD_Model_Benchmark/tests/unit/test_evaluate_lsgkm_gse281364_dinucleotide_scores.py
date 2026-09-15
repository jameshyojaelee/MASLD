from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np

from scripts.evaluate_lsgkm_gse281364_dinucleotide_scores import (
    CANDIDATES,
    fast_spearman,
    fisher_macro,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/lsgkm_gse281364_dinucleotide_development_evaluation.json"


class LSGKMEvaluationTests(unittest.TestCase):
    def test_locked_candidate_roster_and_reporting_limits(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        validate_config(config)
        self.assertEqual(len(CANDIDATES), 3)
        self.assertFalse(config["reporting_limits"]["external_evaluation"])
        self.assertFalse(config["reporting_limits"]["champion_claim"])

    def test_fast_spearman_matches_monotonic_direction(self) -> None:
        observed = np.asarray([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(fast_spearman(observed, observed), 1.0)
        self.assertAlmostEqual(fast_spearman(observed, observed[::-1]), -1.0)

    def test_fisher_macro_is_symmetric(self) -> None:
        self.assertAlmostEqual(fisher_macro([0.2, 0.4]), fisher_macro([0.4, 0.2]))


if __name__ == "__main__":
    unittest.main()
