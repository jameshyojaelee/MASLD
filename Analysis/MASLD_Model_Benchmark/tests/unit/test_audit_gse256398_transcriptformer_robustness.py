from __future__ import annotations

import unittest

import numpy as np

from scripts.audit_gse256398_transcriptformer_robustness import (
    GSE256398RobustnessError,
    cosine_distance,
    summarize,
)


class GSE256398TranscriptFormerRobustnessTests(unittest.TestCase):
    def test_cosine_distance_identity_and_orthogonality(self) -> None:
        self.assertAlmostEqual(cosine_distance(np.array([1.0, 0.0]), np.array([1.0, 0.0])), 0.0)
        self.assertAlmostEqual(cosine_distance(np.array([1.0, 0.0]), np.array([0.0, 1.0])), 1.0)

    def test_zero_norm_is_rejected(self) -> None:
        with self.assertRaises(GSE256398RobustnessError):
            cosine_distance(np.zeros(2), np.ones(2))

    def test_summary_is_descriptive(self) -> None:
        self.assertEqual(
            summarize([1.0, 2.0, 4.0]),
            {"minimum": 1.0, "median": 2.0, "maximum": 4.0},
        )


if __name__ == "__main__":
    unittest.main()
