"""A correction must be aimed at a sentence that is actually sealed.

The failure mode this pins is a correction that drifts off its target: the
wrong reason gets edited somewhere, the correction still runs, and the record
now carries a fix for a defect no artifact has. So the script refuses to build
unless the wrong sentence is still present in the sealed prespecification.

The df arithmetic is pinned directly, because it is the whole claim: the same
correlations must give a larger BH count at an inflated df than at the honest
one, and that must hold as a property of the function rather than as a number
observed once.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "record_bootstrap_df_correction.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


correction = _load(MODULE, "bootstrap_df_correction_tested")


class DfArithmeticTests(unittest.TestCase):
    def test_an_inflated_df_never_reports_fewer_discoveries(self) -> None:
        generator = np.random.default_rng(1)
        for _ in range(30):
            correlations = np.clip(generator.normal(0.0, 0.22, 3000), -0.99, 0.99)
            self.assertGreaterEqual(
                correction.bh_count(correlations, 96),
                correction.bh_count(correlations, 59),
            )

    def test_the_gap_is_material_at_the_realised_duplication_rate(self) -> None:
        """62.5 distinct of 99 is the measured rate; the gap must be real there."""

        generator = np.random.default_rng(2)
        correlations = np.concatenate([
            generator.normal(0.0, 0.12, 5000),
            generator.normal(0.42, 0.05, 200),
        ])
        inflated = correction.bh_count(correlations, 96)
        honest = correction.bh_count(correlations, 62 - 3)
        self.assertGreater(inflated, honest)

    def test_identical_df_gives_identical_counts(self) -> None:
        generator = np.random.default_rng(3)
        correlations = np.clip(generator.normal(0.0, 0.3, 1000), -0.99, 0.99)
        self.assertEqual(
            correction.bh_count(correlations, 96), correction.bh_count(correlations, 96)
        )

    def test_nothing_significant_counts_zero_at_either_df(self) -> None:
        correlations = np.full(500, 0.01)
        self.assertEqual(correction.bh_count(correlations, 96), 0)
        self.assertEqual(correction.bh_count(correlations, 59), 0)


class TargetingTests(unittest.TestCase):
    def test_the_wrong_reason_substring_is_the_one_being_corrected(self) -> None:
        self.assertIn("concordance blocks", correction.WRONG_REASON_SUBSTRING)

    def test_the_corrected_reason_names_the_df_not_duplication_alone(self) -> None:
        self.assertIn("df = distinct - 3", correction.CORRECTED_REASON)
        self.assertIn("degrees of freedom were the bug", correction.CORRECTED_REASON)

    def test_the_standing_rule_states_its_own_scope(self) -> None:
        self.assertIn("df-dependent threshold", correction.STANDING_RULE)
        self.assertIn("does not reach them", correction.STANDING_RULE)

    def test_the_diagnostic_shape_is_recorded(self) -> None:
        self.assertIn("below its own bootstrap", correction.DIAGNOSTIC_SHAPE)


class HelperTests(unittest.TestCase):
    def test_unit_ranks_are_centred_and_unit_norm(self) -> None:
        vector = correction.unit_ranks(np.asarray([0, 0, 1, 2, 2, 3], dtype=float))
        self.assertAlmostEqual(float(vector.sum()), 0.0, places=12)
        self.assertAlmostEqual(float(vector @ vector), 1.0, places=12)

    def test_standardize_flags_a_constant_column(self) -> None:
        block = np.column_stack([np.asarray([1.0, 2.0, 3.0]), np.ones(3)])
        scaled, keep = correction.standardize(block)
        self.assertTrue(keep[0])
        self.assertFalse(keep[1])

    def test_vector_spearman_is_one_on_a_monotone_copy(self) -> None:
        left = np.asarray([0.1, 0.4, 0.2, 0.9, 0.5])
        self.assertAlmostEqual(correction.vector_spearman(left, 3.0 * left), 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
