"""The overlay must stay an overlay, and must reuse the frozen code path.

Two failure modes are pinned. The first is the overlay quietly becoming a gate,
which is how a derived-after number ends up cited as pre-registered. The second
is the overlay reimplementing the partial-correlation machinery instead of
importing it, which would make its numbers incomparable to the four families in
the frozen result while still looking comparable.

The reliability prediction is pinned as a literal because it must have been
recorded before the measurement existed; a test that reads it back from the
output would prove nothing.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "record_axis_count_derived_after_addendum.py"
ANALYSIS = ROOT / "scripts" / "evaluate_gse267145_axis_count.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


addendum = _load(MODULE, "axis_count_addendum_tested")


class OverlayStaysAnOverlayTests(unittest.TestCase):
    def test_the_prediction_is_a_literal_recorded_before_the_measurement(self) -> None:
        self.assertIn("higher reliability than either component", addendum.RELIABILITY_PREDICTION)
        self.assertIn("0.479", addendum.RELIABILITY_PREDICTION)
        self.assertIn("0.518", addendum.RELIABILITY_PREDICTION)

    def test_the_analysis_module_is_imported_not_reimplemented(self) -> None:
        analysis = addendum._load_analysis(ANALYSIS)
        for name in (
            "partial_family",
            "unit_ranks",
            "residualize",
            "unit_columns",
            "average_ranks",
            "centred_ranks",
        ):
            self.assertTrue(hasattr(analysis, name), f"analysis lacks {name}")

    def test_the_addendum_defines_no_partial_correlation_of_its_own(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        for forbidden in ("def partial_family", "def residualize", "def unit_ranks"):
            self.assertNotIn(forbidden, source)


class SplitHalfTests(unittest.TestCase):
    def test_a_perfectly_reproducible_axis_reads_high(self) -> None:
        analysis = addendum._load_analysis(ANALYSIS)
        generator = np.random.default_rng(5)
        n, g = 60, 400
        outcome = generator.integers(0, 5, size=n).astype(float)
        genes = generator.normal(size=(n, g))
        genes[:, :200] += 3.0 * (outcome - outcome.mean())[:, None]
        summary = addendum.split_half_reliability(
            analysis, genes, {"strong": outcome}, n_splits=25, seed=3
        )
        self.assertGreater(summary["strong"]["spearman_brown_full_length_reliability"], 0.8)

    def test_a_pure_noise_axis_reads_near_zero(self) -> None:
        analysis = addendum._load_analysis(ANALYSIS)
        generator = np.random.default_rng(6)
        n, g = 60, 400
        summary = addendum.split_half_reliability(
            analysis,
            generator.normal(size=(n, g)),
            {"noise": generator.integers(0, 5, size=n).astype(float)},
            n_splits=25,
            seed=4,
        )
        self.assertLess(abs(summary["noise"]["median_half_length_reliability"]), 0.25)

    def test_the_surviving_gene_count_is_reported_not_assumed(self) -> None:
        analysis = addendum._load_analysis(ANALYSIS)
        generator = np.random.default_rng(7)
        genes = generator.normal(size=(40, 100))
        genes[:, :10] = 1.0  # constant everywhere, must be dropped
        summary = addendum.split_half_reliability(
            analysis,
            genes,
            {"axis": generator.integers(0, 4, size=40).astype(float)},
            n_splits=10,
            seed=8,
        )
        self.assertLessEqual(
            summary["axis"]["median_genes_non_constant_in_both_halves"], 90.0
        )


if __name__ == "__main__":
    unittest.main()
