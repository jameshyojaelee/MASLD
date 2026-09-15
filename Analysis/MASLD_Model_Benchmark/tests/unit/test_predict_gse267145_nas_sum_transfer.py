from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.predict_gse267145_nas_sum_transfer import (
    NasPredictionError,
    _min_pairwise_correlation,
    load_outcome_blind_roster,
)


def write_roster(path: Path, extra=None) -> None:
    fields = ["participant_index", "participant_id", "rna_source_sample_accession"]
    rows = [{"participant_index": str(i), "participant_id": p,
             "rna_source_sample_accession": f"GSM{i}"} for i, p in enumerate(["B", "A"])]
    if extra:
        fields += list(extra)
        for r in rows:
            r.update(extra)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


class RosterFirewallTests(unittest.TestCase):
    def test_clean_roster_sorts_by_participant(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "r.tsv"
            write_roster(p)
            self.assertEqual([r["participant_id"] for r in load_outcome_blind_roster(p)],
                             ["A", "B"])

    def test_each_outcome_column_raises(self) -> None:
        for col in ("nash_crn_component_sum", "steatosis", "ballooning",
                    "lobular_inflammation", "fibrosis", "stage3", "nas_score"):
            with self.subTest(col=col):
                with tempfile.TemporaryDirectory() as d:
                    p = Path(d) / "r.tsv"
                    write_roster(p, extra={col: "4"})
                    with self.assertRaises(NasPredictionError):
                        load_outcome_blind_roster(p)


class SeedDeterminismTests(unittest.TestCase):
    """The gate turns on this statistic, so it must be measured correctly."""

    def test_identical_seed_vectors_give_correlation_one(self) -> None:
        rng = np.random.default_rng(3)
        one = rng.normal(size=99)
        self.assertAlmostEqual(_min_pairwise_correlation(np.tile(one, (5, 1))), 1.0, places=12)

    def test_a_near_deterministic_fit_clears_the_0999_threshold(self) -> None:
        rng = np.random.default_rng(5)
        base = rng.normal(size=99)
        stack = np.stack([base + rng.normal(scale=1e-6, size=99) for _ in range(5)])
        self.assertGreater(_min_pairwise_correlation(stack), 0.999)

    def test_genuinely_varying_seeds_fall_below_the_threshold(self) -> None:
        rng = np.random.default_rng(7)
        self.assertLess(_min_pairwise_correlation(rng.normal(size=(5, 99))), 0.999)

    def test_the_minimum_not_the_mean_is_reported(self) -> None:
        rng = np.random.default_rng(11)
        base = rng.normal(size=99)
        stack = np.stack([base, base, base, base, rng.normal(size=99)])
        # Four identical vectors must not hide the one that differs.
        self.assertLess(_min_pairwise_correlation(stack), 0.999)


if __name__ == "__main__":
    unittest.main()
