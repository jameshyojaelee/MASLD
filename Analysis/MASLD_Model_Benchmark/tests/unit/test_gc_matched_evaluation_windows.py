from __future__ import annotations

import unittest

from scripts import build_gc_matched_evaluation_windows as builder


class GCMatchedEvaluationWindowTests(unittest.TestCase):
    def test_merge_and_half_open_overlap(self) -> None:
        merged = builder.merge_intervals(((10, 20), (18, 30), (40, 50)))
        self.assertEqual(merged, ((10, 30), (40, 50)))
        starts = tuple(value[0] for value in merged)
        self.assertTrue(builder.overlaps(merged, starts, 0, 11))
        self.assertFalse(builder.overlaps(merged, starts, 30, 40))
        self.assertTrue(builder.overlaps(merged, starts, 49, 60))

    def test_gc_match_is_unique_and_deterministic(self) -> None:
        targets = ((0.21, "a"), (0.51, "b"), (0.79, "c"))
        candidates = ((0.20, 100), (0.50, 200), (0.80, 300), (0.90, 400))
        expected = {"a": (0.20, 100), "b": (0.50, 200), "c": (0.80, 300)}
        self.assertEqual(builder.match_gc_targets(targets, candidates), expected)
        self.assertEqual(
            builder.match_gc_targets(tuple(reversed(targets)), tuple(reversed(candidates))),
            expected,
        )

    def test_gc_fraction_rejects_ambiguous_sequence(self) -> None:
        sequence = "GA" * (builder.INPUT_LENGTH // 2)
        self.assertEqual(builder.gc_fraction(sequence), 0.5)
        with self.assertRaises(builder.BackgroundWindowError):
            builder.gc_fraction("N" + sequence[1:])


if __name__ == "__main__":
    unittest.main()
