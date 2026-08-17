from __future__ import annotations

import unittest

import numpy as np
import torch

from masld_cl.sampling import (
    SamplingError,
    balanced_exact_fraction_indices,
    capped_group_indices,
    concatenate_query_and_replay,
    donor_train_validation_split,
    exact_quantile_replay,
    masked_augmentations,
    stratified_exact_bi_replay,
    stratified_fraction_indices,
)


class TestSampling(unittest.TestCase):
    def test_capped_and_stratified_are_reproducible(self):
        groups = [("D1", "hep")] * 10 + [("D2", "hep")] * 3
        a = capped_group_indices(groups, 4, seed=17)
        b = capped_group_indices(groups, 4, seed=17)
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(len(a), 7)
        strata = [groups[i] for i in a]
        replay = stratified_fraction_indices(a, strata, 0.5, seed=41)
        self.assertEqual(len(np.unique(replay)), len(replay))

    def test_replay_occurs_once(self):
        result = concatenate_query_and_replay([0, 1, 2], [3, 4])
        self.assertEqual(result.tolist(), [0, 1, 2, 3, 4])
        with self.assertRaises(SamplingError):
            concatenate_query_and_replay([0, 1], [1, 2])

    def test_balanced_fraction_is_exact_max_min_and_deterministic(self):
        candidates = np.arange(22)
        strata = [("small",)] * 2 + [("left",)] * 10 + [("right",)] * 10
        first = balanced_exact_fraction_indices(candidates, strata, 0.5, seed=17)
        second = balanced_exact_fraction_indices(candidates, strata, 0.5, seed=17)
        self.assertTrue(np.array_equal(first, second))
        self.assertEqual(len(first), 11)
        counts = {
            name: sum(strata[index][0] == name for index in first)
            for name in ("small", "left", "right")
        }
        self.assertEqual(counts["small"], 2)
        self.assertEqual(counts["left"] + counts["right"], 9)
        self.assertLessEqual(abs(counts["left"] - counts["right"]), 1)

    def test_balanced_fraction_rejects_duplicate_sources(self):
        with self.assertRaises(SamplingError):
            balanced_exact_fraction_indices(
                [0, 0], [("A",), ("B",)], 0.5, seed=17
            )

    def test_donor_split_has_no_leakage(self):
        donors = ["A", "A", "B", "B", "C", "C", "D", "D"]
        train, validation = donor_train_validation_split(donors, 0.25, 17)
        self.assertFalse(set(np.asarray(donors)[train]) & set(np.asarray(donors)[validation]))

    def test_bi_masks_exactly_half_and_selects_exact_size(self):
        x = torch.ones(3, 10)
        generator = torch.Generator().manual_seed(17)
        augmented = masked_augmentations(x, n_augmentations=200, mask_fraction=0.5, generator=generator)
        self.assertEqual(augmented.shape, (200, 3, 10))
        self.assertTrue(torch.equal((augmented == 0).sum(-1), torch.full((200, 3), 5)))
        scores = np.linspace(0, 1, 17)
        for mode in ("top", "bottom", "step"):
            selected = exact_quantile_replay(scores, 7, mode)
            self.assertEqual(len(selected), 7)
            self.assertEqual(len(np.unique(selected)), 7)

    def test_stratified_bi_returns_exact_balanced_buffer(self):
        candidates = np.arange(17)
        strata = [("A",)] * 10 + [("B",)] * 7
        result = stratified_exact_bi_replay(candidates, np.linspace(0, 1, 17), strata, 8, "top")
        self.assertEqual(len(result), 8)
        self.assertEqual(len(np.unique(result)), 8)
        self.assertEqual(np.sum(result < 10), 5)

    def test_stratified_bi_can_match_random_replay_stratum_counts(self):
        candidates = np.arange(17)
        strata = [("A",)] * 10 + [("B",)] * 7
        result = stratified_exact_bi_replay(
            candidates, np.linspace(0, 1, 17), strata, 8, "bottom",
            target_counts={("A",): 4, ("B",): 4},
        )
        self.assertEqual(len(result), 8)
        self.assertEqual(np.sum(result < 10), 4)


if __name__ == "__main__":
    unittest.main()
