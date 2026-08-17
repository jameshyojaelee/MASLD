import unittest

import numpy as np

from masld_cl.orthogonal_bridge_outcome import evaluate_matched_scope
from masld_cl.orthogonal_bridge_invariants import apply_offsets_in_order


class TestOrthogonalBridgeOutcome(unittest.TestCase):
    def test_matched_scope_is_invariant_to_scaled_rotation_and_translation(self):
        raw = np.random.default_rng(17).normal(size=(6, 2))
        controls = np.asarray([True, True, True, False, False, False])
        rotation = np.asarray([[0.0, -1.0], [1.0, 0.0]])
        candidate = 0.37 * raw @ rotation + np.asarray([8.0, -3.0])
        result = evaluate_matched_scope(candidate, raw, controls)
        self.assertIsNotNone(result)
        self.assertTrue(np.isclose(result["retention"], 1.0, rtol=1e-12, atol=1e-12))
        self.assertTrue(np.isclose(result["distance_spearman"], 1.0, rtol=1e-12, atol=1e-12))

    def test_matched_scope_is_untestable_below_three_donors_per_group(self):
        raw = np.arange(10, dtype=float).reshape(5, 2)
        self.assertIsNone(
            evaluate_matched_scope(raw, raw, np.asarray([True, True, False, False, False]))
        )

    def test_dataset_offsets_commute_across_acquisition_orders(self):
        parent = np.arange(12, dtype=np.float32).reshape(6, 2)
        datasets = np.asarray(["A", "A", "B", "B", "R", "R"])
        offsets = {"A": np.asarray([1.0, -1.0]), "B": np.asarray([-2.0, 2.0])}
        forward = apply_offsets_in_order(parent, datasets, offsets, ["A", "B"])
        reverse = apply_offsets_in_order(parent, datasets, offsets, ["B", "A"])
        np.testing.assert_array_equal(forward, reverse)
        np.testing.assert_array_equal(forward[datasets == "R"], parent[datasets == "R"])
