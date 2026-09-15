from __future__ import annotations

import unittest

from scripts.adjudicate_two_axis_external_arms import (
    attainable_ceiling,
    detectable_effect_block,
    minimum_detectable_rho,
)
from scripts.build_showcase_aspect_lineage_substrate import SubstrateError


class AttainableCeilingTests(unittest.TestCase):
    def test_an_untied_outcome_has_a_ceiling_of_one(self) -> None:
        self.assertAlmostEqual(attainable_ceiling([1, 1, 1, 1, 1]), 1.0, places=12)

    def test_ties_lower_the_ceiling_below_one(self) -> None:
        self.assertLess(attainable_ceiling([25, 28, 9, 14, 2]), 1.0)

    def test_a_dominant_level_lowers_the_ceiling_further(self) -> None:
        """GSE267145's 71-of-99 stage-0 mass is the case this must capture."""

        balanced = attainable_ceiling([25, 25, 25, 24])
        dominant = attainable_ceiling([71, 15, 9, 4])
        self.assertLess(dominant, balanced)
        self.assertLess(dominant, 0.85)

    def test_ceiling_is_symmetric_under_reversal(self) -> None:
        self.assertAlmostEqual(
            attainable_ceiling([6, 37, 37, 26]),
            attainable_ceiling([26, 37, 37, 6]),
            places=12,
        )

    def test_a_single_level_outcome_fails_closed(self) -> None:
        with self.assertRaises(SubstrateError):
            attainable_ceiling([40])

    def test_an_empty_level_fails_closed(self) -> None:
        """A zero-count level means the marginal was built wrong, not that a
        level exists with no members."""

        with self.assertRaises(SubstrateError):
            attainable_ceiling([10, 0, 10])


class MinimumDetectableRhoTests(unittest.TestCase):
    def test_more_participants_detect_smaller_effects(self) -> None:
        self.assertLess(minimum_detectable_rho(200, 0.05), minimum_detectable_rho(50, 0.05))

    def test_a_stricter_alpha_needs_a_larger_effect(self) -> None:
        self.assertGreater(
            minimum_detectable_rho(100, 0.05 / 12), minimum_detectable_rho(100, 0.05)
        )

    def test_known_value(self) -> None:
        # n=78, alpha=0.05, power=0.80
        self.assertAlmostEqual(minimum_detectable_rho(78, 0.05), 0.3127, places=3)

    def test_tiny_n_fails_closed(self) -> None:
        with self.assertRaises(SubstrateError):
            minimum_detectable_rho(4, 0.05)


class DetectableEffectBlockTests(unittest.TestCase):
    def test_headroom_uses_the_strictest_family_not_the_loosest(self) -> None:
        block = detectable_effect_block(
            label="t", group_sizes=[25, 28, 9, 14, 2], family_sizes=[4, 12]
        )
        mde = block["minimum_detectable_rho_at_80_percent_power"]
        self.assertEqual(
            block["headroom_at_the_strictest_family"],
            block["attainable_ceiling_rho"] - max(mde.values()),
        )
        self.assertTrue(block["arm_can_detect_an_effect_below_its_ceiling"])

    def test_an_arm_whose_ceiling_sits_under_its_mde_is_reported_as_such(self) -> None:
        """A tiny, heavily tied outcome cannot detect anything, and the record
        must say so rather than emitting a reassuring number."""

        block = detectable_effect_block(
            label="t", group_sizes=[18, 1, 1], family_sizes=[12]
        )
        self.assertFalse(block["arm_can_detect_an_effect_below_its_ceiling"])
        self.assertLess(block["headroom_at_the_strictest_family"], 0.0)

    def test_smallest_level_is_reported(self) -> None:
        block = detectable_effect_block(
            label="t", group_sizes=[25, 28, 9, 14, 2], family_sizes=[12]
        )
        self.assertEqual(block["smallest_level_n"], 2)
        self.assertEqual(block["n"], 78)
        self.assertEqual(block["levels_realised"], 5)


if __name__ == "__main__":
    unittest.main()
