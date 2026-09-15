from __future__ import annotations

import unittest

import numpy as np

from scripts.verify_gse296875_observed_multiome_complete_rectangle import independent_profile_skill, independent_two_way_interval


class VerifyObservedMultiomeCompleteRectangleTests(unittest.TestCase):
    def test_exact_profile_has_unit_skill(self) -> None:
        observed = np.asarray([[8.0, 2.0, 1.0], [1.0, 2.0, 8.0]])
        skill = independent_profile_skill(observed, observed)
        np.testing.assert_allclose(skill, np.ones(2), atol=1e-12, rtol=0)

    def test_two_way_interval_uses_39_donors_and_five_blocks(self) -> None:
        donors = [f"donor-{donor:02d}" for block in range(5) for donor in range(39) for _ in range(5)]
        blocks = [block for block in range(5) for _ in range(39) for _ in range(5)]
        values = [0.1 + 0.001 * lineage for _block in range(5) for _donor in range(39) for lineage in range(5)]
        result = independent_two_way_interval(values, donors, blocks, n_resamples=100, seed=7)
        self.assertEqual(result["n_donors"], 39)
        self.assertEqual(result["n_genomic_blocks"], 5)
        self.assertEqual(result["n_cells"], 195)
        self.assertEqual(result["n_observations"], 975)
        self.assertAlmostEqual(result["estimate"], 0.102)


if __name__ == "__main__":
    unittest.main()
