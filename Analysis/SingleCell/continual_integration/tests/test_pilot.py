from __future__ import annotations

import unittest

from masld_cl.pilot import _pilot_passes


class TestPilot(unittest.TestCase):
    def setUp(self):
        self.config = {"gates": {
            "reference_macro_f1_margin": 0.02,
            "reference_neighborhood_jaccard_loss": 0.05,
            "maximum_stratum_worsening": 0.10,
        }}

    def test_pilot_allows_tuning_after_noncatastrophic_paper_setting(self):
        result = _pilot_passes(self.config, {
            "reference_macro_f1_change_ci_low": -0.01,
            "reference_neighborhood_jaccard_loss": 0.04,
            "alignment_improvement_vs_harmony": -0.05,
            "alignment_improvement_vs_architecture_surgery": 0.02,
        })
        self.assertTrue(all(result.values()))

    def test_pilot_stops_on_reference_failure(self):
        result = _pilot_passes(self.config, {
            "reference_macro_f1_change_ci_low": -0.03,
            "reference_neighborhood_jaccard_loss": 0.04,
            "alignment_improvement_vs_harmony": 0.2,
            "alignment_improvement_vs_architecture_surgery": 0.2,
        })
        self.assertFalse(result["reference_retention"])


if __name__ == "__main__":
    unittest.main()
