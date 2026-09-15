from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from scripts.freeze_lsgkm_gse281364_full_campaign_authorization import (
    FullCampaignAuthorizationError,
    expected_run_fits,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/lsgkm_gse281364_dinucleotide_full_campaign.json"


class FullCampaignAuthorizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_exact_remaining_rectangle_excludes_only_reused_probe(self) -> None:
        validate_config(self.config)
        fits = expected_run_fits(self.config)
        self.assertEqual(len(fits), 24)
        self.assertEqual(len(set(fits)), 24)
        self.assertNotIn(("donor0_genomic0", 1103), fits)

    def test_probe_fit_cannot_be_retrained(self) -> None:
        changed = copy.deepcopy(self.config)
        changed["fit_grid"]["splits"]["donor0_genomic0"]["run_seeds"].append(1103)
        with self.assertRaises(FullCampaignAuthorizationError):
            validate_config(changed)

    def test_outcomes_metrics_controls_aggregation_and_seals_remain_blocked(self) -> None:
        validate_config(self.config)
        firewall = self.config["action_firewall"]
        for field in (
            "outcome_access_authorized",
            "reporter_count_access_authorized",
            "control_prediction_value_access_authorized",
            "sealed_asset_access_authorized",
            "benchmark_metric_calculation_authorized",
            "aggregation_or_evaluation_authorized",
        ):
            self.assertFalse(firewall[field])

    def test_bundles_match_native_single_thread_fits(self) -> None:
        validate_config(self.config)
        execution = self.config["execution"]
        self.assertEqual(
            execution["requested_cpu_count_by_split"],
            {"donor0_genomic0": 4, "donor1_genomic1": 5, "donor2_genomic2": 5, "donor3_genomic3": 5, "donor4_genomic4": 5},
        )


if __name__ == "__main__":
    unittest.main()
