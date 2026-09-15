from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from scripts.run_lsgkm_gse281364_scale_probe import (
    LSGKMScaleProbeError,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/lsgkm_gse281364_dinucleotide_scale_probe.json"


class LSGKMScaleProbeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_only_one_exact_fit_is_authorized(self) -> None:
        validate_config(self.config)
        self.assertEqual(self.config["active_fit"]["split_id"], "donor0_genomic0")
        self.assertEqual(self.config["active_fit"]["model_seed"], 1103)
        self.assertTrue(
            self.config["action_firewall"]["one_fit_production_training_authorized"]
        )
        self.assertFalse(self.config["action_firewall"]["remaining_24_fits_authorized"])

    def test_fit_count_or_seed_cannot_expand_silently(self) -> None:
        changed = copy.deepcopy(self.config)
        changed["active_fit"]["model_seed"] = 2909
        with self.assertRaises(LSGKMScaleProbeError):
            validate_config(changed)

    def test_outcomes_metrics_controls_and_seals_remain_blocked(self) -> None:
        validate_config(self.config)
        firewall = self.config["action_firewall"]
        self.assertFalse(firewall["outcome_access_authorized"])
        self.assertFalse(firewall["reporter_count_access_authorized"])
        self.assertFalse(firewall["control_prediction_value_access_authorized"])
        self.assertFalse(firewall["sealed_asset_access_authorized"])
        self.assertFalse(firewall["benchmark_metric_calculation_authorized"])

    def test_native_command_and_weight_are_frozen(self) -> None:
        validate_config(self.config)
        argv = self.config["runtime"]["gkmtrain_argv"]
        self.assertEqual(argv[argv.index("-w") + 1], "1")
        self.assertEqual(argv[argv.index("-T") + 1], "1")
        self.assertNotIn("-r", argv)
        self.assertNotIn("-x", argv)


if __name__ == "__main__":
    unittest.main()
