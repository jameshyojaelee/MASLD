from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.audit_observed_multiome_target_axis import TargetAxisContractError, validate


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/observed_multiome_target_axis_contract_20260825.json"


class ObservedMultiomeTargetAxisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_contract_passes_and_keeps_production_closed(self) -> None:
        receipt = validate(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["production_target_axis_frozen"])
        self.assertFalse(receipt["synthetic_fixture_production_eligible"])
        self.assertFalse(receipt["production_implementation_ready"])

    def test_held_block_target_parameter_fails(self) -> None:
        config = deepcopy(self.config)
        config["crossed_split"]["target_specific_parameter_learned_from_held_block"] = True
        with self.assertRaises(TargetAxisContractError):
            validate(ROOT, config)

    def test_target_intercept_from_held_block_fails(self) -> None:
        config = deepcopy(self.config)
        config["prediction_contract"]["target_specific_intercept_from_held_block_allowed"] = True
        with self.assertRaises(TargetAxisContractError):
            validate(ROOT, config)


if __name__ == "__main__":
    unittest.main()
