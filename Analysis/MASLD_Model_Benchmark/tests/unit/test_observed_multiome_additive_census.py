from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.freeze_observed_multiome_additive_census import (
    EXPECTED_NEW_BASELINES,
    ObservedMultiomeCensusError,
    validate_revision,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/observed_multiome_additive_census_revision_20260825.json"


class ObservedMultiomeAdditiveCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_revision_passes_without_authorizing_execution(self) -> None:
        receipt = validate_revision(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["census_requirement_satisfied"])
        self.assertFalse(receipt["candidate_fixture_requirement_satisfied"])
        self.assertFalse(receipt["execution_authorized"])
        self.assertEqual(receipt["new_baseline_count"], len(EXPECTED_NEW_BASELINES))

    def test_duplicate_new_baseline_fails(self) -> None:
        config = deepcopy(self.config)
        config["new_mandatory_baselines"].append(
            deepcopy(config["new_mandatory_baselines"][0])
        )
        with self.assertRaises(ObservedMultiomeCensusError):
            validate_revision(ROOT, config)

    def test_missing_as_zero_fails(self) -> None:
        config = deepcopy(self.config)
        config["execution_contract"]["missing_as_zero"] = True
        with self.assertRaises(ObservedMultiomeCensusError):
            validate_revision(ROOT, config)


if __name__ == "__main__":
    unittest.main()
