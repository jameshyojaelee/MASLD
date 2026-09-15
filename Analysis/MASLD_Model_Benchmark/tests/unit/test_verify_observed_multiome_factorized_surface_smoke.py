from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.verify_gse296875_observed_multiome_factorized_surface_smoke import SurfaceSmokeVerificationError, validate_config


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_factorized_surface_smoke_verification_20260825.json"


class SurfaceSmokeVerificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_artifacts_pass_contract_validation(self) -> None:
        self.assertEqual(len(validate_config(ROOT, deepcopy(self.config))), 2)

    def test_open_ranking_is_rejected_before_artifact_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["firewall"]["partial_ranking_authorized"] = True
        with self.assertRaises(SurfaceSmokeVerificationError):
            validate_config(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
