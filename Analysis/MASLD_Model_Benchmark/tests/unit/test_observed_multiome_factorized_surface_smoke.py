from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.run_gse296875_observed_multiome_factorized_surface_smoke import SurfaceSmokeError, validate_config


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_factorized_surface_smoke_20260825.json"


class ObservedMultiomeFactorizedSurfaceSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_surface_and_all_parent_artifacts_pass(self) -> None:
        parents, children = validate_config(ROOT, deepcopy(self.config))
        self.assertEqual(len(parents), 9)
        self.assertEqual(len(children), 5)

    def test_different_surface_is_rejected_before_artifact_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["surface"]["held_donor_fold"] = 1
        with self.assertRaises(SurfaceSmokeError):
            validate_config(ROOT, mutated)

    def test_ranking_authorization_is_rejected_before_artifact_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["evaluation"]["ranking_authorized"] = True
        with self.assertRaises(SurfaceSmokeError):
            validate_config(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
