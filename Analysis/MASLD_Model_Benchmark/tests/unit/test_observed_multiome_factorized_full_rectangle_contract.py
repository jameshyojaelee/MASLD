from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.audit_gse296875_observed_multiome_factorized_full_rectangle_contract import FullRectangleContractError, validate


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_factorized_full_rectangle_contract_20260825.json"


class FullRectangleContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_rectangle_passes_without_outcome_values(self) -> None:
        receipt = validate(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["full_rectangle_execution_authorized"])
        self.assertFalse(receipt["evaluator_artifact_bytes_read"])

    def test_missing_seed_is_rejected_before_artifact_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["rectangle"]["seeds"].pop()
        with self.assertRaises(FullRectangleContractError):
            validate(ROOT, mutated)

    def test_partial_ranking_is_rejected_before_artifact_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["firewall"]["partial_ranking_authorized"] = True
        with self.assertRaises(FullRectangleContractError):
            validate(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
