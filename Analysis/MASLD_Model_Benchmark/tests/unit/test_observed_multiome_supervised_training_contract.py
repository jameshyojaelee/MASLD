from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.audit_gse296875_observed_multiome_supervised_training_contract import (
    SupervisedTrainingContractError,
    validate,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_supervised_training_contract_20260825.json"


class ObservedMultiomeSupervisedTrainingContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_contract_closes_supervised_training_gap(self) -> None:
        receipt = validate(ROOT, deepcopy(self.config))
        self.assertEqual(receipt["training_label_artifacts_planned"], 25)
        self.assertFalse(receipt["model_fit_may_bind_evaluator_outcomes"])
        self.assertFalse(receipt["supervised_fit_authorized"])

    def test_held_donor_training_values_are_rejected(self) -> None:
        mutated = deepcopy(self.config)
        mutated["training_label_artifacts"]["held_donor_rows"] = "allowed"
        with self.assertRaises(SupervisedTrainingContractError):
            validate(ROOT, mutated)

    def test_evaluator_binding_to_fit_is_rejected(self) -> None:
        mutated = deepcopy(self.config)
        mutated["run_binding_v2"]["fit_may_not_bind"].remove("any_evaluator_outcome_artifact_or_verifier")
        with self.assertRaises(SupervisedTrainingContractError):
            validate(ROOT, mutated)

    def test_target_selection_using_counts_is_rejected(self) -> None:
        mutated = deepcopy(self.config)
        mutated["training_target_plan"]["selection_uses_count_values"] = True
        with self.assertRaises(SupervisedTrainingContractError):
            validate(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
