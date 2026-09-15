from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.audit_gse296875_observed_multiome_materialization_contract import (
    ObservedMultiomeMaterializationContractError,
    validate,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_materialization_contract_20260825.json"


class GSE296875ObservedMultiomeMaterializationContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_contract_passes_without_count_access(self) -> None:
        receipt = validate(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["same_fold_held_chromosomes_excluded"])
        self.assertFalse(receipt["biological_count_values_read"])
        self.assertEqual(receipt["model_input_artifacts_planned"], 5)
        self.assertEqual(receipt["evaluator_outcome_artifacts_planned"], 5)

    def test_full_atac_model_load_is_rejected(self) -> None:
        config = deepcopy(self.config)
        config["selective_sparse_read"]["loading_full_atac_csr_data_array_in_model_materializer"] = True
        with self.assertRaises(ObservedMultiomeMaterializationContractError):
            validate(ROOT, config)

    def test_evaluator_binding_to_fit_is_rejected(self) -> None:
        config = deepcopy(self.config)
        config["run_binding"]["fit_and_predict_may_not_bind"].remove(
            "any_evaluator_outcome_artifact"
        )
        with self.assertRaises(ObservedMultiomeMaterializationContractError):
            validate(ROOT, config)

    def test_cell_level_replication_is_rejected(self) -> None:
        config = deepcopy(self.config)
        config["biological_unit"]["cells_or_nuclei_as_replicates"] = True
        with self.assertRaises(ObservedMultiomeMaterializationContractError):
            validate(ROOT, config)


if __name__ == "__main__":
    unittest.main()
