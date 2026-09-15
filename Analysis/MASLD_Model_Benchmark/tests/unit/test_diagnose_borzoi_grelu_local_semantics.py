from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

import torch

from scripts.diagnose_borzoi_grelu_local_semantics import (
    BorzoiSemanticsDiagnosticError,
    fixed_training_length_central_mask,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class BorzoiSemanticsDiagnosticTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(
            (ROOT / "config/borzoi_grelu_local_semantics_diagnostic.json").read_text()
        )

    def test_contract_is_explanatory_only(self) -> None:
        validate_contract(self.contract)
        self.assertFalse(
            self.contract["execution_contract"]["local_port_execution_gate_may_open"]
        )
        self.assertFalse(
            self.contract["claim_boundary"]["native_numeric_parity_established"]
        )

    def test_contract_rejects_post_hoc_promotion(self) -> None:
        changed = deepcopy(self.contract)
        changed["execution_contract"]["local_port_execution_gate_may_open"] = True
        with self.assertRaisesRegex(
            BorzoiSemanticsDiagnosticError, "opened or drifted"
        ):
            validate_contract(changed)

    def test_fixed_training_length_mask_is_length_stable(self) -> None:
        first = torch.zeros((1, 8, 16), dtype=torch.float32)
        second = torch.zeros((2, 8, 16), dtype=torch.float32)
        first_mask = fixed_training_length_central_mask(first, 32, 4096)
        second_mask = fixed_training_length_central_mask(second, 32, 4096)
        self.assertTrue(torch.equal(first_mask, second_mask))
        self.assertEqual(tuple(first_mask.shape), (15, 32))


if __name__ == "__main__":
    unittest.main()
