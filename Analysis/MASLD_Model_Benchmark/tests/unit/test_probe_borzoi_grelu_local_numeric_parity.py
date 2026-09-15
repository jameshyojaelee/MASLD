from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

import numpy as np
import torch

from scripts.probe_borzoi_grelu_local_numeric_parity import (
    BorzoiNumericParityError,
    compare_tensors,
    exponential_linspace_int,
    relative_shift,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class BorzoiNumericParityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(
            (
                ROOT / "config/borzoi_grelu_local_numeric_parity_fixture.json"
            ).read_text()
        )

    def test_contract_is_bounded_and_fail_closed(self) -> None:
        validate_contract(self.contract)
        self.assertEqual(self.contract["fixture"]["input_length_bp"], 1024)
        self.assertFalse(
            self.contract["historical_implementation"][
                "source_revision_is_checkpoint_bound"
            ]
        )
        self.assertFalse(
            self.contract["execution_contract"]["checkpoint_export_allowed"]
        )

    def test_contract_rejects_claim_or_tolerance_drift(self) -> None:
        changed = deepcopy(self.contract)
        changed["claim_boundary"]["open_champion_eligible"] = True
        with self.assertRaisesRegex(BorzoiNumericParityError, "opened or drifted"):
            validate_contract(changed)
        changed = deepcopy(self.contract)
        changed["comparison"]["absolute_tolerance"] = 1.0
        with self.assertRaisesRegex(BorzoiNumericParityError, "opened or drifted"):
            validate_contract(changed)

    def test_enformer_import_primitives(self) -> None:
        self.assertEqual(
            exponential_linspace_int(608, 1536, 6, 32),
            [608, 736, 896, 1056, 1280, 1536],
        )
        value = torch.arange(1 * 2 * 3 * 5, dtype=torch.float32).reshape(1, 2, 3, 5)
        shifted = relative_shift(value)
        self.assertEqual(tuple(shifted.shape), (1, 2, 3, 3))
        self.assertTrue(torch.isfinite(shifted).all())

    def test_comparison_records_positive_and_negative_fixtures(self) -> None:
        base = torch.tensor([[1.0, 2.0, 3.0]], dtype=torch.float32)
        identical = compare_tensors(base, base.clone(), atol=1e-6, rtol=1e-5)
        self.assertTrue(identical["all_elements_close"])
        self.assertTrue(identical["bit_identical"])
        changed = compare_tensors(base, base + 0.1, atol=1e-6, rtol=1e-5)
        self.assertFalse(changed["all_elements_close"])
        self.assertLess(changed["close_fraction"], 1.0)
        self.assertTrue(np.isfinite(changed["pearson_correlation"]))


if __name__ == "__main__":
    unittest.main()
