from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.probe_borzoi_grelu_converted_port_weights_only import (
    WeightsOnlyProbeError,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class BorzoiGreluWeightsOnlyProbeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(
            (
                ROOT
                / "config/borzoi_grelu_converted_port_weights_only_probe.json"
            ).read_text()
        )

    def test_contract_is_weights_only_and_fail_closed(self) -> None:
        validate_contract(self.contract)
        self.assertTrue(self.contract["load_contract"]["torch_load_weights_only"])
        self.assertFalse(self.contract["load_contract"]["torch_load_unrestricted"])
        self.assertFalse(self.contract["load_contract"]["model_forward_allowed"])

    def test_contract_rejects_model_forward(self) -> None:
        changed = deepcopy(self.contract)
        changed["load_contract"]["model_forward_allowed"] = True
        with self.assertRaisesRegex(WeightsOnlyProbeError, "contract opened"):
            validate_contract(changed)

    def test_contract_rejects_unrestricted_load(self) -> None:
        changed = deepcopy(self.contract)
        changed["load_contract"]["torch_load_unrestricted"] = True
        with self.assertRaisesRegex(WeightsOnlyProbeError, "contract opened"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
