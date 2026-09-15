from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import tomllib
import unittest

from scripts.freeze_observed_multiome_development_contract import (
    ObservedMultiomeContractError,
    validate_capabilities,
    validate_contract,
    validate_gate,
)


ROOT = Path(__file__).resolve().parents[2]


def _toml(relative: str) -> dict:
    with (ROOT / relative).open("rb") as handle:
        return tomllib.load(handle)


class ObservedMultiomeDevelopmentContractTests(unittest.TestCase):
    def test_standalone_contract_is_fail_closed(self) -> None:
        receipt = validate_contract(ROOT)
        self.assertFalse(receipt["execution_authorized"])
        self.assertFalse(receipt["external_champion_claim_allowed"])
        self.assertEqual(receipt["primary_eligible_models"], [])

    def test_gate_rejects_model_promotion(self) -> None:
        gate = deepcopy(
            _toml("config/evaluation/observed_multiome_development_gate.toml")
        )
        gate["model_promotion_allowed"] = True
        with self.assertRaisesRegex(ObservedMultiomeContractError, "gate opened"):
            validate_gate(gate)

    def test_capabilities_reject_false_pairing(self) -> None:
        capabilities = deepcopy(
            _toml("config/evaluation/observed_multiome_task_capabilities.toml")
        )
        capabilities["false_cell_pairing_forbidden"] = False
        task = _toml("config/evaluation/observed_multiome_task.toml")
        model_ids = set()
        for path in sorted((ROOT / "config/models").glob("*.toml")):
            for model in _toml(str(path.relative_to(ROOT))).get("models", []):
                model_ids.add(model["model_id"])
        with self.assertRaisesRegex(
            ObservedMultiomeContractError, "safety gate opened"
        ):
            validate_capabilities(
                capabilities,
                task_baselines=tuple(task["baseline_model_ids"]),
                model_ids=model_ids,
            )


if __name__ == "__main__":
    unittest.main()
