from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.freeze_corgi_regular_film_plus_head_contract import (
    CorgiFilmHeadContractError,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/corgi_regular_film_plus_head_smoke.toml"


def mutated_contract(old: str, new: str, temporary: str) -> Path:
    text = CONTRACT.read_text(encoding="utf-8")
    if text.count(old) != 1:
        raise AssertionError(f"mutation target must occur once: {old}")
    path = Path(temporary) / "contract.toml"
    path.write_text(text.replace(old, new), encoding="utf-8")
    return path


class CorgiRegularFilmPlusHeadContractTests(unittest.TestCase):
    def test_contract_is_outcome_blind_and_execution_blocked(self) -> None:
        receipt = validate_contract(ROOT, CONTRACT)
        self.assertEqual(
            receipt["status"], "pass_prospective_contract_execution_blocked"
        )
        self.assertEqual(receipt["crossfit_folds"], 5)
        self.assertEqual(receipt["trainable_parameter_numel"], 5869229)
        self.assertEqual(receipt["film_trainable_parameter_numel"], 5867308)
        self.assertEqual(receipt["replacement_head_parameter_numel"], 1921)
        self.assertIs(receipt["development_outcomes_read"], False)
        self.assertIs(receipt["training_execution_authorized"], False)
        self.assertIs(receipt["evaluation_execution_authorized"], False)

    def test_rejects_open_training_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = mutated_contract(
                "training_execution_authorized = false",
                "training_execution_authorized = true",
                temporary,
            )
            with self.assertRaises(CorgiFilmHeadContractError):
                validate_contract(ROOT, path)

    def test_rejects_changed_trainable_surface(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = mutated_contract(
                'trainable_parameter_prefixes = ["film_mlp_conv.", "film_mlp_transformer.", "masld_atac_head."]',
                'trainable_parameter_prefixes = ["film_mlp_conv.", "final_conv.", "masld_atac_head."]',
                temporary,
            )
            with self.assertRaises(CorgiFilmHeadContractError):
                validate_contract(ROOT, path)

    def test_rejects_outer_fold_in_inner_training(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = mutated_contract(
                "inner_training_folds = [2, 3, 4]",
                "inner_training_folds = [0, 3, 4]",
                temporary,
            )
            with self.assertRaises(CorgiFilmHeadContractError):
                validate_contract(ROOT, path)

    def test_rejects_context_specific_calibration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = mutated_contract(
                "arm_specific_calibration_forbidden = true",
                "arm_specific_calibration_forbidden = false",
                temporary,
            )
            with self.assertRaises(CorgiFilmHeadContractError):
                validate_contract(ROOT, path)

    def test_rejects_global_promotion_gate_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = mutated_contract(
                "global_promotion_gate_revision = false",
                "global_promotion_gate_revision = true",
                temporary,
            )
            with self.assertRaises(CorgiFilmHeadContractError):
                validate_contract(ROOT, path)


if __name__ == "__main__":
    unittest.main()
