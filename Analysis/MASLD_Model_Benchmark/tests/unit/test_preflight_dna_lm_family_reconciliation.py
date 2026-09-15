from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.preflight_dna_lm_family_reconciliation import (
    DnaLmReconciliationError,
    MODELS,
    preflight,
    validate_config,
    validate_embedding_contracts,
    validate_evaluations,
    validate_model_authorities,
    validate_shared_evaluator_code,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/dna_lm_family_reconciliation.json"


class DnaLmFamilyReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_family_contract_is_task_specific_and_fail_closed(self) -> None:
        validate_config(self.config)
        self.assertEqual(tuple(self.config["model_order"]), MODELS)
        task = self.config["task_contract"]
        self.assertEqual(task["dataset_id"], "gse281364")
        self.assertEqual(task["biological_donors"], 0)
        self.assertEqual(task["experimental_replicates_per_context"], 4)
        self.assertFalse(task["replicates_used_as_independent_donors"])
        self.assertFalse(task["external_evaluation"])
        self.assertFalse(task["champion_eligible"])
        self.assertEqual(
            task["context_policy"],
            "family_native_inputs_with_identical_downstream_heads",
        )
        self.assertEqual(
            task["head_contract"]["heads"],
            ["linear_ridge", "two_layer_gelu"],
        )

    def test_nucleotide_transformer_remains_restricted(self) -> None:
        models = self.config["models"]
        self.assertTrue(models["nucleotide_transformer"]["restricted_comparator"])
        self.assertFalse(
            models["nucleotide_transformer"][
                "open_champion_eligible_after_task_and_external_gates"
            ]
        )
        self.assertEqual(
            models["nucleotide_transformer"]["weight_license"],
            "CC-BY-NC-SA-4.0",
        )
        for model in ("dnabert2", "hyenadna", "caduceus"):
            self.assertFalse(models[model]["restricted_comparator"])

    def test_context_lengths_are_family_native_not_matched(self) -> None:
        models = self.config["models"]
        self.assertEqual(
            {models[model]["task_input_bp"] for model in MODELS[:3]}, {4096}
        )
        self.assertEqual(models["caduceus"]["task_input_bp"], 131072)
        self.assertIn(
            "head_comparison_is_identical_but_representation_context_is_not_matched",
            self.config["task_contract"]["context_comparability_warning"],
        )

    def test_claim_expansion_is_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["claim_boundaries"]["universal_model_claim_supported"] = True
        with self.assertRaisesRegex(DnaLmReconciliationError, "claim boundary"):
            validate_config(changed)

    def test_exact_model_and_evaluator_authorities_reconcile(self) -> None:
        validate_shared_evaluator_code(ROOT, self.config)
        census = validate_model_authorities(ROOT, self.config)
        self.assertEqual(tuple(row["model_id"] for row in census), MODELS)
        self.assertEqual(
            [row["model_id"] for row in census if row["restricted_comparator"]],
            ["nucleotide_transformer"],
        )
        embedding = validate_embedding_contracts(ROOT, self.config)
        self.assertEqual(embedding["projection_width"], 256)
        self.assertEqual(embedding["head_input_width"], 1024)
        self.assertFalse(
            embedding["reporter_outcomes_read_during_embedding_or_projection"]
        )
        self.assertFalse(embedding["sealed_labels_read"])
        promotions, audit = validate_evaluations(ROOT, self.config)
        self.assertEqual(tuple(row["model_id"] for row in promotions), MODELS)
        self.assertTrue(audit["baseline_rows_identical_between_evaluators"])
        self.assertTrue(audit["identical_head_functions_reused"])
        self.assertTrue(
            audit["all_model_gains_over_allele_identity_ridge_are_negative"]
        )
        self.assertEqual(audit["three_seed_screening_recommended"], [])
        self.assertEqual(audit["family_promotion"], "none")

    def test_integrated_receipt_preserves_negative_result(self) -> None:
        receipt = preflight(root=ROOT, config_path=CONFIG_PATH)
        self.assertEqual(
            receipt["status"],
            "pass_reconciled_development_smoke_no_model_promoted",
        )
        self.assertFalse(receipt["representation_context_matched"])
        self.assertTrue(receipt["downstream_head_contract_matched"])
        self.assertTrue(receipt["development_summary_metrics_read"])
        self.assertFalse(receipt["raw_reporter_outcomes_parsed"])
        self.assertFalse(receipt["sealed_assets_read"])
        self.assertFalse(receipt["supervised_eQTL_or_ieQTL_head_fit"])
        self.assertFalse(receipt["external_evaluation"])
        self.assertFalse(receipt["champion_eligible"])
        self.assertFalse(receipt["universal_claim_supported"])
        self.assertEqual(
            receipt["evaluation_audit"]["best_model_by_mean_spearman"],
            "nucleotide_transformer",
        )
        self.assertTrue(receipt["evaluation_audit"]["best_model_is_restricted"])
        self.assertEqual(receipt["evaluation_audit"]["family_promotion"], "none")


if __name__ == "__main__":
    unittest.main()
