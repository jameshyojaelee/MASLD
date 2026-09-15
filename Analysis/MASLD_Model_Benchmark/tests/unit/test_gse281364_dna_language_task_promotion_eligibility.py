from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.freeze_gse281364_dna_language_task_promotion_eligibility import (
    EligibilityAuditError,
    REQUIRED_COMPARATORS,
    eligibility_tsv,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/gse281364_dna_language_task_promotion_eligibility.json"


class GSE281364DNALanguageTaskPromotionEligibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_current_performance_blind_contract_passes(self) -> None:
        validate_config(self.config)
        self.assertFalse(self.config["scope"]["performance_read"])
        self.assertFalse(self.config["scope"]["prediction_values_read"])
        self.assertFalse(self.config["scope"]["reporter_outcomes_read"])
        self.assertFalse(self.config["scope"]["sealed_assets_read"])
        self.assertFalse(self.config["mandatory_comparator_gate"]["promotion_allowed"])

    def test_all_mandatory_task_native_comparators_remain_blockers(self) -> None:
        self.assertEqual(
            self.config["mandatory_comparator_gate"][
                "required_comparators_not_yet_complete_in_one_authoritative_same_universe_comparison"
            ],
            REQUIRED_COMPARATORS,
        )
        mutated = deepcopy(self.config)
        mutated["mandatory_comparator_gate"][
            "required_comparators_not_yet_complete_in_one_authoritative_same_universe_comparison"
        ].pop()
        with self.assertRaisesRegex(EligibilityAuditError, "comparator census"):
            validate_config(mutated)

    def test_restricted_nucleotide_transformer_cannot_be_open_champion(self) -> None:
        nt = self.config["models"]["nucleotide_transformer"]
        self.assertEqual(nt["terms"]["checkpoint_repository_declared_license"], "CC-BY-NC-SA-4.0")
        self.assertEqual(nt["action_eligibility"]["open_champion"], "ineligible_under_current_terms")
        mutated = deepcopy(self.config)
        mutated["models"]["nucleotide_transformer"]["action_eligibility"]["open_champion"] = "eligible"
        with self.assertRaisesRegex(EligibilityAuditError, "open champion"):
            validate_config(mutated)

    def test_hyenadna_release_fails_closed_without_notice_and_data_review(self) -> None:
        hyena = self.config["models"]["hyenadna"]
        self.assertTrue(hyena["terms"]["upstream_checkpoint_redistribution"].startswith("blocked_"))
        self.assertTrue(hyena["terms"]["project_derivative_head_redistribution"].startswith("blocked_"))
        mutated = deepcopy(self.config)
        mutated["models"]["hyenadna"]["terms"]["upstream_checkpoint_redistribution"] = "permitted"
        with self.assertRaisesRegex(EligibilityAuditError, "redistribution"):
            validate_config(mutated)

    def test_base_and_adapted_exposure_are_not_conflated(self) -> None:
        for model_id in ("hyenadna", "nucleotide_transformer"):
            contamination = self.config["models"][model_id]["contamination"]
            self.assertEqual(contamination["base_checkpoint_GSE281364_target_labels"], "target_label_unexposed")
            self.assertEqual(contamination["adapted_GSE281364_head"], "continual_seen_training_folds_only")
        mutated = deepcopy(self.config)
        mutated["models"]["hyenadna"]["contamination"]["adapted_GSE281364_head"] = "clean_declared"
        with self.assertRaisesRegex(EligibilityAuditError, "adapted-head exposure"):
            validate_config(mutated)

    def test_pinned_model_card_resolves_nt_token_count_conflict(self) -> None:
        pretraining = self.config["models"]["nucleotide_transformer"]["declared_pretraining"]
        self.assertEqual(pretraining["exact_pinned_model_card_pretraining_tokens"], 300_000_000_000)
        self.assertEqual(pretraining["prior_generic_manifest_pretraining_tokens"], 50_000_000_000)
        self.assertIn("do_not_cite_prior_50B_field", pretraining["token_count_disposition"])

    def test_performance_or_sealed_bindings_fail_closed(self) -> None:
        mutated = deepcopy(self.config)
        mutated["scope"]["performance_read"] = True
        with self.assertRaisesRegex(EligibilityAuditError, "firewall"):
            validate_config(mutated)

        mutated = deepcopy(self.config)
        mutated["frozen_authorities"]["source_admission_tree"]["path"] = "executions/model-cpu-train-605-21099008"
        with self.assertRaisesRegex(EligibilityAuditError, "forbidden"):
            validate_config(mutated)

    def test_eligibility_table_contains_actions_not_performance(self) -> None:
        receipt = {
            "models": {
                model_id: {"action_eligibility": record["action_eligibility"]}
                for model_id, record in self.config["models"].items()
            }
        }
        table = eligibility_tsv(receipt)
        self.assertIn("task_specific_tournament_promotion", table)
        self.assertIn("hyenadna", table)
        self.assertIn("nucleotide_transformer", table)
        self.assertNotIn("spearman", table.casefold())
        self.assertNotIn("macro_f1", table.casefold())


if __name__ == "__main__":
    unittest.main()
