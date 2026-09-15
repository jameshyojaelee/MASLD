from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.preflight_local_variant_family_reconciliation import (
    LocalVariantReconciliationError,
    MODELS,
    STATUS,
    TASKS,
    preflight,
    validate_authorities,
    validate_capability_and_task,
    validate_config,
    validate_frozen_evidence,
    validate_registry,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/local_variant_family_reconciliation.json"


class LocalVariantFamilyReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_roster_is_partitioned_by_native_endpoint(self) -> None:
        validate_config(self.config)
        self.assertEqual(tuple(self.config["model_order"]), MODELS)
        self.assertEqual(tuple(self.config["native_task_order"]), TASKS)
        tasks = self.config["native_tasks"]
        self.assertEqual(
            tasks["local_profile_accessibility"]["models"], list(MODELS[:4])
        )
        self.assertEqual(tasks["mpra_reporter_activity"]["models"], ["mpralegnet"])
        self.assertEqual(
            tasks["enhancer_gene_link"]["models"], ["abc", "re2g", "epinformer"]
        )
        self.assertEqual(
            tasks["observed_atac_native_link"]["models"], ["abc", "re2g"]
        )
        self.assertFalse(
            tasks["mpra_reporter_activity"]["replicates_used_as_independent_donors"]
        )

    def test_license_exposure_and_assay_boundaries_fail_closed(self) -> None:
        models = self.config["models"]
        self.assertEqual(models["gkmsvm"]["code_license"], "GPL-3.0")
        self.assertEqual(models["deltasvm"]["code_license"], "GPL-3.0")
        self.assertEqual(models["abc"]["exposure"], "reference_only")
        self.assertEqual(models["re2g"]["exposure"], "clean_declared")
        self.assertTrue(models["abc"]["query_time_observed_atac"])
        self.assertTrue(models["re2g"]["query_time_observed_atac"])
        self.assertTrue(models["epinformer"]["runtime"].startswith("blocked"))
        self.assertIn("assay_QC_diagnostic_only", models["chrombpnet"]["full_head_role"])
        self.assertTrue(all(not row["signed_gene_effect"] for row in models.values()))
        self.assertTrue(
            all(not row["open_champion_eligible_now"] for row in models.values())
        )

    def test_hashes_registry_and_frozen_receipts_reconcile(self) -> None:
        authorities = validate_authorities(ROOT, self.config)
        self.assertEqual(len(authorities), 38)
        registry = validate_registry(ROOT, self.config)
        self.assertEqual(len(registry), 8)
        validate_capability_and_task(ROOT, self.config)
        self.assertEqual(
            [row["model_id"] for row in registry if row["registry_status"] == "blocked_terms"],
            ["epinformer"],
        )
        evidence = validate_frozen_evidence(ROOT, self.config)
        self.assertEqual(evidence["local_profile_evaluation_groups"], 7)
        self.assertFalse(evidence["local_profile_matrix_rankable"])
        self.assertEqual(evidence["bpnet_frozen_logical_tasks"], 8)
        self.assertEqual(evidence["sequence_cnn_frozen_logical_tasks"], 8)
        self.assertEqual(evidence["mpralegnet_biological_donors"], 0)

    def test_claim_expansion_and_premature_screen_are_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["family_boundaries"]["champion_claim_allowed"] = True
        with self.assertRaises(LocalVariantReconciliationError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["native_tasks"]["enhancer_gene_link"][
            "additional_development_screen_justified_now"
        ] = True
        with self.assertRaises(LocalVariantReconciliationError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["models"]["epinformer"]["runtime"] = "ready"
        with self.assertRaises(LocalVariantReconciliationError):
            validate_config(changed)

    def test_integrated_receipt_permits_only_conditional_next_screens(self) -> None:
        receipt = preflight(root=ROOT, config_path=CONFIG_PATH)
        self.assertEqual(receipt["status"], STATUS)
        self.assertEqual(receipt["model_rows"], 10)
        self.assertEqual(receipt["native_task_rows"], 5)
        self.assertEqual(receipt["authority_rows"], 38)
        self.assertEqual(receipt["additional_development_screens_allowed_now"], [])
        self.assertEqual(len(receipt["conditional_next_screens"]), 3)
        self.assertFalse(receipt["active_campaign_outputs_bound_or_read"])
        self.assertFalse(receipt["sealed_assets_read"])
        self.assertFalse(receipt["test_outcomes_read"])
        self.assertFalse(receipt["new_training_scoring_or_model_import_performed"])
        self.assertFalse(receipt["primary_signed_eQTL_or_ieQTL_candidate_present"])
        self.assertFalse(receipt["champion_claim_allowed"])
        self.assertFalse(receipt["universal_model_claim_allowed"])
        self.assertFalse(receipt["clinical_claim_allowed"])
        self.assertFalse(receipt["bundled_cpu_validation_ask"]["submitted"])


if __name__ == "__main__":
    unittest.main()
