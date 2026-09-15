from __future__ import annotations

import json
from pathlib import Path
import tomllib
import unittest

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.registry import load_dataset_manifest, load_task_spec


ROOT = Path(__file__).resolve().parents[2]
ACTIVATION = ROOT / "executions/gse267145-coordinate-activation-21070593"
ACTIVATION_SHA256 = "e5042eb43dddd84cda81c5aeed5d6dbeb0378ab1c942a0183e2ef1b41cd18b20"


class GSE267145CoordinateActivationContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        verify_frozen_tree(ACTIVATION)
        cls.artifacts = json.loads(
            (ACTIVATION / "ARTIFACTS.json").read_text(encoding="utf-8")
        )
        cls.coordinate = json.loads(
            (ACTIVATION / "coordinate_audit.json").read_text(encoding="utf-8")
        )
        cls.activation = json.loads(
            (ACTIVATION / "activation_contract.json").read_text(encoding="utf-8")
        )
        cls.task = load_task_spec(
            ROOT / "config/evaluation/paired_bulk_rna_h3k27ac_task.toml"
        )
        with (ROOT / "config/evaluation/cross_cohort_expansion.toml").open(
            "rb"
        ) as handle:
            cls.expansion = tomllib.load(handle)

    def test_frozen_receipt_activates_only_coordinate_independent_work(self) -> None:
        self.assertEqual(
            self.artifacts["metadata"]["task_id"], "paired_bulk_rna_h3k27ac"
        )
        self.assertTrue(
            self.artifacts["metadata"]
            ["task_specific_coordinate_independent_training_allowed"]
        )
        self.assertFalse(self.artifacts["metadata"]["sequence_extraction_allowed"])
        self.assertFalse(self.artifacts["metadata"]["champion_eligible"])

    def test_primary_sources_do_not_resolve_coordinate_semantics(self) -> None:
        self.assertFalse(self.coordinate["coordinate_semantics_resolved"])
        self.assertEqual(self.coordinate["coordinate_convention"], "unresolved")
        self.assertEqual(self.coordinate["supplement_liver_regions"], 14_348)
        self.assertTrue(
            self.coordinate["supplement_regions_are_exact_matrix_subset"]
        )
        for source in ("article", "geo_subseries_metadata", "workbook"):
            self.assertEqual(
                self.coordinate["direct_coordinate_claims"][source]
                ["zero_based_half_open"],
                [],
            )
            self.assertEqual(
                self.coordinate["direct_coordinate_claims"][source]
                ["one_based_inclusive"],
                [],
            )

    def test_task_is_participant_safe_and_outcome_firewalled(self) -> None:
        self.assertEqual(self.task.unit_of_inference, "participant")
        self.assertEqual(
            self.task.required_pairing_levels[0].value,
            "same_sample_different_aliquot",
        )
        self.assertEqual(self.activation["participants"], 99)
        self.assertEqual(
            self.activation["preprocessing"]
            ["all_learned_transforms_and_feature_selection_fit"],
            "outer_training_participants_only",
        )
        self.assertFalse(self.activation["single_cell_model_families_allowed"])
        self.assertFalse(
            self.activation["source_outcome_firewall"]
            ["model_selection_or_error_selector_allowed"]
        )

    def test_global_dataset_and_sequence_admission_remain_blocked(self) -> None:
        dataset = load_dataset_manifest(
            ROOT / "config/datasets/gse267145_znf469_human_liver.toml"
        )
        self.assertEqual(dataset.role.value, "blocked")
        self.assertEqual(dataset.status.value, "blocked")
        self.assertTrue(dataset.admission_blocking)
        self.assertFalse(self.activation["dataset_wide_or_sequence_activation"])
        self.assertIn("sequence_extraction", self.activation["blocked_actions"])
        self.assertIn(
            "derivative_weight_release_without_model_specific_terms_and_legal_review",
            self.activation["blocked_actions"],
        )

    def test_expansion_registry_pins_exact_activation(self) -> None:
        record = next(
            value
            for value in self.expansion["cohort_family"]
            if value["family_id"] == "gse267145_znf469_human_liver"
        )
        self.assertEqual(record["coordinate_activation_artifacts_sha256"], ACTIVATION_SHA256)
        self.assertTrue(record["task_specific_coordinate_independent_training_active"])
        self.assertFalse(record["coordinate_dependent_training_active"])


if __name__ == "__main__":
    unittest.main()
