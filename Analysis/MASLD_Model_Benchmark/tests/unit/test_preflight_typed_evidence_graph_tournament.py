from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.preflight_typed_evidence_graph_tournament import (
    BLOCKED_STATUS,
    GraphTournamentPreflightError,
    MANDATORY_BASELINES,
    MODEL_ORDER,
    preflight,
    validate_authorities,
    validate_config,
    validate_dataset_dispositions,
    validate_model_registries,
    validate_task_split_and_promotion,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/typed_evidence_graph_tournament_preflight.json"


class TypedEvidenceGraphTournamentPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_roster_and_task_are_blocked_fail_closed(self) -> None:
        validate_config(self.config)
        self.assertEqual(tuple(self.config["tournament_model_order"]), MODEL_ORDER)
        self.assertEqual(
            tuple(
                model_id
                for model_id in MODEL_ORDER
                if self.config["models"][model_id]["mandatory"]
            ),
            MANDATORY_BASELINES,
        )
        self.assertIsNone(self.config["strongest_currently_executable_task"])
        self.assertEqual(
            self.config["first_activation_candidate"]["candidate_id"],
            "current_coloc_gwas_leave_one_study_family_out",
        )
        self.assertEqual(
            self.config["first_activation_candidate"]["status"], "blocked"
        )
        self.assertTrue(
            all(
                not row["currently_executable"]
                for row in self.config["candidate_source_audit"].values()
            )
        )
        self.assertFalse(
            self.config["authority_lineage_audit"][
                "current_registry_matches_recorded_registry"
            ]
        )

    def test_topology_boundaries_are_not_overclaimed(self) -> None:
        models = self.config["models"]
        self.assertIn("common_support_only", models["node2vec"]["evaluation_support"])
        self.assertIn("common_support_only", models["bionic"]["evaluation_support"])
        self.assertIn(
            "source_independent_root_features",
            models["graphsage"]["evaluation_support"],
        )
        self.assertEqual(
            models["rgcn"]["evaluation_support"],
            "cannot_directly_complete_an_unseen_relation_type",
        )
        self.assertEqual(
            models["hgt"]["evaluation_support"],
            "cannot_score_untrained_unseen_node_or_edge_types",
        )
        self.assertTrue(
            models["typed_transformer"]["current_status"].startswith(
                "blocked_no_implementation_identity"
            )
        )
        self.assertIn(
            "identical_source_removal", models["late_fusion"]["evaluation_support"]
        )

    def test_authorities_registries_and_data_roles_are_consistent(self) -> None:
        authorities = validate_authorities(ROOT, self.config)
        self.assertEqual(len(authorities), 38)
        model_rows = validate_model_registries(ROOT, self.config)
        self.assertEqual(tuple(row["model_id"] for row in model_rows), MODEL_ORDER)
        self.assertTrue(all(row["admission_blocking"] for row in model_rows))
        validate_task_split_and_promotion(ROOT, self.config)
        validate_dataset_dispositions(ROOT, self.config)

    def test_config_rejects_activation_without_frozen_artifacts(self) -> None:
        changed = deepcopy(self.config)
        changed["strongest_currently_executable_task"] = {
            "candidate_id": "unfrozen_graph"
        }
        with self.assertRaises(GraphTournamentPreflightError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["candidate_source_audit"]["current_coloc_gwas"][
            "currently_executable"
        ] = True
        with self.assertRaises(GraphTournamentPreflightError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["required_artifact_bindings"]["graph_snapshot"] = "unfrozen.tsv"
        with self.assertRaises(GraphTournamentPreflightError):
            validate_config(changed)

    def test_config_rejects_missing_degree_or_observability_guards(self) -> None:
        changed = deepcopy(self.config)
        changed["tournament_model_order"].remove("graph_degree")
        with self.assertRaises(GraphTournamentPreflightError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["activation_gates"] = [
            row
            for row in changed["activation_gates"]
            if row["gate_id"]
            != "degree_source_count_and_observability_matched_controls"
        ]
        with self.assertRaises(GraphTournamentPreflightError):
            validate_config(changed)

    def test_integrated_receipt_contains_no_outcomes_or_execution(self) -> None:
        receipt = preflight(root=ROOT, config_path=CONFIG_PATH)
        self.assertEqual(receipt["status"], BLOCKED_STATUS)
        self.assertIsNone(receipt["strongest_currently_executable_task"])
        self.assertEqual(receipt["model_census_rows"], 12)
        self.assertEqual(receipt["mandatory_baseline_rows"], 5)
        self.assertEqual(receipt["candidate_source_rows"], 3)
        self.assertEqual(receipt["activation_gate_rows"], 17)
        self.assertEqual(receipt["required_gates_unresolved"], 16)
        self.assertFalse(receipt["sealed_graph_source_registered"])
        self.assertFalse(receipt["graph_snapshot_bound"])
        self.assertFalse(receipt["wholly_removed_source_bound"])
        self.assertFalse(receipt["explicit_tested_negative_denominator_bound"])
        self.assertFalse(receipt["training_or_scoring_allowed"])
        self.assertFalse(receipt["champion_claim_allowed"])
        self.assertFalse(receipt["raw_or_summary_labels_read"])
        self.assertFalse(receipt["sealed_assets_read"])
        self.assertFalse(receipt["model_import_fit_predict_or_score"])
        self.assertFalse(receipt["network_download_or_install"])
        self.assertFalse(receipt["bundled_cpu_validation_ask"]["submitted"])


if __name__ == "__main__":
    unittest.main()
