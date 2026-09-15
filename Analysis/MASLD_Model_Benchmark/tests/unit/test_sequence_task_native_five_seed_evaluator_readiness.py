from __future__ import annotations

from inspect import signature
import importlib.util
import json
from pathlib import Path
import unittest

from masld_bench.evaluators.stats import two_way_donor_block_bootstrap


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts/audit_sequence_task_native_five_seed_evaluator_readiness.py"
SPEC = importlib.util.spec_from_file_location("five_seed_evaluator_readiness", SOURCE)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class FiveSeedRectangularEvaluatorReadinessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = AUDIT.load_contract(
            ROOT / "config/evaluation/sequence_task_native_five_seed_rectangular_evaluator_20260825.json"
        )
        cls.keys = AUDIT.expected_keys(cls.contract)

    @staticmethod
    def artifact(key: tuple[str, str, int, int]) -> dict[str, object]:
        return {
            "artifact_path": f"/frozen/{key[0]}/{key[1]}/{key[2]}/{key[3]}",
            "artifacts_sha256": "a" * 64,
        }

    def test_exact_rectangle_and_chrombpnet_view_expansion(self) -> None:
        self.assertEqual(len(self.keys), 500)
        compatible = {key: self.artifact(key) for key in self.keys}
        rows, status = AUDIT.classify_rectangle(self.keys, compatible, {}, {})
        self.assertEqual(status, "ready_full_rectangle_locked")
        views = AUDIT.expand_views(rows, {
            "bpnet": ["bpnet"],
            "chrombpnet": ["chrombpnet_full", "chrombpnet_nobias"],
            "sequence_cnn_control": ["sequence_cnn_control"],
            "sequence_transformer_control": ["sequence_transformer_control"],
        })
        summary = AUDIT.readiness_summary(
            fit_rows=rows,
            view_rows=views,
            gate_status=status,
        )
        self.assertEqual(len(views), 625)
        self.assertEqual(summary["complete_units_by_candidate_view"]["chrombpnet_full"], 125)
        self.assertEqual(summary["complete_units_by_candidate_view"]["chrombpnet_nobias"], 125)
        self.assertTrue(summary["outcome_access_authorized"])
        self.assertTrue(summary["development_cross_model_ranking_authorized"])

    def test_one_missing_unit_forbids_outcome_access_and_any_partial_ranking(self) -> None:
        compatible = {key: self.artifact(key) for key in self.keys[:-1]}
        rows, status = AUDIT.classify_rectangle(self.keys, compatible, {}, {})
        views = AUDIT.expand_views(rows, {
            "bpnet": ["bpnet"],
            "chrombpnet": ["chrombpnet_full", "chrombpnet_nobias"],
            "sequence_cnn_control": ["sequence_cnn_control"],
            "sequence_transformer_control": ["sequence_transformer_control"],
        })
        summary = AUDIT.readiness_summary(
            fit_rows=rows,
            view_rows=views,
            gate_status=status,
        )
        self.assertEqual(status, "waiting_for_full_rectangle")
        self.assertFalse(summary["outcome_access_authorized"])
        self.assertFalse(summary["development_cross_model_ranking_authorized"])
        self.assertFalse(summary["partial_cross_model_ranking_authorized"])

    def test_terminal_failure_and_incompatibility_fail_closed(self) -> None:
        missing = self.keys[-1]
        compatible = {key: self.artifact(key) for key in self.keys[:-1]}
        receipts = {missing: [{"status": "failed_missing_artifact"}]}
        _rows, status = AUDIT.classify_rectangle(self.keys, compatible, receipts, {})
        self.assertEqual(status, "blocked_terminal_failure_requires_new_immutable_retry")

        receipts = {missing: [{"status": "passed"}]}
        _rows, status = AUDIT.classify_rectangle(self.keys, compatible, receipts, {})
        self.assertEqual(status, "blocked_incompatible_artifact")
        _rows, status = AUDIT.classify_rectangle(
            self.keys,
            compatible,
            {},
            {missing: [{"reason": "required_fit_or_prediction_artifact_missing"}]},
        )
        self.assertEqual(status, "blocked_incompatible_artifact")

    def test_contract_locks_five_seed_ensemble_and_two_way_uncertainty(self) -> None:
        ensemble = self.contract["five_seed_ensemble"]
        uncertainty = self.contract["uncertainty"]
        self.assertIn("exactly_five", ensemble["profile_probability_aggregation"])
        self.assertFalse(ensemble["seeds_are_biological_replicates"])
        self.assertEqual(uncertainty["bootstrap_replicates"], 10_000)
        self.assertEqual(uncertainty["method"], "paired_two_way_donor_genomic_block_bootstrap")
        self.assertFalse(uncertainty["cells_windows_and_seeds_are_independent_units"])
        self.assertEqual(signature(two_way_donor_block_bootstrap).parameters["n_resamples"].default, 10_000)

    def test_contract_forbids_test_sealed_and_partial_evaluation(self) -> None:
        firewall = self.contract["firewall"]
        gate = self.contract["readiness_gate"]
        self.assertFalse(firewall["test_role_authorized"])
        self.assertFalse(firewall["sealed_data_authorized"])
        self.assertFalse(gate["model_specific_partial_ranking_allowed"])
        self.assertFalse(gate["missing_or_failed_units_may_be_zero_filled"])
        self.assertFalse(gate["missing_or_failed_units_may_be_dropped"])


if __name__ == "__main__":
    unittest.main()
