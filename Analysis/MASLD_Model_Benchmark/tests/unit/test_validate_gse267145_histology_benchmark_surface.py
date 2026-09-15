from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.validate_gse267145_histology_benchmark_surface import (
    HistologySurfaceError,
    build_endpoint_safe_inner_contract,
    read_tsv,
    validate_surface,
)


ROOT = Path(__file__).resolve().parents[2]
SURFACE = ROOT / "config/evaluation/gse267145_histology_benchmark_surface.json"


def surface() -> dict:
    return json.loads(SURFACE.read_text(encoding="utf-8"))


class GSE267145HistologyBenchmarkSurfaceTests(unittest.TestCase):
    def test_frozen_surface_has_three_distinct_lanes_and_no_fit_gate(self) -> None:
        result = validate_surface(surface(), ROOT)
        self.assertEqual(result["status"], "passed_design_only_no_fit")
        self.assertEqual(result["lanes"], ["rna_only", "h3k27ac_only", "observed_pair_multimodal"])
        self.assertEqual(result["baseline_count"], 11)
        self.assertFalse(result["production_fit_authorized"])
        self.assertFalse(result["model_fitted"])
        self.assertTrue(result["one_standard_error_selection"])
        self.assertTrue(result["endpoint_safe_inner_splits"])

    def test_outer_test_outcomes_and_recorded_sex_are_fit_worker_blocked(self) -> None:
        value = surface()
        blocked = value["firewall"]["fit_worker_never_receives"]
        self.assertIn("outer_test_endpoint_values", blocked)
        self.assertIn("recorded_sex", blocked)
        self.assertEqual(value["firewall"]["participant_age_state"], "structurally_missing")
        self.assertFalse(value["firewall"]["single_cell_methods_allowed"])

    def test_global_feature_selection_and_coordinate_inference_fail_closed(self) -> None:
        value = surface()
        value["firewall"]["h3_coordinate_inference_allowed"] = True
        with self.assertRaises(HistologySurfaceError):
            validate_surface(value, ROOT)
        value = surface()
        value["firewall"]["published_outcome_selected_features_allowed"] = True
        with self.assertRaises(HistologySurfaceError):
            validate_surface(value, ROOT)

    def test_lane_modality_crossover_fails_closed(self) -> None:
        value = surface()
        value["lanes"]["rna_only"]["query_modalities"].append("h3k27ac_cutrun")
        with self.assertRaises(HistologySurfaceError):
            validate_surface(value, ROOT)

    def test_male_macro_f1_cannot_be_claimed_with_only_nash_males(self) -> None:
        value = surface()
        audit = value["independent_evaluator"]["sex_error_audit"]
        self.assertEqual(audit["male_census"]["stage3"], {"NOR": 0, "NAFL": 0, "NASH": 14})
        self.assertEqual(audit["male_stage3_macro_f1"], "not_applicable_missing_NOR_and_NAFL")
        changed = deepcopy(value)
        changed["independent_evaluator"]["sex_error_audit"]["male_stage3_macro_f1"] = "report"
        with self.assertRaises(HistologySurfaceError):
            validate_surface(changed, ROOT)

    def test_surface_binds_exact_molecular_and_evaluator_artifacts(self) -> None:
        value = surface()
        value["inputs"]["molecular_fixture"]["artifacts_sha256"] = "0" * 64
        with self.assertRaises(HistologySurfaceError):
            validate_surface(value, ROOT)

    def test_every_inner_training_partition_retains_sparse_fibrosis_classes(self) -> None:
        _, endpoints = read_tsv(
            ROOT / "executions/model-data-061-21079623/activation/participant_endpoints.tsv"
        )
        roster, census = build_endpoint_safe_inner_contract(endpoints)
        self.assertEqual(len(roster), 396)
        self.assertTrue(census["every_inner_training_endpoint_safe"])
        self.assertGreaterEqual(census["minimum_inner_training_fibrosis3"], 1)
        self.assertGreaterEqual(census["minimum_inner_training_threshold_class_count"], 1)
        self.assertNotIn("recorded_sex", census["algorithm_uses"])
        self.assertNotIn("source_stage5", census["algorithm_uses"])

    def test_exact_tie_only_rule_and_missing_regression_head_fail_closed(self) -> None:
        value = surface()
        value["resampling"]["inner_hyperparameter_selection"]["exact_tie_rule_only"] = True
        with self.assertRaises(HistologySurfaceError):
            validate_surface(value, ROOT)
        value = surface()
        del value["shared_task_native_heads"]["fibrosis_exact_regression"]
        with self.assertRaises(HistologySurfaceError):
            validate_surface(value, ROOT)


if __name__ == "__main__":
    unittest.main()
