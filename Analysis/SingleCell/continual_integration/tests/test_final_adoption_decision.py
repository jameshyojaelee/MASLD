import tempfile
import unittest
from pathlib import Path

from masld_cl.contracts import ContractError
from masld_cl.final_adoption_decision import _build_gates, parse_unittest_log


def artifact(value):
    return Path("unused"), value


class TestFinalAdoptionDecision(unittest.TestCase):
    def test_unittest_log_is_parsed(self):
        policy = {"software_tests": {"minimum_tests": 130, "maximum_skips": 1}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tests.log"
            path.write_text("Ran 134 tests in 8.2s\n\nOK (skipped=1)\n")
            observed = parse_unittest_log(path, policy)
        self.assertTrue(observed["passed"])
        self.assertEqual(observed["tests_run"], 134)
        self.assertEqual(observed["skipped"], 1)

    def test_failed_test_log_is_rejected(self):
        policy = {"software_tests": {"minimum_tests": 1, "maximum_skips": 1}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tests.log"
            path.write_text("Ran 1 test in 0.1s\n\nFAILED (failures=1)\n")
            with self.assertRaises(ContractError):
                parse_unittest_log(path, policy)

    def test_external_routing_is_a_required_failed_gate(self):
        artifacts = {
            "contract": artifact({"production_ready": True, "observed_contract": {
                "strict_reference": {"cells": 216957, "donors": 7, "libraries": 29},
                "primary_query": {"cells": 687559, "donors": 64, "libraries": 184, "controls": 9, "cases": 55},
                "analyzed": {"cells": 1232285, "donors": 102, "libraries": 273},
                "descriptive": {"cells": 1232318, "donors": 104, "libraries": 275},
            }}),
            "internal_confirmation": artifact({"confirmation_pass": True}),
            "selection_lock": artifact({"selection_frozen": True, "outcomes_unlocked": True}),
            "outcome_decision": artifact({"disease_preservation_pass": True, "decision": "continue_remaining_gates"}),
            "reference_decision": artifact({"decision": "retain_common_strict7"}),
            "full_atlas_projection": artifact({
                "contract": {"descriptive_cells": 1232318, "analyzed_cells": 1232285},
                "v32_coordinates_bitwise_identical": True,
                "secondary_and_descriptive_cells_used_for_fit": False,
            }),
            "production_expansion": artifact({
                "model_roster": ["all_lineage", "hepatocytes", "macrophages", "fibroblasts", "cholangiocytes", "t_cells"],
                "models": {name: {} for name in ["all_lineage", "hepatocytes", "macrophages", "fibroblasts", "cholangiocytes", "t_cells"]},
                "conditions_available_to_fit": False,
            }),
            "secondary_stress": artifact({"technical_stress_pass": True, "independent_secondary_confirmation_pass": False}),
            "program_firewall": artifact({"passed": True, "registry": {"programs": 117}}),
            "independent_external": artifact({
                "gates": {"all_required_control_alignment_gates": True, "all_required_disease_preservation_gates": True},
                "routing_audit": {"mean_donor_macro_f1": 0.276, "minimum": 0.7, "pass": False},
            }),
            "gpu_tolerance": artifact({
                "same_source_identity": True, "same_input_artifacts": True,
                "same_arguments_except_output": True, "same_cell_roster_and_order": True,
                "full_embedding": {"within_absolute_1e_6": True},
                "training_embedding": {"within_absolute_1e_6": True},
                "checkpoint": {"bitwise_identical_tensors": True},
                "observed_gpu_tolerance": 0.0,
            }),
        }
        gates = _build_gates(artifacts, {"passed": True})
        failed = [gate["gate_id"] for gate in gates if not gate["passed"]]
        self.assertEqual(failed, ["independent_external_routing"])


if __name__ == "__main__":
    unittest.main()
