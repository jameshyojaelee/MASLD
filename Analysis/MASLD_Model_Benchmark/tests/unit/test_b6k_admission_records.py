"""Contract tests for the frozen B6K hardware admission records.

These tests read configuration and execution manifests only. They never open a
probe fixture, a benchmark outcome, or any biological data.
"""

from __future__ import annotations

import hashlib
import json
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AUTHORIZATION = ROOT / "config" / "b6k_hardware_authorization_20260825.json"
ROUTING = ROOT / "config" / "b6k_routing_guidance_20260825.json"
PRESPECIFICATION = ROOT / "config" / "b6k_numeric_tolerance_prespecification_20260825.json"
DISPATCHER_POLICY = (
    ROOT / "config" / "campaigns" / "gpu_dependency_dispatcher_policy_20260825.json"
)
RUNTIMES = ROOT / "config" / "runtimes"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


class B6KAdmissionRecordTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authorization = load_json(AUTHORIZATION)
        cls.routing = load_json(ROUTING)
        cls.prespecification = load_json(PRESPECIFICATION)

    def test_records_are_record_only_and_not_consumed(self) -> None:
        for label, record in (
            ("authorization", self.authorization),
            ("routing", self.routing),
        ):
            self.assertIs(record["record_only"], True, label)
            self.assertIs(record["consumed_by_runtime"], False, label)

    def test_prespecification_identity_is_unchanged(self) -> None:
        self.assertEqual(
            digest(PRESPECIFICATION),
            "a26a6f3fe0d77ef12706038eaede32f222421154c0972e228b5238febeaa722f",
        )
        self.assertIs(
            self.prespecification["frozen_before_any_gpu_cell_executed"], True
        )
        self.assertIs(
            self.prespecification["decision_rule"]["no_post_hoc_threshold_revision"],
            True,
        )
        bound = self.authorization["probe_result"]["tolerance_prespecification"]
        self.assertEqual(bound["sha256"], digest(PRESPECIFICATION))
        self.assertEqual(self.routing["numeric_findings"]["prespecification_sha256"], bound["sha256"])

    def test_retired_blocker_clauses_match_the_runtime_files(self) -> None:
        for entry in self.authorization["supersedes_blocker_clauses"]:
            path = ROOT / entry["path"]
            self.assertTrue(path.is_file(), entry["path"])
            self.assertEqual(digest(path), entry["sha256_after"], entry["path"])
            self.assertNotEqual(entry["sha256_before"], entry["sha256_after"], entry["path"])
            runtime = tomllib.loads(path.read_text(encoding="utf-8"))
            clause = entry["retired_clause"]
            self.assertIn(clause, runtime["blockers_resolved"], entry["path"])
            for remaining in runtime["blockers"]:
                self.assertNotIn("B6K", remaining, entry["path"])

    def test_runtime_files_still_satisfy_the_registry_contract(self) -> None:
        for name in ("gpu_finetune", "gpu_inference"):
            path = RUNTIMES / f"{name}.toml"
            runtime = tomllib.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(runtime["schema_version"], "masld-bench-runtime-v1")
            self.assertEqual(runtime["runtime_id"], name)
            self.assertEqual(runtime["resource_profile"], "gpu_single")
            # Environment, container digest, modules, and command prefix stay
            # unresolved, so the registry contract still requires a blocker.
            self.assertIs(runtime["admission_blocking"], True, name)
            self.assertTrue(runtime["blockers"], name)
            self.assertEqual(runtime["environment_lock"], "UNRESOLVED", name)
            self.assertEqual(runtime["b6k_admission_status"], "admitted", name)
            self.assertEqual(runtime["blockers_resolved_on"], "2026-08-25", name)
            self.assertEqual(
                runtime["b6k_admission_record"],
                "config/b6k_hardware_authorization_20260825.json",
                name,
            )
            self.assertEqual(
                runtime["b6k_routing_record"],
                "config/b6k_routing_guidance_20260825.json",
                name,
            )
            for evidence in runtime["blockers_resolved_evidence"]:
                self.assertTrue((ROOT / evidence).is_dir(), evidence)

    def test_shared_resource_profile_still_defaults_to_l40s(self) -> None:
        resources = tomllib.loads(
            (ROOT / "config" / "resources.toml").read_text(encoding="utf-8")
        )
        self.assertEqual(resources["profiles"]["gpu_single"]["accelerator"], "l40s")

    def test_every_bound_execution_exists_and_hashes_as_recorded(self) -> None:
        bindings = [
            (cell["path"], cell["artifacts_sha256"]) for cell in self.authorization["probe_cells"]
        ]
        bindings.append(
            (
                self.authorization["numeric_comparison"]["path"],
                self.authorization["numeric_comparison"]["artifacts_sha256"],
            )
        )
        bindings.extend(
            (item["path"], item["artifacts_sha256"])
            for item in self.authorization["preflights"]
        )
        both = self.routing["pools_are_runtime_disjoint"][
            "only_runtime_that_executes_on_both_pools"
        ]
        bindings.append((both["environment_path"], both["artifacts_sha256"]))
        self.assertEqual(len(bindings), 10)
        for relative, expected in bindings:
            manifest = ROOT / relative / "ARTIFACTS.json"
            self.assertTrue(manifest.is_file(), relative)
            self.assertEqual(digest(manifest), expected, relative)

    def test_probe_cell_verdicts_match_the_probe_manifests(self) -> None:
        expected = {
            "cu124-l40s": ("pass", "minor_version_binary_compatibility"),
            "cu124-b6k": ("fail", "no_matching_sass"),
            "cu130-l40s": ("fail", "no_cuda_device_visible"),
            "cu130-b6k": ("pass", "native_sass"),
            "cu128-l40s": ("pass", "minor_version_binary_compatibility"),
            "cu128-b6k": ("pass", "native_sass"),
        }
        seen = set()
        for cell in self.authorization["probe_cells"]:
            cell_id = cell["cell_id"]
            seen.add(cell_id)
            status, execution_path = expected[cell_id]
            self.assertEqual(cell["status"], status, cell_id)
            self.assertEqual(cell["execution_path"], execution_path, cell_id)
            probe = load_json(ROOT / cell["path"] / "probe" / "probe.json")
            self.assertEqual(probe["cell_id"], cell_id)
            self.assertEqual(probe["status"], status, cell_id)
            self.assertEqual(probe["slurm_job_id"], cell["job_id"], cell_id)
            first = probe["probes"]["probe_1_import_architecture"]
            observed = first.get("execution_path", first.get("status"))
            self.assertEqual(observed, execution_path, cell_id)
            layer_guard = load_json(ROOT / cell["path"] / "layer_guard.json")
            self.assertIs(layer_guard["layer_guard_passed"], True, cell_id)
            self.assertEqual(layer_guard["torch_version"], cell["torch_version"], cell_id)
            self.assertIn(layer_guard["expected_layer"], layer_guard["torch_file"], cell_id)
        self.assertEqual(seen, set(expected))

    def test_all_six_probes_pass_on_b6k_under_both_admitted_runtimes(self) -> None:
        required = self.authorization["probe_result"]["required_probes"]
        self.assertEqual(len(required), 6)
        for cell_id in ("cu130-b6k", "cu128-b6k"):
            cell = next(c for c in self.authorization["probe_cells"] if c["cell_id"] == cell_id)
            probe = load_json(ROOT / cell["path"] / "probe" / "probe.json")
            for name in required:
                self.assertEqual(probe["probes"][name]["status"], "pass", f"{cell_id}/{name}")

    def test_numeric_findings_match_the_comparison_artifact(self) -> None:
        summary = load_json(
            ROOT / self.routing["numeric_findings"]["evidence"]["path"] / "comparison_summary.json"
        )
        recorded = self.routing["numeric_findings"]["results_on_the_prespecified_metric"]
        self.assertEqual(set(recorded), set(summary["results"]))
        for name, values in recorded.items():
            observed = summary["results"][name]
            self.assertEqual(values["max_rel_difference"], observed["max_rel_difference"], name)
            self.assertEqual(values["max_abs_difference"], observed["max_abs_difference"], name)
            self.assertEqual(values["verdict"], observed["verdict"], name)
            self.assertEqual(values["verdict"], "materially_different", name)
        self.assertEqual(
            summary["tolerance_prespecification_sha256"], digest(PRESPECIFICATION)
        )

    def test_caveat_arithmetic_is_reproducible_from_the_recorded_numbers(self) -> None:
        caveat = self.routing["numeric_findings"]["honest_caveat"]
        results = self.routing["numeric_findings"]["results_on_the_prespecified_metric"]
        architecture = results["architecture_effect__cu128_l40s_vs_cu128_b6k"]
        version = results["version_effect__cu124_l40s_vs_cu128_l40s"]

        self.assertIs(self.routing["numeric_findings"]["threshold_frozen_before_data"], True)
        self.assertIs(self.routing["numeric_findings"]["post_hoc_threshold_revision"], False)
        self.assertEqual(
            caveat["the_absolute_differences"]["architecture_effect"],
            architecture["max_abs_difference"],
        )
        self.assertEqual(
            caveat["the_absolute_differences"]["version_effect"],
            version["max_abs_difference"],
        )

        derived = caveat["derived_denominators"]
        self.assertAlmostEqual(
            derived["architecture_effect_worst_element_magnitude"],
            architecture["max_abs_difference"] / architecture["max_rel_difference"],
            places=8,
        )
        self.assertAlmostEqual(
            derived["version_effect_worst_element_magnitude"],
            version["max_abs_difference"] / version["max_rel_difference"],
            places=8,
        )

        ordering = caveat["what_remains_informative"]
        self.assertAlmostEqual(
            ordering["relative_ratio_version_over_architecture"],
            version["max_rel_difference"] / architecture["max_rel_difference"],
            places=3,
        )
        self.assertAlmostEqual(
            ordering["absolute_ratio_version_over_architecture"],
            version["max_abs_difference"] / architecture["max_abs_difference"],
            places=3,
        )
        # The version change must dominate the architecture change on both scales,
        # which is the only claim the metric still supports.
        self.assertGreater(ordering["relative_ratio_version_over_architecture"], 1.0)
        self.assertGreater(ordering["absolute_ratio_version_over_architecture"], 1.0)

    def test_eligibility_excludes_the_existing_l40s_rectangle_on_feasibility(self) -> None:
        excluded = self.routing["eligibility"]["not_eligible"]["named_exclusion"]
        self.assertEqual(excluded["ground"], "feasibility, not tolerance")
        self.assertIn("may NOT join", excluded["verdict"])
        self.assertIn("cu124", excluded["explanation"])
        self.assertIs(
            self.prespecification["constraints"][
                "existing_500_fit_rectangle_stays_l40s_homogeneous_by_default"
            ],
            True,
        )
        conditions = self.routing["eligibility"]["eligible"]["conditions"]
        self.assertEqual(len(conditions), 3)
        self.assertTrue(any("zero fits" in text for text in conditions))
        self.assertTrue(any("entirely on B6K" in text for text in conditions))
        self.assertTrue(any("single frozen runtime" in text for text in conditions))

    def test_layer_guard_requirement_is_recorded(self) -> None:
        guard = self.routing["mandatory_torch_layer_guard"]
        self.assertIn("torch.__file__", guard["requirement"])
        self.assertIn("PYTHONPATH", guard["why"])
        self.assertGreaterEqual(len(guard["assert_at_minimum"]), 4)
        reference = ROOT / guard["reference_implementation"]
        self.assertTrue(reference.is_file(), guard["reference_implementation"])
        self.assertIs(load_json(reference)["layer_guard_passed"], True)

    def test_dispatcher_gates_are_recorded_and_still_present_in_the_source(self) -> None:
        gates = self.authorization["dispatcher_gates_bypassed_by_manual_sbatch"]["gates"]
        self.assertEqual(
            {gate["gate"] for gate in gates},
            {"validate_item_route", "_submit", "load_policy", "audit_live_campaign"},
        )
        for gate in gates:
            source = (ROOT / gate["path"]).read_text(encoding="utf-8").splitlines()
            self.assertIn(f"def {gate['gate']}(", source[gate["line"] - 1], gate["gate"])

    def test_dispatcher_source_bindings_have_not_drifted(self) -> None:
        # This work must not perturb the running dispatcher. Every binding it
        # re-hashes must still match the policy.
        policy = load_json(DISPATCHER_POLICY)
        bindings = policy["source_bindings"]
        self.assertEqual(
            set(bindings),
            {
                "dependency_dispatcher",
                "legacy_queue_parser",
                "unit_test",
                "validation_wrapper",
                "runtime_wrapper",
            },
        )
        for label, record in bindings.items():
            path = ROOT / record["path"]
            self.assertTrue(path.is_file(), label)
            self.assertEqual(digest(path), record["sha256"], label)
        self.assertEqual(policy["routing"]["gres"], "gpu:l40s:1")
        self.assertEqual(policy["campaign_id"], "gpu_dependency_dispatcher_20260825")
        recorded = self.authorization["running_dispatcher_unaffected"]
        self.assertEqual(recorded["job_id"], "21100574")
        self.assertEqual(
            set(recorded["source_bindings_untouched"]),
            {record["path"] for record in bindings.values()},
        )


if __name__ == "__main__":
    unittest.main()
