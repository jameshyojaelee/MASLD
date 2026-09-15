"""Contract tests binding the protein_transport TaskSpec to its activation evidence."""

from __future__ import annotations

import json
import unittest
from hashlib import sha256
from pathlib import Path

from masld_bench.registry import load_task_spec

ROOT = Path(__file__).resolve().parents[2]
TASK = ROOT / "config" / "evaluation" / "protein_transport_task.toml"
GATE = ROOT / "config" / "evaluation" / "protein_transport_promotion_gate.json"
POWER = ROOT / "config" / "campaigns" / "protein_transport_power_position_20260825.json"
REGISTRATION = ROOT / "config" / "campaigns" / "protein_transport_task_registration_20260825.json"
ACTIVATION_ARTIFACTS_SHA256 = "8488568d468f4fe15dff152cd1529f069b69b9ef6f27708a518fdfc757b681a7"


class ProteinTransportTaskSpec(unittest.TestCase):
    def setUp(self) -> None:
        self.task = load_task_spec(TASK)
        self.gate = json.loads(GATE.read_text())
        self.power = json.loads(POWER.read_text())
        self.registration = json.loads(REGISTRATION.read_text())

    def test_schema_version_and_identity(self) -> None:
        self.assertEqual(self.task.schema_version, "masld-bench-task-v1")
        self.assertEqual(self.task.task_id, "protein_transport")
        self.assertEqual(self.task.unit_of_inference, "participant")

    def test_promotion_gate_hash_binds(self) -> None:
        self.assertEqual(
            self.task.promotion_gate_config_sha256, sha256(GATE.read_bytes()).hexdigest()
        )
        self.assertEqual(self.task.promotion_gate_id, self.gate["promotion_gate_id"])

    def test_no_champion_is_reachable(self) -> None:
        self.assertFalse(self.gate["champion_eligible"])
        self.assertEqual(self.gate["claim_mode"], "development_only")
        self.assertEqual(self.gate["external_holdout_dataset_ids"], [])
        self.assertEqual(
            self.registration["champion_claim"]["state"], "structurally_impossible_not_pending"
        )
        self.assertEqual(self.registration["champion_claim"]["eligible_external_protein_cohorts"], 0)

    def test_baselines_match_the_registration_record(self) -> None:
        self.assertEqual(
            sorted(self.task.baseline_model_ids),
            sorted(self.registration["task_contract"]["mandatory_baselines"]),
        )

    def test_pairing_level_is_same_donor_different_tissue(self) -> None:
        self.assertEqual(
            [level.value for level in self.task.required_pairing_levels],
            ["same_donor_different_tissue"],
        )

    def test_the_only_dataset_is_pxd051911_and_partitions_are_disjoint(self) -> None:
        self.assertEqual(list(self.task.datasets_train), ["pxd051911"])
        self.assertEqual(list(self.task.datasets_development), [])
        self.assertEqual(list(self.task.datasets_sealed), [])

    def test_missingness_is_never_zero(self) -> None:
        self.assertIn("never_encode_missing_as_zero", self.task.missingness_policy)
        self.assertIn(
            "encoding_missing_protein_or_fat_pct_as_zero", self.registration["prohibited"]
        )

    def test_no_rna_protein_pair_is_assumed(self) -> None:
        self.assertFalse(self.registration["task_contract"]["rna_protein_pair_assumed"])
        self.assertFalse(self.task.evaluator_parameters["rna_protein_pair_assumed"])
        self.assertTrue(
            any("RNA-protein" in gate or "RNA measurement" in gate for gate in self.task.admission_gates)
        )

    def test_resampling_unit_is_the_participant(self) -> None:
        self.assertEqual(list(self.task.resampling_units), ["participant"])
        self.assertEqual(self.task.bootstrap_replicates, 10000)


class ProteinTransportPowerPosition(unittest.TestCase):
    def setUp(self) -> None:
        self.power = json.loads(POWER.read_text())
        self.task = load_task_spec(TASK)

    def test_prespecified_before_any_fitting(self) -> None:
        self.assertEqual(self.power["status"], "PRESPECIFIED_BEFORE_ANY_FITTING")
        self.assertEqual(self.power["models_fit"], [])
        self.assertEqual(self.power["models_scored"], [])
        self.assertFalse(self.power["metrics_calculated"])
        self.assertFalse(self.power["sealed_outcomes_read"])

    def test_binds_the_activation_artifact_by_hash(self) -> None:
        bound = {item["path"]: item["sha256"] for item in self.power["bound_evidence"]}
        self.assertEqual(
            bound["executions/pxd051911-activation-readiness-21109461/ARTIFACTS.json"],
            ACTIVATION_ARTIFACTS_SHA256,
        )
        for path, digest in bound.items():
            self.assertEqual(sha256((ROOT / path).read_bytes()).hexdigest(), digest, path)

    def test_split_agrees_with_the_task_spec(self) -> None:
        self.assertEqual(self.power["split"]["outer_folds"], self.task.evaluator_parameters["outer_folds"])
        self.assertEqual(self.power["split"]["seed"], self.task.evaluator_parameters["fold_seed"])
        self.assertEqual(self.power["split"]["fold_sizes"], [10, 17, 10, 14, 7])
        self.assertEqual(sum(self.power["split"]["fold_sizes"]), 58)
        self.assertFalse(self.power["split"]["stratified"])

    def test_primary_endpoint_is_defined_in_every_fold(self) -> None:
        primary = self.power["primary_endpoint"]["name"]
        self.assertEqual(primary, self.task.evaluator_parameters["primary_endpoint"])
        endpoint = self.power["endpoints"][primary]
        # The primary is a rank endpoint. Spearman needs a non-constant held-out vector,
        # not every level in every fold, so that is the property under test. NAS levels
        # 6, 7 and 8 hold one participant each and are absent from most folds by
        # construction; asserting no empty level here would be the wrong criterion.
        self.assertEqual(endpoint["endpoint_kind"], "rank")
        self.assertEqual(endpoint["definedness_criterion"], "non_constant_within_every_held_out_fold")
        self.assertTrue(endpoint["non_constant_in_every_fold"])
        self.assertTrue(all(count >= 2 for count in endpoint["distinct_values_per_fold"]))
        self.assertTrue(endpoint["defined_pooled_out_of_fold"])

    def test_saf_is_admissible_only_pooled_because_a_level_is_empty_in_a_fold(self) -> None:
        saf = self.power["endpoints"]["saf_diagnosis_three_state"]
        self.assertEqual(saf["endpoint_kind"], "classification")
        self.assertFalse(saf["every_level_present_in_every_fold"])
        self.assertIn("No_MASLD", saf["levels_with_an_empty_outer_fold"])
        self.assertTrue(saf["defined_pooled_out_of_fold"])

    def test_scoring_is_pooled_out_of_fold_never_averaged_per_fold(self) -> None:
        self.assertEqual(
            self.task.evaluator_parameters["scoring_pooling"],
            "pooled_out_of_fold_never_average_of_per_fold",
        )
        self.assertEqual(
            self.power["primary_endpoint"]["pooling"],
            "pooled out-of-fold over all 58 participants, never an average of per-fold values",
        )

    def test_underpowered_endpoints_are_named_and_excluded(self) -> None:
        excluded = {item["name"] for item in self.power["endpoints_prespecified_as_not_used"]}
        self.assertIn("kleiner_fibrosis_grade, four state", excluded)
        self.assertIn("ballooning, three state", excluded)

    def test_a_negative_is_an_acceptable_outcome(self) -> None:
        policy = self.power["negative_result_policy"]
        self.assertTrue(policy["a_well_characterised_negative_is_an_acceptable_outcome"])
        self.assertIn("switching the primary endpoint after seeing any result",
                      policy["prohibited_responses_to_a_negative"])

    def test_prevalence_is_prohibited_as_a_baseline(self) -> None:
        floor = self.power["metric_floor_binding"]
        self.assertTrue(floor["prevalence_as_baseline_prohibited"])
        self.assertTrue(floor["both_references_required"])
        for field in ("n", "positive_count", "negative_count",
                      "distinct_score_values_emitted_by_the_scorer"):
            self.assertIn(field, floor["receipt_must_record"])
        self.assertTrue(self.task.evaluator_parameters["prevalence_as_baseline_prohibited"])

    def test_cross_tissue_scopes_are_separate_and_sized(self) -> None:
        scopes = {s["scope_id"]: s for s in self.power["cross_tissue_scopes"]}
        self.assertEqual(scopes["S-LIVER"]["n"], 58)
        self.assertEqual(scopes["S-SCWAT"]["n"], 56)
        self.assertEqual(scopes["S-OWAT"]["n"], 25)
        self.assertTrue(scopes["S-OWAT"]["thin"])
        for key in ("S-SCWAT", "S-OWAT"):
            self.assertEqual(scopes[key]["pairing_level"], "same_donor_different_tissue")
            self.assertTrue(scopes[key]["never_same_sample"])


if __name__ == "__main__":
    unittest.main()
