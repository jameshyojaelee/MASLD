from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import tomllib
import unittest

from masld_bench.registry import load_task_spec


ROOT = Path(__file__).resolve().parents[2]
EVALUATION = ROOT / "config/evaluation"


def sha256_file(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class MicroarrayTransferContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (EVALUATION / "microarray_transfer_preprocessing.toml").open("rb") as handle:
            cls.preprocessing = tomllib.load(handle)
        with (EVALUATION / "cross_cohort_expansion.toml").open("rb") as handle:
            cls.expansion = tomllib.load(handle)
        with (EVALUATION / "microarray_transfer_label_semantics.toml").open("rb") as handle:
            cls.labels = tomllib.load(handle)
        cls.fibrosis = load_task_spec(
            EVALUATION / "gse49541_fibrosis_transfer_task.toml"
        )
        cls.baseline = load_task_spec(
            EVALUATION / "gse83452_baseline_nash_transfer_task.toml"
        )
        cls.paired = load_task_spec(
            EVALUATION / "gse83452_paired_intervention_stress_task.toml"
        )

    def test_gse49541_parent_alias_is_one_family(self) -> None:
        families = {
            value["family_id"]: value for value in self.expansion["cohort_family"]
        }
        family = families["gse31803_gse49541_fibrosis_array"]
        self.assertEqual(set(family["accessions"]), {"GSE31803", "GSE49541"})
        self.assertNotIn("gse49541_fibrosis_array", families)
        self.assertIn(
            "gse31803_gse49541_fibrosis_array",
            self.expansion["activation_order"]["priority_3"],
        )

    def test_scored_preprocessing_is_inductive(self) -> None:
        self.assertFalse(
            self.preprocessing["ordinary_all_sample_RMA_allowed_for_scored_transfer"]
        )
        self.assertFalse(
            self.preprocessing["GEO_series_matrix_allowed_for_scored_transfer"]
        )
        self.assertFalse(
            self.preprocessing["cross_platform_pooling_before_gene_mapping"]
        )
        firewall = self.preprocessing["inductive_firewall"]
        self.assertFalse(firewall["held_cohort_intensity_distribution_visible_to_training"])
        self.assertEqual(
            firewall["held_array_application"],
            "one_array_at_a_time_with_frozen_training_object",
        )
        self.assertIn("all_sample_RMA", firewall["blocked_paths"])
        self.assertIn("probe_effect_or_summarization_parameters", firewall["outer_training_fit_only"])

    def test_antwerp_accessions_are_one_project_exposed_family(self) -> None:
        families = {
            value["family_id"]: value for value in self.expansion["cohort_family"]
        }
        family = families["antwerp_inserm_shared"]
        self.assertEqual(set(family["accessions"]), {"GSE106737", "GSE83452"})
        self.assertEqual(family["cross_accession_fingerprint_matched_arrays"], 78)
        self.assertEqual(family["cross_accession_fingerprint_matched_participants"], 41)
        self.assertEqual(family["exact_family_biological_units"], "join_unresolved")
        self.assertEqual(family["family_union_participant_range"], [171, 200])
        self.assertEqual(family["timepoint_exact_token_pairs"], 38)
        self.assertEqual(family["timepoint_punctuation_only_pairs"], 40)
        self.assertEqual(family["timepoint_semantic_mismatches"], 0)
        self.assertEqual(
            family["cross_accession_overlap_artifacts_sha256"],
            "e577e1f581a484095748a888cfaebcfb83548d452bab523eab8ef19dbb0c89eb",
        )
        self.assertEqual(
            family["runtime_inventory_artifacts_sha256"],
            "2672e60569ef2eb8b50048e019c1897e15d67a7e3c7a3bcd4f89c19a4b04c522",
        )
        self.assertEqual(family["gpl16686_mapping_state"], "unavailable_permission")
        self.assertFalse(family["gpl16686_single_array_summarization_runtime_ready"])
        self.assertNotIn("reported_biological_units", family)
        self.assertNotIn("gse83452_longitudinal_array", families)
        self.assertIn(
            "antwerp_inserm_shared",
            self.expansion["activation_order"]["priority_3"],
        )

    def test_platforms_remain_separate_and_assay_native(self) -> None:
        platforms = self.preprocessing["platform"]
        self.assertEqual(platforms["GPL570"]["CEL_records"], 72)
        self.assertEqual(platforms["GPL16686"]["CEL_records"], 231)
        self.assertEqual(platforms["GPL570"]["CEL_format"], "calvin_v1")
        self.assertEqual(platforms["GPL16686"]["CEL_format"], "calvin_v1")
        self.assertIn("hgu133plus2frmavecs", platforms["GPL570"]["required_runtime_packages"])
        self.assertIn("pd.hugene.2.0.st", platforms["GPL16686"]["required_runtime_packages"])

    def test_tasks_preserve_participant_units_and_claim_limits(self) -> None:
        for task in (self.fibrosis, self.baseline, self.paired):
            self.assertEqual(task.unit_of_inference, "participant")
            self.assertEqual(task.bootstrap_replicates, 10000)
            self.assertIn("champion", task.claim_gate.lower())
        self.assertEqual(self.fibrosis.evaluator_parameters["participants"], 72)
        self.assertFalse(self.fibrosis.evaluator_parameters["age_available"])
        self.assertEqual(
            self.baseline.evaluator_parameters["participants_endpoint_evaluable"],
            148,
        )
        self.assertEqual(self.baseline.evaluator_parameters["undefined"], 4)
        self.assertEqual(
            self.paired.evaluator_parameters["paired_participants_expression"],
            60,
        )
        self.assertEqual(
            self.paired.evaluator_parameters[
                "paired_participants_defined_nash_transition"
            ],
            54,
        )
        self.assertEqual(self.baseline.datasets_development, ("antwerp_inserm_shared",))
        self.assertEqual(self.paired.datasets_train, ("antwerp_inserm_shared",))
        self.assertEqual(self.paired.required_pairing_levels, ())
        self.assertEqual(
            self.paired.evaluator_parameters["longitudinal_topology"],
            "same_participant_repeated_baseline_followup",
        )
        self.assertFalse(
            self.paired.evaluator_parameters["followup_expression_allowed_as_input"]
        )

    def test_source_native_label_semantics_fail_closed(self) -> None:
        self.assertEqual(self.labels["cohort_family_id"], "antwerp_inserm_shared")
        self.assertFalse(
            self.labels["generic_disease_binary_is_compatible_training_label"]
        )
        self.assertEqual(self.labels["target"]["positive_source_value"], "NASH")
        self.assertEqual(self.labels["target"]["negative_source_value"], "no NASH")
        self.assertEqual(self.labels["target"]["missing_source_value"], "undefined")
        self.assertEqual(self.labels["source_training"]["eligible_source_roster"], [])
        self.assertFalse(self.labels["source_training"]["transfer_training_allowed_now"])

    def test_promotion_gate_hashes_match_task_specs(self) -> None:
        pairs = (
            (
                self.fibrosis,
                "gse49541_fibrosis_transfer_promotion_gate.json",
            ),
            (
                self.baseline,
                "gse83452_baseline_nash_transfer_promotion_gate.json",
            ),
            (
                self.paired,
                "gse83452_paired_intervention_stress_promotion_gate.json",
            ),
        )
        for task, name in pairs:
            with self.subTest(task=task.task_id):
                self.assertEqual(
                    task.promotion_gate_config_sha256,
                    sha256_file(EVALUATION / name),
                )


if __name__ == "__main__":
    unittest.main()
