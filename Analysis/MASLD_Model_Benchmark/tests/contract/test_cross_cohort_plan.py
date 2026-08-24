"""Fail-closed checks for the plan-only cross-cohort tournament expansion."""

from __future__ import annotations

from pathlib import Path
import tomllib
import unittest

from masld_bench.registry import load_dataset_manifest


ROOT = Path(__file__).resolve().parents[2]
EVALUATION = ROOT / "config" / "evaluation"


def _load(name: str) -> dict:
    with (EVALUATION / name).open("rb") as handle:
        return tomllib.load(handle)


class CrossCohortPlanContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.expansion = _load("cross_cohort_expansion.toml")
        cls.tournament = _load("family_native_tournament.toml")

    def test_expansion_is_plan_only_and_gse296875_is_not_universal(self) -> None:
        self.assertEqual(
            self.expansion["status"], "plan_only_until_dataset_manifest_activation"
        )
        self.assertFalse(self.expansion["candidate_record_is_active_dataset"])
        self.assertFalse(self.expansion["gse296875_is_sole_training_source"])
        self.assertFalse(self.expansion["gse296875_is_sole_selection_source"])
        self.assertFalse(self.expansion["gse296875_is_masld_phenotype_source"])

    def test_cross_study_and_missingness_guards_are_locked(self) -> None:
        contract = self.expansion["generalization_contract"]
        forbidden = set(contract["forbidden"])
        self.assertEqual(
            contract["three_or_more_compatible_families"],
            "leave_one_study_out_outer_evaluation",
        )
        self.assertIn("random_cell_outer_split", forbidden)
        self.assertIn("false_cell_pairing", forbidden)
        self.assertIn("missing_as_zero", forbidden)
        self.assertTrue(self.expansion["same_person_all_modalities_same_outer_fold"])
        self.assertFalse(self.expansion["missing_modality_encoded_as_zero"])

    def test_family_ids_and_source_accessions_are_unique(self) -> None:
        families = self.expansion["cohort_family"]
        family_ids = [family["family_id"] for family in families]
        self.assertEqual(len(family_ids), len(set(family_ids)))
        znf = next(
            family
            for family in families
            if family["family_id"] == "gse267145_znf469_human_liver"
        )
        self.assertEqual(
            set(znf["accessions"]), {"GSE267145", "GSE267119", "GSE269412"}
        )
        self.assertEqual(znf["reported_participants_in_source_paper"], 108)
        self.assertEqual(znf["reported_h3k27ac_sample_records"], 99)
        self.assertEqual(znf["candidate_title_join_exact_one_to_one"], 99)
        self.assertEqual(znf["candidate_title_join_h3k27ac_only"], 0)
        self.assertTrue(znf["candidate_title_join_authoritative"])
        self.assertRegex(znf["metadata_audit_artifacts_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(znf["matrix_audit_artifacts_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(
            znf["authoritative_join_artifacts_sha256"], r"^[0-9a-f]{64}$"
        )
        self.assertRegex(znf["rna_measurement_artifacts_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(
            znf["reference_crosswalk_artifacts_sha256"], r"^[0-9a-f]{64}$"
        )
        self.assertEqual(
            znf["qc_rights_artifacts_sha256"],
            "ac7c92fd01b2b4240b380e4f6a58d98f1a6f853f8d0ad98d7167d7af34baaa73",
        )
        self.assertEqual(znf["authoritatively_joined_participants"], 99)
        self.assertEqual(znf["assay_native_qc_passed_participants"], 99)
        self.assertEqual(znf["assay_native_qc_hard_failures"], 0)
        self.assertFalse(znf["new_dua_or_controlled_access_required"])
        self.assertTrue(znf["internal_nonclinical_research_use_allowed"])
        self.assertTrue(znf["open_derivative_weight_review_required"])
        self.assertFalse(znf["model_training_active"])
        self.assertTrue(znf["continuous_rna_input_allowed"])
        self.assertFalse(znf["raw_integer_rna_count_likelihood_allowed"])
        self.assertEqual(znf["rna_genes_exact_v98_v49"], 42163)
        self.assertEqual(znf["rna_genes_masked_retired_or_absent_v49"], 1122)
        self.assertEqual(znf["h3_regions_on_grch38p14_primary_contigs"], 96460)
        self.assertFalse(znf["h3_coordinate_base_resolved"])
        self.assertFalse(znf["h3_sequence_extraction_allowed"])
        self.assertEqual(znf["matrix_paired_candidate_ids"], 99)
        self.assertEqual(znf["rna_matrix_feature_rows"], 43285)
        self.assertEqual(znf["h3k27ac_matrix_feature_rows"], 96460)
        self.assertFalse(znf["rna_matrix_nonnegative_integer"])
        self.assertTrue(znf["h3k27ac_matrix_nonnegative_integer"])
        self.assertIn("topology_activation_passed", znf["status"])
        rrbs = next(
            family
            for family in families
            if family["family_id"] == "gse105127_zonated_rna_rrbs"
        )
        self.assertEqual(
            rrbs["authoritative_join_rights_artifacts_sha256"],
            "574eeee1119054a2d7b55066804accf29fb45f23b591dd26f9d61d907b27a75f",
        )
        self.assertEqual(rrbs["authoritatively_joined_participants"], 19)
        self.assertEqual(rrbs["authoritatively_joined_assay_records"], 114)
        self.assertEqual(rrbs["participant_safe_outer_folds"], 5)
        self.assertTrue(rrbs["candidate_title_join_authoritative"])
        self.assertEqual(rrbs["topology"], "adjacent_section")
        self.assertFalse(rrbs["rrbs_bed_qc_complete"])
        self.assertFalse(rrbs["hg19_to_grch38p14_crosswalk_complete"])
        self.assertFalse(rrbs["model_training_active"])

    def test_existing_atlas_sources_cannot_be_double_counted(self) -> None:
        atlas = next(
            family
            for family in self.expansion["cohort_family"]
            if family["family_id"] == "resource_atlas_current"
        )
        self.assertIn("GSE202379", atlas["accessions"])
        self.assertIn("GSE244832", atlas["accessions"])
        self.assertTrue(self.expansion["cohort_family_dedup_precedes_split"])
        self.assertFalse(
            self.tournament["dataset_portfolio"][
                "same_accession_may_count_as_independent_cohort_twice"
            ]
        )

    def test_independent_atac_sources_require_raw_fragment_rebuild(self) -> None:
        families = {
            family["family_id"]: family
            for family in self.expansion["cohort_family"]
        }
        expected_sha = (
            "ef43d3c3908c5c81554d750a928adf42f4bec0c7a90f519cdaf062e145e89392"
        )
        gse244 = families["gse244832"]
        gse281 = families["gse281367"]
        self.assertEqual(gse244["source_audit_artifacts_sha256"], expected_sha)
        self.assertEqual(gse281["source_audit_artifacts_sha256"], expected_sha)
        self.assertEqual(gse244["audited_atac_cells"], 88814)
        self.assertEqual(gse281["audited_atac_cells"], 226224)
        self.assertFalse(gse244["legacy_pseudobulk_usable"])
        self.assertFalse(gse281["legacy_pseudobulk_usable"])
        outcome_sha = (
            "6f5092f1343b91ff592851d7d265ab45f8e154d210984c19f96bf05421025be3"
        )
        membership_sha = (
            "32e1d7b3ba147aea6493a0fcdabe862c1eb887a9a8ec84a1d972b5adbc5b9db6"
        )
        for family in (gse244, gse281):
            self.assertEqual(family["transport_outcomes_artifacts_sha256"], outcome_sha)
            self.assertEqual(family["membership_artifacts_sha256"], membership_sha)
            self.assertTrue(family["development_transport_active"])
            self.assertFalse(family["model_training_active"])

    def test_new_modalities_have_assay_native_planned_tasks(self) -> None:
        planned = set(self.tournament["task_ids_to_register"])
        self.assertTrue(
            {"histone_regulatory", "methylation_regulatory", "protein_transport"}
            <= planned
        )
        self.assertFalse(
            self.tournament["task"]["histone_regulatory"][
                "single_cohort_external_claim_allowed"
            ]
        )
        self.assertFalse(
            self.tournament["task"]["methylation_regulatory"][
                "single_cohort_external_claim_allowed"
            ]
        )
        self.assertTrue(
            self.tournament["dataset_portfolio"]["assay_native_heads_required"]
        )

    def test_high_priority_sources_are_blocked_dataset_manifests(self) -> None:
        datasets = ROOT / "config" / "datasets"
        znf = load_dataset_manifest(
            datasets / "gse267145_znf469_human_liver.toml"
        )
        rrbs = load_dataset_manifest(datasets / "gse105127_zonated_rna_rrbs.toml")
        self.assertEqual(znf.role.value, "blocked")
        self.assertEqual(znf.status.value, "blocked")
        self.assertEqual(znf.expected_biological_units, 108)
        self.assertEqual(znf.modality_status_default.value, "observed")
        self.assertEqual(rrbs.role.value, "blocked")
        self.assertEqual(rrbs.status.value, "blocked")
        self.assertEqual(rrbs.expected_biological_units, 19)
        self.assertEqual(rrbs.modality_status_default.value, "observed")
        self.assertEqual(rrbs.pairing_levels[0].value, "adjacent_section")
        self.assertEqual(rrbs.split.label_visibility, "development_visible")
        self.assertIn("rrbs_methylation", rrbs.modalities)

    def test_required_model_roster_remains_complete(self) -> None:
        observed = {
            record["model_id"]
            for record in self.tournament["required_family_disposition"]
        }
        required = {
            "alphagenome",
            "borzoi_ensemble",
            "caduceus",
            "context_borzoi",
            "corgi_plus",
            "corgi_regular",
            "dnabert2",
            "enformer",
            "epcotv2",
            "epibert",
            "evo",
            "hyenadna",
            "nucleotide_transformer",
            "scooby_epicardioids",
            "scooby_neurips",
            "scooby_onek1k",
            "sei",
        }
        self.assertEqual(observed, required)

    def test_human_plan_preserves_goal_and_cross_cohort_requirements(self) -> None:
        plan = " ".join(
            (ROOT / "OVERALL_PLAN.md").read_text(encoding="utf-8").split()
        )
        required_phrases = {
            "The benchmark is the selection and proof mechanism, not the scientific end.",
            "No champion may be selected from GSE296875 alone",
            "leave-one-study-out outer evaluation",
            "Observed-ATAC models are not excluded.",
            "Conditional context model",
            "one-time sealed evaluation",
            "Borzoi",
            "Corgi+",
            "AlphaGenome",
            "EpiBERT",
            "EPCOTv2",
            "scooby",
        }
        for phrase in required_phrases:
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, plan)


if __name__ == "__main__":
    unittest.main()
