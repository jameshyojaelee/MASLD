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
        self.assertEqual(
            znf["coordinate_activation_artifacts_sha256"],
            "e5042eb43dddd84cda81c5aeed5d6dbeb0378ab1c942a0183e2ef1b41cd18b20",
        )
        self.assertEqual(
            znf["coordinate_audit_sha256"],
            "8d9f01618ff51183a24d70741c8803048552705a166c015c26762ea221523691",
        )
        self.assertEqual(
            znf["coordinate_independent_task_spec_sha256"],
            "8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28",
        )
        self.assertEqual(znf["authoritatively_joined_participants"], 99)
        self.assertEqual(znf["assay_native_qc_passed_participants"], 99)
        self.assertEqual(znf["assay_native_qc_hard_failures"], 0)
        self.assertFalse(znf["new_dua_or_controlled_access_required"])
        self.assertTrue(znf["internal_nonclinical_research_use_allowed"])
        self.assertTrue(znf["open_derivative_weight_review_required"])
        self.assertFalse(znf["model_training_active"])
        self.assertTrue(znf["task_specific_coordinate_independent_training_active"])
        self.assertFalse(znf["coordinate_dependent_training_active"])
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
        self.assertTrue(rrbs["rrbs_bed_qc_complete"])
        self.assertEqual(rrbs["rrbs_files_exhaustively_audited"], 57)
        self.assertEqual(rrbs["rrbs_cytosine_rows_exhaustively_audited"], 800713381)
        self.assertTrue(rrbs["rna_raw_run_manifest_complete"])
        self.assertFalse(rrbs["rna_quantification_complete"])
        self.assertFalse(rrbs["hg19_to_grch38p14_crosswalk_complete"])
        self.assertFalse(rrbs["model_training_active"])
        for key in (
            "activation_execution_artifacts_sha256",
            "activation_fixture_artifacts_sha256",
            "task_spec_sha256",
        ):
            self.assertRegex(rrbs[key], r"^[0-9a-f]{64}$")

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

    def test_gse256398_is_audited_but_full_activation_is_qc_blocked(self) -> None:
        families = {
            family["family_id"]: family
            for family in self.expansion["cohort_family"]
        }
        cohort = families["gse256398"]
        self.assertEqual(
            cohort["status"],
            "compatibility_smoke_registered_full_activation_qc_blocked",
        )
        self.assertEqual(cohort["reported_biological_units"], 26)
        self.assertEqual(cohort["audited_nuclei"], 197942)
        self.assertEqual(cohort["project_reconstructed_source_exact_qc_nuclei"], 165372)
        self.assertEqual(cohort["prospective_qc_sensitivity_nuclei"], 172997)
        self.assertEqual(cohort["compatibility_smoke_nuclei"], 1040)
        self.assertEqual(cohort["compatibility_smoke_nuclei_per_donor"], 40)
        self.assertEqual(cohort["masld_relevant_nuclei"], 135928)
        self.assertEqual(cohort["alcohol_ood_nuclei"], 62014)
        self.assertEqual(cohort["source_features"], 36601)
        self.assertEqual(cohort["gencode49_stable_id_exact_features"], 35455)
        self.assertEqual(cohort["gencode49_masked_features"], 1146)
        self.assertEqual(cohort["age_observed_donors"], 26)
        self.assertEqual(cohort["recorded_sex_observed_donors"], 26)
        self.assertFalse(cohort["authoritative_barcode_cell_labels_available"])
        self.assertTrue(cohort["compatibility_smoke_active"])
        self.assertFalse(cohort["model_training_active"])
        self.assertEqual(cohort["topology"], "same_study_unpaired")

    def test_gse260666_external_fixture_is_project_exposed_and_locked(self) -> None:
        families = {
            family["family_id"]: family
            for family in self.expansion["cohort_family"]
        }
        cohort = families["gse260666_bulk_rna"]
        self.assertEqual(cohort["authoritatively_joined_participants"], 16)
        self.assertEqual(cohort["label_counts"], {"control": 6, "NAFL": 6, "NASH": 4})
        self.assertEqual(cohort["mapped_genes"], 23576)
        self.assertEqual(cohort["genes_shared_with_gse267145"], 20410)
        self.assertEqual(
            cohort["fixture_campaign_artifacts_sha256"],
            "00bed1de86ddf805fde54a9f963fd05b1384aa1d4ee02d427f123a07ee83bc4c",
        )
        self.assertFalse(cohort["project_sealed"])
        self.assertFalse(cohort["model_fitted"])
        self.assertFalse(cohort["outcomes_scored"])
        self.assertIn("selection_lock", cohort["status"])

    def test_gse268273_external_fixture_is_participant_safe_and_raw_blocked(self) -> None:
        families = {
            family["family_id"]: family
            for family in self.expansion["cohort_family"]
        }
        cohort = families["gse268273_imid_masld"]
        self.assertEqual(cohort["authoritatively_joined_participants"], 109)
        self.assertEqual(cohort["technical_runs"], 824)
        self.assertEqual(cohort["effective_single_end_runs"], 680)
        self.assertEqual(cohort["effective_paired_end_runs"], 144)
        self.assertEqual(cohort["effective_single_end_participants"], 73)
        self.assertEqual(cohort["effective_paired_end_participants"], 36)
        self.assertTrue(cohort["declared_single_layout_discrepancy_admitted"])
        self.assertEqual(cohort["source_group_counts"], {"imid_masld": 69, "classic_masld": 40})
        self.assertEqual(sum(cohort["fibrosis_counts"].values()), 109)
        self.assertEqual(sum(cohort["recorded_sex_counts"].values()), 109)
        self.assertEqual(cohort["candidate_unique_gencode_v49_genes"], 14078)
        self.assertEqual(cohort["candidate_aggregated_gencode_v49_genes"], 11)
        self.assertEqual(cohort["raw_fastq_bytes"], 512881099727)
        self.assertFalse(cohort["project_sealed"])
        self.assertFalse(cohort["raw_campaign_submitted"])
        self.assertFalse(cohort["raw_fastq_materialized"])
        self.assertFalse(cohort["model_fitted"])
        self.assertFalse(cohort["outcomes_scored"])
        self.assertFalse(cohort["deposited_voom_scored_input_allowed"])
        self.assertIn("raw_campaign_plan_passed", cohort["status"])
        for key in (
            "fixture_campaign_artifacts_sha256",
            "fixture_artifacts_sha256",
            "source_evidence_artifacts_sha256",
            "preprocessing_input_artifacts_sha256",
            "model_input_artifacts_sha256",
            "evaluator_only_artifacts_sha256",
            "contamination_audit_artifacts_sha256",
            "raw_campaign_plan_execution_artifacts_sha256",
            "raw_campaign_plan_artifacts_sha256",
        ):
            self.assertRegex(cohort[key], r"^[0-9a-f]{64}$")

    def test_gse274114_activation_is_within_instrument_only(self) -> None:
        cohort = next(
            family
            for family in self.expansion["cohort_family"]
            if family["family_id"] == "gse274114_mash_hbv"
        )
        self.assertEqual(cohort["reported_biological_units"], 39)
        self.assertEqual(cohort["technical_runs"], 71)
        self.assertEqual(sum(cohort["group_counts"].values()), 39)
        self.assertEqual(cohort["model_input_genes_gencode_v49"], 60324)
        self.assertEqual(cohort["unmapped_source_rows_masked"], 1230)
        self.assertFalse(cohort["project_sealed"])
        self.assertFalse(cohort["model_fitted"])
        self.assertFalse(cohort["outcomes_scored"])
        self.assertEqual(cohort["tasks"], ["within_instrument_etiology_contrast", "comorbidity_ood_diagnostic"])
        joined_rules = " ".join(cohort["rules"]).lower()
        self.assertIn("perfectly confounded", joined_rules)
        self.assertIn("four-class accuracy", joined_rules)
        self.assertIn("hbv-only participants are not masld controls", joined_rules)
        for key in ("activation_artifacts_sha256", "task_spec_sha256", "promotion_gate_sha256"):
            self.assertRegex(cohort[key], r"^[0-9a-f]{64}$")

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
        histone = self.tournament["task"]["histone_regulatory"]
        self.assertEqual(
            histone["coordinate_independent_task_id"],
            "paired_bulk_rna_h3k27ac",
        )
        self.assertTrue(histone["coordinate_independent_bulk_lane_active"])
        self.assertFalse(histone["coordinate_dependent_sequence_lane_active"])
        self.assertFalse(histone["single_cell_models_allowed_on_bulk_lane"])
        self.assertFalse(histone["source_owned_outcomes_available_for_selection"])
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
