"""The model census must carry a native-task disposition for every family.

OVERALL_PLAN completion criterion 1 requires an exact-checkpoint, license,
exposure, runtime and native-task disposition per family.  Fourteen registry
model_ids across nine families carried no native-task disposition at all, and
ten more carried none for ``sequence_native_regulatory``.  These tests pin the
additive overlay that closes that gap and, critically, pin the fact that it is
additive: the base registry under ``config/models/`` is never consulted as a
mutable surface and every disposition is re-derived from it here.
"""
import hashlib
import json
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "config" / "artifacts" / "models"
NAME = "native_task_disposition_20260825.json"
OVERLAY = "model_census_native_task_disposition_20260825"

GROUP_A = {
    "cellplm": ["cellplm_85m", "cellplm_85m_vae_20231027"],
    "geneformer": [
        "geneformer_v1_10m",
        "geneformer_v2_104m",
        "geneformer_v2_316m",
    ],
    "langcell": ["langcell"],
    "regformer": ["regformer"],
    "scbert": ["scbert"],
    "sccello": ["sccello"],
    "scgpt": ["scgpt_whole_human", "scgpt_continual"],
    "scimilarity": ["scimilarity_v1_1"],
    "uce": ["uce_4l", "uce_33l"],
}

GROUP_B = {
    "alphagenome": "alphagenome",
    "borzoi": "borzoi_ensemble",
    "caduceus": "caduceus",
    "dnabert2": "dnabert2",
    "enformer": "enformer",
    "evo": "evo",
    "hyenadna": "hyenadna",
    "nucleotide_transformer": "nucleotide_transformer",
    "sei": "sei",
    "scooby_onek1k": "scooby_onek1k",
}

GROUP_C = ["epibert", "scooby_epicardioids", "scooby_neurips"]

ALL_FAMILIES = sorted(set(GROUP_A) | set(GROUP_B) | set(GROUP_C))


def load(family: str) -> dict:
    return json.loads((MODELS / family / NAME).read_text(encoding="utf-8"))


def registry_rows() -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for path in sorted((ROOT / "config" / "models").glob("*.toml")):
        with path.open("rb") as handle:
            for model in tomllib.load(handle).get("models", []):
                rows[model["model_id"]] = model
    return rows


class TestNativeTaskDispositionOverlay(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.records = {family: load(family) for family in ALL_FAMILIES}
        cls.rows = registry_rows()

    # -- the overlay exists and is record-only -----------------------------

    def test_every_gap_family_now_carries_a_native_task_disposition(self) -> None:
        self.assertEqual(len(self.records), 22)
        for family, record in self.records.items():
            self.assertEqual(record["family_id"], family)
            self.assertEqual(record["overlay_id"], OVERLAY)
            self.assertEqual(
                record["schema_version"],
                "masld-bench-native-task-disposition-v1",
            )
            self.assertTrue(record["disposition"], family)
            self.assertTrue(record["disposition_reason"], family)

    def test_every_record_is_record_only_and_additive(self) -> None:
        for family, record in self.records.items():
            self.assertTrue(record["record_only"], family)
            self.assertFalse(record["consumed_by_runtime"], family)
            self.assertTrue(record["base_registry_unchanged"], family)
            self.assertTrue(record["supersedes_nothing"], family)

    def test_no_record_touched_an_outcome_metric_or_checkpoint(self) -> None:
        forbidden = (
            "outcomes_read",
            "metrics_calculated",
            "sealed_data_read",
            "biological_data_read",
            "checkpoint_opened",
            "unrestricted_deserializer_used",
            "license_gate_bypassed",
        )
        for family, record in self.records.items():
            for key in forbidden:
                self.assertFalse(record[key], f"{family}:{key}")

    def test_no_record_authorizes_scoring_ranking_or_a_champion(self) -> None:
        for family, record in self.records.items():
            consequences = record["binding_consequences"]
            for key in (
                "sealed_scoring_authorized",
                "sealed_champion_eligible",
                "open_champion_eligible",
                "may_enter_a_sealed_shortlist_or_ranking",
            ):
                self.assertFalse(consequences[key], f"{family}:{key}")

    # -- every binding is live ---------------------------------------------

    def test_every_bound_evidence_digest_matches_the_live_file(self) -> None:
        checked = 0
        for family, record in self.records.items():
            self.assertTrue(record["bound_evidence"], family)
            for item in record["bound_evidence"]:
                path = ROOT / item["path"]
                self.assertTrue(path.is_file(), item["path"])
                actual = hashlib.sha256(path.read_bytes()).hexdigest()
                self.assertEqual(actual, item["sha256"], item["path"])
                checked += 1
        self.assertGreater(checked, 40)

    def test_the_base_registry_is_unchanged_at_the_recorded_digests(self) -> None:
        for family, record in self.records.items():
            self.assertTrue(record["base_registry_files"], family)
            for item in record["base_registry_files"]:
                path = ROOT / item["registry_file"]
                self.assertTrue(path.is_file(), item["registry_file"])
                self.assertEqual(
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    item["registry_file_sha256"],
                    f"{family}:{item['registry_file']}",
                )

    def test_recorded_registry_facts_still_match_the_live_registry(self) -> None:
        """A disposition is only valid for the registry state it was read from."""
        for family, record in self.records.items():
            for model_id, facts in record["registry_facts_at_audit"].items():
                row = self.rows[model_id]
                for key in (
                    "license_status",
                    "exposure_status",
                    "checkpoint_sha256",
                    "status",
                    "admission_blocking",
                ):
                    self.assertEqual(row[key], facts[key], f"{model_id}:{key}")
                self.assertEqual(
                    list(row["supported_tasks"]),
                    facts["supported_tasks"],
                    model_id,
                )

    # -- group A: cell_state_mapping ---------------------------------------

    def test_group_a_covers_all_fourteen_cell_state_model_ids(self) -> None:
        covered: list[str] = []
        for family, expected in GROUP_A.items():
            record = self.records[family]
            self.assertEqual(record["task_id"], "cell_state_mapping")
            self.assertEqual(sorted(record["model_ids"]), sorted(expected))
            covered.extend(record["model_ids"])
        self.assertEqual(len(covered), 14)
        self.assertEqual(len(set(covered)), 14)

    def test_group_a_records_the_missing_capability_roster_as_root_cause(self) -> None:
        capabilities = tomllib.loads(
            (
                ROOT / "config/evaluation/cell_state_mapping_capabilities.toml"
            ).read_text(encoding="utf-8")
        )
        variant = tomllib.loads(
            (
                ROOT / "config/evaluation/variant_to_regulation_capabilities.toml"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(capabilities.get("models", []), [])
        self.assertEqual(len(variant["models"]), 40)
        for family in GROUP_A:
            root_cause = self.records[family]["root_cause"]
            self.assertEqual(root_cause["cell_state_capability_rows"], 0)
            self.assertEqual(root_cause["variant_capability_rows"], 40)

    def test_blocked_is_recorded_with_its_exact_reason(self) -> None:
        expected = {
            "cellplm": "blocked",
            "scbert": "blocked",
            "sccello": "blocked",
            "langcell": "blocked",
            "scgpt": "blocked_pending_weight_terms",
            "regformer": "blocked_exposure_unknown",
        }
        for family, disposition in expected.items():
            record = self.records[family]
            self.assertEqual(record["disposition"], disposition, family)
            self.assertTrue(record["disposition_reason"], family)
            for entry in record["model_dispositions"]:
                self.assertEqual(entry["disposition"], disposition, family)
                self.assertGreater(len(entry["reason"]), 80, entry["model_id"])

    def test_the_four_blocked_families_have_no_declared_weight_terms(self) -> None:
        for family in ("cellplm", "scbert", "langcell", "scgpt"):
            for model_id in GROUP_A[family]:
                self.assertIn(
                    "UNDECLARED", self.rows[model_id]["license_status"], model_id
                )
        self.assertEqual(
            self.rows["sccello"]["license_status"],
            "blocked_no_detected_code_or_weight_license",
        )

    def test_cellplm_records_its_unresolved_checkpoint_release_conflict(self) -> None:
        bundle = json.loads(
            (MODELS / "cellplm" / "checkpoints.json").read_text(encoding="utf-8")
        )
        self.assertIn("checkpoint_release_conflict", bundle)
        self.assertEqual(
            bundle["registered_checkpoint"]["status"],
            "undocumented_release_variant",
        )
        blockers = self.records["cellplm"]["family_specific_blockers"]
        self.assertTrue(blockers["checkpoint_release_conflict_unresolved"])
        self.assertEqual(
            blockers["registered_checkpoint_status"],
            bundle["registered_checkpoint"]["status"],
        )

    def test_scgpt_records_the_missing_runtime_contract_and_gdrive_pin(self) -> None:
        bundle = json.loads(
            (MODELS / "scgpt" / "checkpoints.json").read_text(encoding="utf-8")
        )
        self.assertNotIn("runtime_contract", bundle)
        self.assertFalse(
            self.records["scgpt"]["family_specific_blockers"][
                "runtime_contract_block_present"
            ]
        )
        for model_id in GROUP_A["scgpt"]:
            row = self.rows[model_id]
            self.assertIn("gdrive-folder:", row["checkpoint_revision"])
            self.assertEqual(row["checkpoint_sha256"], "UNRESOLVED")
            self.assertEqual(row["license_status"], "code_MIT_weights_UNDECLARED")

    def test_regformer_is_blocked_on_exposure_not_on_licence(self) -> None:
        row = self.rows["regformer"]
        self.assertEqual(row["exposure_status"], "unknown")
        self.assertNotIn("UNDECLARED", row["license_status"])
        bundle = json.loads(
            (MODELS / "regformer" / "checkpoints.json").read_text(encoding="utf-8")
        )
        self.assertEqual(bundle["canonical_checkpoint"]["sha256"], "UNRESOLVED")
        self.assertEqual(
            self.records["regformer"]["disposition"], "blocked_exposure_unknown"
        )

    # -- group A: the exposure rule and the executed-development split ------

    def test_scimilarity_ineligibility_is_the_frozen_rule_not_a_judgement(self) -> None:
        record = self.records["scimilarity"]
        rule = record["champion_ineligibility_is_a_consequence_of_the_frozen_exposure_rule"]
        self.assertEqual(self.rows["scimilarity_v1_1"]["exposure_status"], "encoder_seen")
        self.assertEqual(rule["model_exposure_state"], "encoder_seen")
        self.assertEqual(
            rule["champion_permitting_states"],
            ["clean_declared", "target_label_unexposed"],
        )
        self.assertNotIn("encoder_seen", rule["champion_permitting_states"])
        self.assertTrue(rule["not_a_quality_judgement"])
        audit = json.loads(
            (MODELS / "scimilarity" / "exposure_audit.json").read_text(encoding="utf-8")
        )
        self.assertFalse(audit["sealed_champion_eligible"])
        self.assertEqual(
            audit["sealed_champion_eligibility"],
            "ineligible_aggregate_development_encoder_seen",
        )

    def test_scimilarity_is_the_only_encoder_seen_row_in_the_census(self) -> None:
        seen = [
            model_id
            for model_id, row in self.rows.items()
            if row["exposure_status"] == "encoder_seen"
        ]
        self.assertEqual(seen, ["scimilarity_v1_1"])

    def test_the_champion_permitting_exposure_subset_is_frozen_upstream(self) -> None:
        record = self.records["scimilarity"]["exposure_rule"]
        self.assertEqual(
            record["champion_permitting_subset"],
            ["clean_declared", "target_label_unexposed"],
        )
        universe = tomllib.loads(
            (ROOT / "config/gse281364_five_seed_row_universe.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            universe["admission_gates"]["allowed_checkpoint_exposure_states"],
            record["champion_permitting_subset"],
        )
        resources = tomllib.loads(
            (ROOT / "config/resources.toml").read_text(encoding="utf-8")
        )
        for state in record["champion_permitting_subset"]:
            self.assertIn(state, resources["semantics"]["exposure_statuses"])
        # The subset is prespecified for the GSE281364 lane; the record must
        # say so rather than assert a rule the cell-state TaskSpec never made.
        self.assertIn("GSE281364", record["scope_caveat"])

    def test_geneformer_exposure_scope_is_stated_not_papered_over(self) -> None:
        scope = self.records["geneformer"]["exposure_scope_resolution"]
        self.assertTrue(scope["resolved"])
        target = scope["sealed_target_the_claim_is_scoped_to"]
        audit = json.loads(
            (MODELS / "geneformer" / "exposure_audit.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            target["sealed_dataset_id"], audit["audit_scope"]["sealed_dataset_id"]
        )
        self.assertEqual(
            target["sealed_endpoint_labels"],
            audit["audit_scope"]["sealed_endpoint_labels"],
        )
        # The registry field and the blockers text describe different axes.
        for token in ("GSE136103", "Liver_Atlas", "DEVELOPMENT", "sealed"):
            self.assertIn(token, scope["statement"], token)
        self.assertTrue(scope["unresolved_remainder"])

    def test_geneformer_and_uce_are_not_called_exposure_ineligible(self) -> None:
        """Only SCimilarity is champion-ineligible on the exposure axis."""
        for family in ("geneformer", "uce"):
            audit = json.loads(
                (MODELS / family / "exposure_audit.json").read_text(encoding="utf-8")
            )
            self.assertTrue(
                audit["sealed_champion_eligibility"].startswith("eligible_on_exposure"),
                family,
            )
            self.assertIn("NOT", self.records[family]["champion_eligibility_axis"])

    def test_only_executed_models_are_called_executed(self) -> None:
        executed = {
            ("geneformer", "geneformer_v2_316m"),
            ("scimilarity", "scimilarity_v1_1"),
        }
        for family in GROUP_A:
            for entry in self.records[family]["model_dispositions"]:
                is_executed = entry["disposition"].startswith("executed_development")
                self.assertEqual(
                    is_executed,
                    (family, entry["model_id"]) in executed,
                    f"{family}:{entry['model_id']}",
                )

    def test_unexecuted_geneformer_and_uce_checkpoints_say_so(self) -> None:
        by_id = {
            entry["model_id"]: entry["disposition"]
            for family in ("geneformer", "uce")
            for entry in self.records[family]["model_dispositions"]
        }
        self.assertEqual(by_id["geneformer_v1_10m"], "acquired_not_executed")
        self.assertEqual(by_id["geneformer_v2_104m"], "acquired_not_executed")
        self.assertEqual(by_id["uce_4l"], "staged_not_executed")
        self.assertEqual(by_id["uce_33l"], "staged_not_executed")

    def test_uce_runtime_axis_is_bound_to_the_existing_frozen_screen(self) -> None:
        bound = self.records["uce"]["runtime_axis_bound_from_existing_artifact"]
        screen = json.loads(
            (MODELS / "uce" / "frozen_screen_50000.json").read_text(encoding="utf-8")
        )
        self.assertEqual(bound["runtime_root"], screen["runtime"]["root"])
        self.assertEqual(
            bound["runtime_artifacts_sha256"], screen["runtime"]["artifacts_sha256"]
        )
        self.assertEqual(bound["python"], screen["runtime"]["python"])
        self.assertEqual(bound["torch"], screen["runtime"]["torch"])
        self.assertEqual(bound["cuda_build"], screen["runtime"]["cuda_build"])
        self.assertEqual(bound["split_id"], screen["split_id"])
        self.assertEqual(bound["dataset_view_id"], screen["dataset_view_id"])
        self.assertEqual(bound["rows"], screen["rows"])
        self.assertIn("not an executed forward pass", bound["statement"])

    def test_uce_registry_artifact_checkpoint_disagreement_is_recorded(self) -> None:
        bundle = json.loads(
            (MODELS / "uce" / "checkpoints.json").read_text(encoding="utf-8")
        )
        self.assertTrue(bundle["weight_content_downloaded"])
        for model_id in ("uce_4l", "uce_33l"):
            self.assertEqual(self.rows[model_id]["checkpoint_sha256"], "UNRESOLVED")
        recorded = self.records["uce"]["recorded_inconsistencies"]
        observed = " ".join(item["observed"] for item in recorded)
        for model_id in ("uce_4l", "uce_33l"):
            self.assertIn(bundle["weight_artifacts"][model_id]["sha256"], observed)

    def test_geneformer_stale_download_flag_is_recorded(self) -> None:
        """The recorded defect stays on record after it is corrected.

        This overlay recorded weight_content_downloaded=false as stale while
        all three bodies sat on disk, and deferred the fix to the owning lane.
        That lane corrected the bundle to true on 2026-08-25.  The assertion
        tracks the resolution rather than the old value, but still pins that
        the record preserves what was originally wrong -- a debt record that
        forgets its paid items loses the history.
        """
        bundle = json.loads(
            (MODELS / "geneformer" / "checkpoints.json").read_text(encoding="utf-8")
        )
        recorded = self.records["geneformer"]["recorded_inconsistencies"]
        self.assertEqual(len(recorded), 1)
        entry = recorded[0]
        self.assertIn("weight_content_downloaded", entry["field"])
        self.assertFalse(entry["recorded_value"], "the original defect must stay on record")
        if bundle["weight_content_downloaded"]:
            self.assertEqual(entry["action"], "corrected_2026_08_25")
            self.assertIn("weight_content_downloaded set to true", entry["correction"])
        else:
            self.assertEqual(
                entry["action"], "recorded_only_not_corrected_by_this_overlay"
            )

    # -- group B: sequence_native_regulatory --------------------------------

    def test_group_b_covers_all_ten_sequence_model_ids(self) -> None:
        tournament = tomllib.loads(
            (ROOT / "config/evaluation/family_native_tournament.toml").read_text(
                encoding="utf-8"
            )
        )
        required = sorted(
            binding["model_id"]
            for binding in tournament["required_family_disposition"]
            if "sequence_native_regulatory" in binding["required_tasks"]
        )
        self.assertEqual(len(required), 10)
        self.assertEqual(sorted(GROUP_B.values()), required)
        for family, model_id in GROUP_B.items():
            record = self.records[family]
            self.assertEqual(record["task_id"], "sequence_native_regulatory")
            self.assertEqual(record["model_ids"], [model_id])
            self.assertEqual(record["disposition"], "blocked_task_not_registered")
            self.assertEqual(
                record["model_dispositions"][0]["tournament_disposition_rule"],
                next(
                    b["disposition_rule"]
                    for b in tournament["required_family_disposition"]
                    if b["model_id"] == model_id
                ),
            )

    def test_the_sequence_task_is_genuinely_unregistered(self) -> None:
        tournament = tomllib.loads(
            (ROOT / "config/evaluation/family_native_tournament.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            tournament["task"]["sequence_native_regulatory"]["status"],
            "planned_registration",
        )
        self.assertIn(
            "sequence_native_regulatory", tournament["task_ids_to_register"]
        )
        self.assertNotIn(
            "sequence_native_regulatory", tournament["active_task_ids"]
        )
        registered = sorted(
            tomllib.loads(path.read_text(encoding="utf-8"))["task_id"]
            for path in (ROOT / "config" / "tasks").glob("*.toml")
        )
        self.assertNotIn("sequence_native_regulatory", registered)
        for family in GROUP_B:
            cause = self.records[family]["root_cause"]["reason_1_task_not_registered"]
            self.assertFalse(cause["task_spec_exists"])
            self.assertEqual(cause["registered_task_ids"], registered)

    def test_the_sequence_artifacts_are_bound_to_rna_conditioned_atac(self) -> None:
        """The mislabelling that cost this campaign time must stay on record."""
        evaluation = ROOT / "config" / "evaluation"
        declared = {
            "config/evaluation/sequence_terminal_admission.json": json.loads(
                (evaluation / "sequence_terminal_admission.json").read_text(
                    encoding="utf-8"
                )
            ).get("task_id"),
            "config/evaluation/sequence_task_native_five_seed_rectangular_evaluator_20260825.json": json.loads(
                (
                    evaluation
                    / "sequence_task_native_five_seed_rectangular_evaluator_20260825.json"
                ).read_text(encoding="utf-8")
            ).get("task_id"),
            "config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json": json.loads(
                (
                    evaluation
                    / "sequence_regulatory_full_rectangle_completion_20260825.json"
                ).read_text(encoding="utf-8")
            ).get("task_id"),
        }
        self.assertEqual(
            declared[
                "config/evaluation/sequence_terminal_admission.json"
            ],
            "rna_conditioned_atac",
        )
        self.assertEqual(
            declared[
                "config/evaluation/sequence_task_native_five_seed_rectangular_evaluator_20260825.json"
            ],
            "rna_conditioned_atac",
        )
        self.assertIsNone(
            declared[
                "config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json"
            ]
        )
        for family in GROUP_B:
            cause = self.records[family]["root_cause"][
                "reason_2_sequence_artifacts_are_bound_to_a_different_task"
            ]
            for artifact in cause["artifacts"]:
                self.assertEqual(
                    artifact["declared_task_id"], declared[artifact["path"]]
                )

    def test_no_bundle_anywhere_declares_a_sequence_native_task_fit(self) -> None:
        for path in sorted(MODELS.glob("*/checkpoints.json")):
            fit = json.loads(path.read_text(encoding="utf-8")).get("task_fit")
            if isinstance(fit, dict):
                self.assertNotIn("sequence_native_regulatory", fit, path.parent.name)

    def test_the_sequence_lane_is_recorded_as_l40s_bound(self) -> None:
        bpnet = json.loads(
            (MODELS / "bpnet" / "checkpoints.json").read_text(encoding="utf-8")
        )["runtime_contract"]
        chrombpnet = json.loads(
            (MODELS / "chrombpnet" / "checkpoints.json").read_text(encoding="utf-8")
        )["runtime_contract"]
        self.assertEqual(bpnet["declared_tensorflow"], "2.4.1")
        self.assertEqual(bpnet["compatibility_runtime"], "tensorflow_2.8.2_gpu")
        self.assertEqual(bpnet["device_class"], "NVIDIA_L40S_sm89")
        self.assertEqual(
            bpnet["environment_status"],
            "PASSED_COMPATIBILITY_NOT_NATIVE_VERSION_PARITY",
        )
        self.assertEqual(chrombpnet["declared_tensorflow"], "2.8.0")
        for family in GROUP_B:
            constraint = self.records[family]["sequence_lane_runtime_constraint"]
            self.assertEqual(
                constraint["constraint"],
                "the_task_native_sequence_rectangle_is_L40S_bound",
            )
            self.assertEqual(
                constraint["bpnet_declared_tensorflow"],
                bpnet["declared_tensorflow"],
            )
            self.assertEqual(
                constraint["bpnet_device_class"], bpnet["device_class"]
            )
            self.assertEqual(
                constraint["bpnet_environment_status"], bpnet["environment_status"]
            )
            self.assertEqual(
                constraint["chrombpnet_declared_tensorflow"],
                chrombpnet["declared_tensorflow"],
            )
            self.assertIn("sm_120", constraint["statement"])

    def test_the_b6k_comparison_receipt_agrees_with_the_recorded_constraint(self) -> None:
        summary = json.loads(
            (
                ROOT
                / "executions/b6k-admission-comparison-21106761/comparison_summary.json"
            ).read_text(encoding="utf-8")
        )
        self.assertFalse(summary["metrics_calculated"])
        self.assertFalse(summary["outcomes_read"])
        self.assertTrue(summary["synthetic_inputs_only"])
        for result in summary["results"].values():
            self.assertEqual(result["verdict"], "materially_different")
        constraint = self.records["borzoi"]["sequence_lane_runtime_constraint"]
        self.assertTrue(constraint["b6k_admissible_only_for_rectangles_fit_entirely_on_b6k"])
        for token in ("21106572", "21106573", "no_matching_sass"):
            self.assertIn(token, constraint["cluster_fact"], token)

    # -- group C: observed_multiome ----------------------------------------

    def test_group_c_covers_the_three_models_with_no_observed_multiome_fit(self) -> None:
        for family in GROUP_C:
            fit = json.loads(
                (MODELS / family / "checkpoints.json").read_text(encoding="utf-8")
            )["task_fit"]
            self.assertNotIn("observed_multiome", fit, family)
            record = self.records[family]
            self.assertEqual(record["task_id"], "observed_multiome")
            self.assertEqual(
                record["disposition"],
                "registered_candidate_execution_not_authorized",
            )

    def test_group_c_disposition_matches_the_capability_registry(self) -> None:
        capabilities = tomllib.loads(
            (
                ROOT / "config/evaluation/observed_multiome_task_capabilities.toml"
            ).read_text(encoding="utf-8")
        )
        self.assertFalse(capabilities["execution_authorized"])
        self.assertEqual(capabilities["profile_fixture_passed_models"], [])
        self.assertEqual(capabilities["primary_eligible_models"], [])
        registered = set(capabilities["observed_atac_conditioned_track_candidates"]) | set(
            capabilities["sequence_plus_rna_track_candidates"]
        )
        for family in GROUP_C:
            record = self.records[family]
            self.assertIn(record["model_ids"][0], registered, family)
            self.assertIn(
                record["model_dispositions"][0]["capability_track"], capabilities
            )
            self.assertFalse(record["binding_consequences"]["execution_authorized"])
            self.assertFalse(
                record["binding_consequences"]["rna_conditioned_atac_claim_allowed"]
            )
            self.assertEqual(
                record["required_future_baseline_ids"],
                capabilities["required_future_baseline_ids"],
            )

    def test_group_c_disposition_matches_the_additive_census(self) -> None:
        base = ROOT / "executions/model-check-288-21098630"
        receipt = json.loads((base / "receipt.json").read_text(encoding="utf-8"))
        self.assertFalse(receipt["execution_authorized"])
        self.assertFalse(receipt["candidate_fixture_requirement_satisfied"])
        self.assertFalse(receipt["external_seal_bound"])
        rows = {}
        lines = (base / "candidate_models.tsv").read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        for line in lines[1:]:
            if line.strip():
                row = dict(zip(header, line.split("\t")))
                rows[row["model_id"]] = row
        for family in GROUP_C:
            entry = self.records[family]["model_dispositions"][0]
            row = rows[entry["model_id"]]
            self.assertEqual(entry["census_role"], row["role"])
            self.assertEqual(entry["census_output_family"], row["output_family"])
            self.assertEqual(
                entry["census_query_modalities"],
                row["query_modalities"].split(";"),
            )
            census = self.records[family]["root_cause"]["census_facts"]
            self.assertEqual(
                census["candidate_model_count"], receipt["candidate_model_count"]
            )
            self.assertFalse(census["execution_authorized"])


if __name__ == "__main__":
    unittest.main()
