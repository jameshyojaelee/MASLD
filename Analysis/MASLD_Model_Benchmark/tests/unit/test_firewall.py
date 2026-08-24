from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import canonical_hash, write_json_exclusive
from masld_bench.contracts import DatasetManifest
from masld_bench.firewall import (
    CampaignPhase,
    DEFAULT_RESOURCE_AUTHORITIES,
    FirewallError,
    SealedArtifact,
    UNIT_TEST_RESAMPLING_POLICY,
    assert_dataset_access,
    authorize_one_time_outcome_join,
    commit_predictions,
    freeze_development_power_evidence,
    freeze_power_decision,
    assert_resource_unchanged,
    resource_snapshot,
    verify_development_power_evidence,
    verify_outcome_consumption,
    verify_power_decision,
    verify_prediction_commit,
    _validate_resampling_policy,
)

from tests.unit._sealed_fixtures import (
    TASK_IDS,
    commit_prediction_set,
    create_joint_outcome_bundle,
    create_selection,
    freeze_power_set,
)


PACKAGE_ROOT = Path(__file__).resolve().parents[2]


class DatasetAccessTests(unittest.TestCase):
    def test_access_is_derived_from_a_strict_manifest(self) -> None:
        withheld = DatasetManifest.load_toml(
            PACKAGE_ROOT / "config" / "datasets" / "gse289173.toml"
        )
        with self.assertRaises(FirewallError):
            assert_dataset_access(
                dataset=withheld,
                phase=CampaignPhase.SEALED_INFERENCE,
                artifact=SealedArtifact.LABELS,
            )
        assert_dataset_access(
            dataset=withheld.to_dict(),
            phase=CampaignPhase.SEALED_EVALUATION,
            artifact=SealedArtifact.LABELS,
        )

        available = DatasetManifest.load_toml(
            PACKAGE_ROOT / "config" / "datasets" / "gse256398.toml"
        )
        assert_dataset_access(
            dataset=available,
            phase=CampaignPhase.DEVELOPMENT,
            artifact=SealedArtifact.LABELS,
        )

    def test_blocked_or_deferred_manifest_is_never_accessible(self) -> None:
        deferred = DatasetManifest.load_toml(
            PACKAGE_ROOT / "config" / "datasets" / "steatosite.toml"
        )
        with self.assertRaisesRegex(FirewallError, "deferred"):
            assert_dataset_access(
                dataset=deferred,
                phase=CampaignPhase.SEALED_EVALUATION,
                artifact=SealedArtifact.FEATURES,
            )

    def test_role_string_cannot_impersonate_a_manifest(self) -> None:
        with self.assertRaises(TypeError):
            assert_dataset_access(  # type: ignore[call-arg]
                dataset_role="project_sealed_final",
                phase=CampaignPhase.SEALED_EVALUATION,
                artifact=SealedArtifact.LABELS,
            )


class ResourceAuthorityFirewallTests(unittest.TestCase):
    def test_default_authorities_cover_exact_resource_surfaces(self) -> None:
        authorities = set(DEFAULT_RESOURCE_AUTHORITIES)
        self.assertIn("figures/main", authorities)
        self.assertIn("masld-atlas-v2/src", authorities)
        self.assertIn("masld-atlas-v2/public/data", authorities)
        self.assertIn(
            "Analysis/Multimodal_Program_Projection/candidates/"
            "gene-catalog-v2-contract-2026-08-11",
            authorities,
        )
        self.assertNotIn("figures/figures_manifest.tsv", authorities)
        self.assertNotIn("scripts/portal", authorities)

    def test_directory_membership_addition_and_deletion_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            authority = root / "authority"
            authority.mkdir()
            original = authority / "original.txt"
            original.write_text("original\n", encoding="utf-8")
            baseline = resource_snapshot(root, ("authority",))
            observed = assert_resource_unchanged(baseline, repository_root=root)
            self.assertEqual(
                observed["snapshot_sha256"], baseline["snapshot_sha256"]
            )

            added = authority / "added.txt"
            added.write_text("added\n", encoding="utf-8")
            with self.assertRaisesRegex(FirewallError, "changed"):
                assert_resource_unchanged(baseline, repository_root=root)
            added.unlink()
            original.unlink()
            with self.assertRaisesRegex(FirewallError, "changed"):
                assert_resource_unchanged(baseline, repository_root=root)


class SealedEvaluationTests(unittest.TestCase):
    def _ready_state(
        self, root: Path
    ) -> tuple[Path, Path, dict[str, object], Path]:
        lock_dir = create_selection(root)
        state = root / "evaluator"
        commits = commit_prediction_set(
            root,
            selection_lock_dir=lock_dir,
            evaluator_state_dir=state,
        )
        powers = freeze_power_set(evaluator_state_dir=state, commits=commits)
        outcome = create_joint_outcome_bundle(root)
        return lock_dir, state, powers, outcome

    def test_joint_bundle_is_consumed_once_after_all_tasks_are_ready(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, outcome = self._ready_state(root)
            ledger = root / "independent-ledger"
            ledger.mkdir()
            consumption = authorize_one_time_outcome_join(
                evaluator_state_dir=state,
                power_decision_sha256s=(
                    power.power_decision_sha256 for power in powers.values()
                ),
                sealed_outcome_bundle_dir=outcome,
                consumption_ledger_dir=ledger,
            )
            self.assertTrue(Path(consumption.marker_path).is_file())
            self.assertEqual(
                verify_outcome_consumption(consumption.marker_path).to_dict(),
                consumption.to_dict(),
            )
            self.assertEqual(
                {binding["task_id"] for binding in consumption.task_bindings},
                set(TASK_IDS),
            )
            variant_power = powers["variant_to_regulation"]
            self.assertEqual(
                sum(
                    binding["model_role"] == "secondary_comparator"
                    for binding in variant_power.prediction_bindings
                ),
                3,
            )
            with self.assertRaisesRegex(FirewallError, "already consumed"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=(
                        power.power_decision_sha256 for power in powers.values()
                    ),
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=ledger,
                )

    def test_joint_bundle_cannot_open_one_task_early(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, outcome = self._ready_state(root)
            ledger = root / "independent-ledger"
            ledger.mkdir()
            one_power = next(iter(powers.values()))
            with self.assertRaisesRegex(FirewallError, "every task"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=(one_power.power_decision_sha256,),
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=ledger,
                )
            self.assertEqual(list(ledger.iterdir()), [])

    def test_joint_bundle_must_use_the_locked_prediction_unit_identifier(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, _ = self._ready_state(root)
            outcome = create_joint_outcome_bundle(
                root / "wrong-unit", unit_id_field="cell_id"
            )
            ledger = root / "independent-ledger"
            ledger.mkdir()
            with self.assertRaisesRegex(FirewallError, "unit identifiers"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=(
                        power.power_decision_sha256 for power in powers.values()
                    ),
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=ledger,
                )
            self.assertEqual(list(ledger.iterdir()), [])

    def test_outcome_row_to_unit_assignments_must_exactly_match_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, _ = self._ready_state(root)
            outcome = create_joint_outcome_bundle(
                root / "mismatched-outcome",
                unit_offsets={TASK_IDS[0]: 1},
            )
            ledger = root / "independent-ledger"
            ledger.mkdir()
            with self.assertRaisesRegex(FirewallError, "unit or row set differs"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=(
                        power.power_decision_sha256 for power in powers.values()
                    ),
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=ledger,
                )
            self.assertEqual(len(list(ledger.iterdir())), 1)

    def test_valid_json_file_cannot_impersonate_frozen_selection_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_lock = root / "selection_lock.json"
            bundle = root / "prediction_bundle.json"
            write_json_exclusive(fake_lock, {"locked": True})
            write_json_exclusive(bundle, {"schema_version": "forged"})
            with self.assertRaisesRegex(FirewallError, "frozen SelectionLock"):
                commit_predictions(
                    selection_lock_dir=fake_lock,
                    prediction_bundle_path=bundle,
                    evaluator_state_dir=root / "evaluator",
                )

    def test_prediction_artifact_mutation_blocks_power_freeze(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            task_commits = commits[TASK_IDS[0]]
            artifact = (
                Path(task_commits[0].prediction_bundle_path).parent
                / task_commits[0].standardized_table.path
            )
            artifact.write_text('{"forged":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(FirewallError, "reverification"):
                freeze_power_decision(
                    evaluator_state_dir=state,
                    prediction_commit_sha256s=(
                        commit.prediction_commit_sha256 for commit in task_commits
                    ),
                    development_power_evidence_dir=root / "unused-evidence",
                )

    def test_fixture_power_evidence_cannot_authorize_sealed_predictions(self) -> None:
        self.assertEqual(
            _validate_resampling_policy(
                UNIT_TEST_RESAMPLING_POLICY, 100, allow_fixture_policy=True
            ),
            (UNIT_TEST_RESAMPLING_POLICY, 100),
        )
        with self.assertRaisesRegex(FirewallError, "explicit fixture-policy opt-in"):
            _validate_resampling_policy(
                UNIT_TEST_RESAMPLING_POLICY, 100, allow_fixture_policy=False
            )

    def test_power_evidence_rehashes_its_endpoint_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            finalists = sorted(
                (root / "finalist_metric_bundles" / TASK_IDS[0]).iterdir()
            )
            self.assertEqual(len(finalists), 1)
            evidence_dir = freeze_development_power_evidence(
                finalist_development_metric_bundle_dir=finalists[0],
                output_root=state,
            )
            distribution = evidence_dir / "bootstrap_distribution.json"
            distribution.chmod(0o640)
            distribution.write_text(
                distribution.read_text(encoding="utf-8") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(FirewallError, "(?:size|checksum) mismatch"):
                verify_development_power_evidence(evidence_dir)

    def test_prediction_bundle_cannot_switch_its_locked_dataset(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            valid_bundle = Path(commits[TASK_IDS[0]][0].prediction_bundle_path)
            payload = json.loads(valid_bundle.read_text(encoding="utf-8"))
            payload["dataset_ids"] = ["invented_seal"]
            payload["source_join_key_sha256"] = canonical_hash(
                {
                    "task_id": payload["task_id"],
                    "dataset_ids": payload["dataset_ids"],
                    "split_id": payload["split_id"],
                    "row_id_field": payload["row_id_field"],
                    "unit_id_field": payload["unit_id_field"],
                    "unit_id_namespace": payload["unit_id_namespace"],
                    "biological_unit": payload["biological_unit"],
                }
            )
            forged_bundle = valid_bundle.parent / "forged_dataset_bundle.json"
            write_json_exclusive(forged_bundle, payload)
            with self.assertRaisesRegex(FirewallError, "sealed binding"):
                commit_predictions(
                    selection_lock_dir=lock_dir,
                    prediction_bundle_path=forged_bundle,
                    evaluator_state_dir=state,
                )

    def test_power_decision_requires_every_locked_seed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            task_commits = commits[TASK_IDS[0]][:-1]
            with self.assertRaisesRegex(FirewallError, "selected and baseline run"):
                freeze_power_decision(
                    evaluator_state_dir=state,
                    prediction_commit_sha256s=(
                        commit.prediction_commit_sha256 for commit in task_commits
                    ),
                    development_power_evidence_dir=root / "unused-evidence",
                )

    def test_variant_power_decision_requires_every_secondary_comparator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            task_id = "variant_to_regulation"
            secondary = [
                commit
                for commit in commits[task_id]
                if commit.model_role == "secondary_comparator"
            ]
            self.assertEqual(len(secondary), 3)
            incomplete = tuple(
                commit
                for commit in commits[task_id]
                if commit.prediction_commit_sha256
                != secondary[-1].prediction_commit_sha256
            )
            with self.assertRaisesRegex(FirewallError, "secondary comparator run"):
                freeze_power_set(
                    evaluator_state_dir=state,
                    commits={task_id: incomplete},
                )

    def test_secondary_comparator_bundle_cannot_claim_primary_endpoint_scope(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            comparator = next(
                commit
                for commit in commits["variant_to_regulation"]
                if commit.model_role == "secondary_comparator"
            )
            bundle_path = Path(comparator.prediction_bundle_path)
            payload = json.loads(bundle_path.read_text(encoding="utf-8"))
            payload["metadata"]["primary_endpoint_scoring_allowed"] = True
            forged_bundle = bundle_path.parent / "forged_primary_scope_bundle.json"
            write_json_exclusive(forged_bundle, payload)
            with self.assertRaisesRegex(FirewallError, "endpoint metadata"):
                commit_predictions(
                    selection_lock_dir=lock_dir,
                    prediction_bundle_path=forged_bundle,
                    evaluator_state_dir=state,
                )

    def test_unlocked_secondary_run_model_cannot_be_committed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            comparator = next(
                commit
                for commit in commits["variant_to_regulation"]
                if commit.model_role == "secondary_comparator"
            )
            bundle_path = Path(comparator.prediction_bundle_path)
            payload = json.loads(bundle_path.read_text(encoding="utf-8"))
            payload["model_id"] = "unlocked_secondary_model"
            forged_bundle = bundle_path.parent / "forged_secondary_model_bundle.json"
            write_json_exclusive(forged_bundle, payload)
            with self.assertRaisesRegex(FirewallError, "locked secondary comparator"):
                commit_predictions(
                    selection_lock_dir=lock_dir,
                    prediction_bundle_path=forged_bundle,
                    evaluator_state_dir=state,
                )

    def test_secondary_comparator_bundle_and_commit_are_tamper_evident(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            comparator = next(
                commit
                for commit in commits["variant_to_regulation"]
                if commit.model_role == "secondary_comparator"
            )
            commit_dir = (
                state
                / "prediction_commits"
                / comparator.prediction_commit_sha256
            )
            verified = verify_prediction_commit(commit_dir, reverify_sources=True)
            self.assertEqual(verified.model_role, "secondary_comparator")

            bundle_path = Path(comparator.prediction_bundle_path)
            payload = json.loads(bundle_path.read_text(encoding="utf-8"))
            payload["metadata"]["variant_score_transform_id"] = "forged_transform"
            bundle_path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(FirewallError, "changed after prediction commit"):
                verify_prediction_commit(commit_dir, reverify_sources=True)

    def test_power_decision_rejects_different_selected_and_baseline_splits(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            task_commits = commits[TASK_IDS[0]]
            replaced = task_commits[-1]
            valid_bundle = Path(replaced.prediction_bundle_path)
            payload = json.loads(valid_bundle.read_text(encoding="utf-8"))
            payload["split_id"] = "different-sealed-split"
            payload["source_join_key_sha256"] = canonical_hash(
                {
                    "task_id": payload["task_id"],
                    "dataset_ids": payload["dataset_ids"],
                    "split_id": payload["split_id"],
                    "row_id_field": payload["row_id_field"],
                    "unit_id_field": payload["unit_id_field"],
                    "unit_id_namespace": payload["unit_id_namespace"],
                    "biological_unit": payload["biological_unit"],
                }
            )
            forged_bundle = valid_bundle.parent / "different_split_bundle.json"
            write_json_exclusive(forged_bundle, payload)
            forged_commit = commit_predictions(
                selection_lock_dir=lock_dir,
                prediction_bundle_path=forged_bundle,
                evaluator_state_dir=state,
            )
            with self.assertRaisesRegex(FirewallError, "different sealed splits"):
                freeze_power_decision(
                    evaluator_state_dir=state,
                    prediction_commit_sha256s=(
                        commit.prediction_commit_sha256
                        for commit in (*task_commits[:-1], forged_commit)
                    ),
                    development_power_evidence_dir=root / "unused-evidence",
                )

    def test_underpowered_recomputed_evidence_stays_terminal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock_dir = create_selection(root)
            state = root / "evaluator"
            commits = commit_prediction_set(
                root,
                selection_lock_dir=lock_dir,
                evaluator_state_dir=state,
            )
            with self.assertRaisesRegex(
                ValueError, "cannot override recursively verified unit counts"
            ):
                freeze_power_set(
                    evaluator_state_dir=state,
                    commits=commits,
                    underpowered_tasks={TASK_IDS[0]},
                )

    def test_tampered_frozen_power_decision_fails_reverification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, _ = self._ready_state(root)
            power = powers[TASK_IDS[0]]
            power_dir = state / "power_decisions" / power.power_decision_sha256
            document = power_dir / "power_decision.json"
            payload = json.loads(document.read_text(encoding="utf-8"))
            payload["power_result"]["passed"] = False
            os.chmod(document, 0o640)
            document.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(FirewallError, "(?:size|checksum) mismatch"):
                verify_power_decision(power_dir)

    def test_consumption_ledger_must_be_independent_and_preexisting(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, outcome = self._ready_state(root)
            with self.assertRaisesRegex(FirewallError, "pre-existing"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=(
                        power.power_decision_sha256 for power in powers.values()
                    ),
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=root / "missing-ledger",
                )
            inside = state / "ledger"
            inside.mkdir(parents=True)
            with self.assertRaisesRegex(FirewallError, "outside evaluator state"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=(
                        power.power_decision_sha256 for power in powers.values()
                    ),
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=inside,
                )

    def test_corrupt_outcome_is_consumed_before_content_verification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, state, powers, outcome = self._ready_state(root)
            artifact = outcome / f"{TASK_IDS[0]}.tsv"
            payload = bytearray(artifact.read_bytes())
            payload[-3] = ord("X") if payload[-3] != ord("X") else ord("Y")
            artifact.write_bytes(bytes(payload))
            ledger = root / "independent-ledger"
            ledger.mkdir()
            power_ids = tuple(
                power.power_decision_sha256 for power in powers.values()
            )
            with self.assertRaisesRegex(FirewallError, "checksum mismatch"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=power_ids,
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=ledger,
                )
            self.assertEqual(len(list(ledger.iterdir())), 1)
            with self.assertRaisesRegex(FirewallError, "already consumed"):
                authorize_one_time_outcome_join(
                    evaluator_state_dir=state,
                    power_decision_sha256s=power_ids,
                    sealed_outcome_bundle_dir=outcome,
                    consumption_ledger_dir=ledger,
                )


if __name__ == "__main__":
    unittest.main()
