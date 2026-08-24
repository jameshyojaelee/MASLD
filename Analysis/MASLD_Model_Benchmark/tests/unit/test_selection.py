from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
from typing import Any, Sequence
import unittest

from masld_bench.artifacts import canonical_hash, freeze_tree, write_json_exclusive
from masld_bench.contracts import SelectionLock, _as_metadata
from masld_bench.selection import (
    _VARIANT_MANDATORY_SECONDARY_MODELS,
    _eligible_five_seed_candidates,
    CANDIDATE_STANDARD_ERROR_POLICY,
    DEVELOPMENT_ENSEMBLE_POLICY_ID,
    FINALIST_LEDGER_AUTHORITY,
    FINALIST_SEED_COUNT,
    FIT_STATE_ARTIFACT_POLICY,
    SelectionError,
    build_selection_lock,
    freeze_selection_lock,
    verify_selection_lock,
)
from masld_bench.tournament import (
    _require_comparable_score_scale,
    LEDGER_AUTHORITY_SCREENING,
    TournamentError,
)

from tests.unit._sealed_fixtures import (
    BASELINE_MODEL,
    PREDICTION_FIRST_STRESS_DATASET,
    SELECTED_MODEL,
    TASK_IDS,
    create_candidate,
    decision_document,
    digest,
)


class LockedVariantSecondaryRosterTests(unittest.TestCase):
    """Guard a comparison that could never succeed.

    ``SelectionLock`` normalizes every nested array inside a task decision to a
    tuple (``contracts._as_metadata``), so a validator that compares such a
    field against ``list(...)`` always raises, no matter what the roster holds.
    That defect blocked every fixture-dependent test in the suite while looking
    like a legitimate roster mismatch.
    """

    def test_contract_normalization_turns_nested_arrays_into_tuples(self) -> None:
        normalized = _as_metadata(
            {"mandatory_model_ids": list(_VARIANT_MANDATORY_SECONDARY_MODELS)},
            "probe",
        )
        self.assertIsInstance(normalized["mandatory_model_ids"], tuple)

    def test_normalized_roster_is_accepted_and_a_changed_one_is_not(self) -> None:
        normalized = _as_metadata(
            {"mandatory_model_ids": list(_VARIANT_MANDATORY_SECONDARY_MODELS)},
            "probe",
        )
        self.assertEqual(
            tuple(normalized["mandatory_model_ids"]),
            _VARIANT_MANDATORY_SECONDARY_MODELS,
        )
        changed = _as_metadata({"mandatory_model_ids": ["abc", "re2g"]}, "probe")
        self.assertNotEqual(
            tuple(changed["mandatory_model_ids"]),
            _VARIANT_MANDATORY_SECONDARY_MODELS,
        )


class SelectionLockTests(unittest.TestCase):
    def test_builds_strict_aggregate_lock_for_every_schedulable_task(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            lock = build_selection_lock(
                candidate=candidate, decisions=decision_document()
            )
            self.assertIsInstance(lock, SelectionLock)
            self.assertEqual(
                {decision["task_id"] for decision in lock.task_decisions},
                set(TASK_IDS),
            )
            self.assertEqual(len(lock.selected_run_ids), 5 * len(TASK_IDS))
            self.assertTrue(all(len(run_id) == 64 for run_id in lock.candidate_run_ids))
            self.assertTrue(
                all(
                    PREDICTION_FIRST_STRESS_DATASET
                    not in decision["candidate_dataset_ids"]
                    for decision in lock.task_decisions
                )
            )
            variant = next(
                decision
                for decision in lock.task_decisions
                if decision["task_id"] == "variant_to_regulation"
            )
            self.assertTrue(
                variant["variant_primary_capability"]["selected"]["capability"][
                    "primary_eligible"
                ]
            )
            self.assertTrue(
                variant["variant_primary_capability"]["baseline"]["capability"][
                    "primary_eligible"
                ]
            )
            secondary = variant["variant_secondary_evaluation"]
            # SelectionLock normalizes nested arrays to tuples, so compare
            # type-agnostically rather than against a list literal.
            self.assertEqual(
                tuple(secondary["mandatory_model_ids"]),
                ("abc", "nearest_gene", "re2g"),
            )
            self.assertEqual(
                [item["model_id"] for item in secondary["comparators"]],
                ["abc", "nearest_gene", "re2g"],
            )
            self.assertEqual(
                secondary["candidate_transform_id"],
                "absolute_signed_effect_v1",
            )
            self.assertEqual(
                secondary["primary_baseline_transform_id"],
                "absolute_signed_effect_v1",
            )
            self.assertTrue(
                all(
                    item["capability"]["primary_eligible"] is False
                    for item in secondary["comparators"]
                )
            )
            identity = lock.to_dict()
            claimed = identity.pop("lock_id")
            self.assertEqual(claimed, canonical_hash(identity))

    def test_variant_secondary_comparator_input_is_explicit_exact_and_task_local(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)

            missing = decision_document()
            variant = next(
                item
                for item in missing["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            variant.pop("variant_secondary_comparators")
            with self.assertRaisesRegex(SelectionError, "must explicitly bind"):
                build_selection_lock(candidate=candidate, decisions=missing)

            reordered = decision_document()
            variant = next(
                item
                for item in reordered["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            variant["variant_secondary_comparators"].reverse()
            with self.assertRaisesRegex(SelectionError, "canonical abc"):
                build_selection_lock(candidate=candidate, decisions=reordered)

            unknown_run = decision_document()
            variant = next(
                item
                for item in unknown_run["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            variant["variant_secondary_comparators"][0]["run_id"] = digest(
                "unknown-secondary-run"
            )
            with self.assertRaisesRegex(SelectionError, "absent from the candidate"):
                build_selection_lock(candidate=candidate, decisions=unknown_run)

            duplicate_run = decision_document()
            variant = next(
                item
                for item in duplicate_run["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            variant["variant_secondary_comparators"][1]["run_id"] = variant[
                "variant_secondary_comparators"
            ][0]["run_id"]
            with self.assertRaisesRegex(SelectionError, "repeats a run_id"):
                build_selection_lock(candidate=candidate, decisions=duplicate_run)

            wrong_task = decision_document()
            variant = next(
                item
                for item in wrong_task["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            cell = next(
                item
                for item in wrong_task["task_decisions"]
                if item["task_id"] == "cell_state_mapping"
            )
            cell["variant_secondary_comparators"] = [
                dict(item) for item in variant["variant_secondary_comparators"]
            ]
            with self.assertRaisesRegex(SelectionError, "restricted to the variant task"):
                build_selection_lock(candidate=candidate, decisions=wrong_task)

    def test_refrozen_secondary_primary_eligibility_forgery_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            lock = build_selection_lock(
                candidate=candidate, decisions=decision_document()
            )
            payload = lock.to_dict()
            variant = next(
                item
                for item in payload["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            comparator = variant["variant_secondary_evaluation"]["comparators"][0]
            comparator["capability"]["primary_eligible"] = True
            comparator["capability_sha256"] = canonical_hash(comparator["capability"])
            identity = dict(payload)
            identity.pop("lock_id")
            payload["lock_id"] = canonical_hash(identity)

            forged = root / "forged-secondary-lock"
            forged.mkdir()
            write_json_exclusive(forged / "selection_lock.json", payload)
            freeze_tree(
                forged,
                {"artifact_class": "selection_lock", "lock_id": payload["lock_id"]},
            )
            with self.assertRaisesRegex(
                SelectionError, "primary_eligible must remain exactly false"
            ):
                verify_selection_lock(forged)

    def test_refrozen_self_comparison_forgery_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            lock = build_selection_lock(
                candidate=candidate, decisions=decision_document()
            )
            payload = lock.to_dict()
            decision = payload["task_decisions"][0]
            decision["baseline_model_id"] = decision["selected_model_id"]
            identity = dict(payload)
            identity.pop("lock_id")
            payload["lock_id"] = canonical_hash(identity)

            forged = root / "forged-self-comparison-lock"
            forged.mkdir()
            write_json_exclusive(forged / "selection_lock.json", payload)
            freeze_tree(
                forged,
                {"artifact_class": "selection_lock", "lock_id": payload["lock_id"]},
            )
            with self.assertRaisesRegex(SelectionError, "must be distinct"):
                verify_selection_lock(forged)

            negative_payload = lock.to_dict()
            negative_payload["task_decisions"][0]["selected_model_id"] = (
                "shuffled_outcome"
            )
            negative_identity = dict(negative_payload)
            negative_identity.pop("lock_id")
            negative_payload["lock_id"] = canonical_hash(negative_identity)
            negative_forged = root / "forged-negative-control-lock"
            negative_forged.mkdir()
            write_json_exclusive(
                negative_forged / "selection_lock.json", negative_payload
            )
            freeze_tree(
                negative_forged,
                {
                    "artifact_class": "selection_lock",
                    "lock_id": negative_payload["lock_id"],
                },
            )
            with self.assertRaisesRegex(SelectionError, "cannot be selected"):
                verify_selection_lock(negative_forged)

            shuffled_context_payload = lock.to_dict()
            shuffled_context_payload["task_decisions"][0]["selected_model_id"] = (
                "shuffled_context"
            )
            shuffled_context_identity = dict(shuffled_context_payload)
            shuffled_context_identity.pop("lock_id")
            shuffled_context_payload["lock_id"] = canonical_hash(
                shuffled_context_identity
            )
            shuffled_context_forged = root / "forged-shuffled-context-champion"
            shuffled_context_forged.mkdir()
            write_json_exclusive(
                shuffled_context_forged / "selection_lock.json",
                shuffled_context_payload,
            )
            freeze_tree(
                shuffled_context_forged,
                {
                    "artifact_class": "selection_lock",
                    "lock_id": shuffled_context_payload["lock_id"],
                },
            )
            with self.assertRaisesRegex(SelectionError, "cannot be selected"):
                verify_selection_lock(shuffled_context_forged)

    def test_rejects_partial_task_selection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            decisions = decision_document()
            decisions["task_decisions"] = decisions["task_decisions"][:1]
            with self.assertRaisesRegex(SelectionError, "every and only schedulable"):
                build_selection_lock(candidate=candidate, decisions=decisions)

    def test_rejects_model_seed_regime_tuple_absent_from_plan(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            decisions = decision_document()
            decisions["task_decisions"][0]["selected_adaptation_regime"] = "invented"
            with self.assertRaisesRegex(SelectionError, "exactly one frozen run"):
                build_selection_lock(candidate=candidate, decisions=decisions)

            decisions = decision_document()
            decisions["task_decisions"][0]["selected_model_id"] = "invented-model"
            with self.assertRaisesRegex(SelectionError, "was not admitted"):
                build_selection_lock(candidate=candidate, decisions=decisions)

    def test_baseline_comparator_is_prespecified_and_excludes_the_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)

            self_comparison = decision_document()
            self_comparison["task_decisions"][0]["selected_model_id"] = (
                self_comparison["task_decisions"][0]["baseline_model_id"]
            )
            with self.assertRaisesRegex(SelectionError, "cannot be its own baseline"):
                build_selection_lock(candidate=candidate, decisions=self_comparison)

            unregistered_comparator = decision_document()
            unregistered_comparator["task_decisions"][0]["baseline_model_id"] = (
                SELECTED_MODEL
            )
            with self.assertRaisesRegex(
                SelectionError, "not a prespecified TaskSpec baseline"
            ):
                build_selection_lock(
                    candidate=candidate, decisions=unregistered_comparator
                )

            negative_control = decision_document()
            negative_control["task_decisions"][0]["selected_model_id"] = (
                "shuffled_outcome"
            )
            with self.assertRaisesRegex(
                SelectionError, "cannot be a selected scientific model"
            ):
                build_selection_lock(candidate=candidate, decisions=negative_control)

            shuffled_baseline = decision_document()
            shuffled_variant = next(
                item
                for item in shuffled_baseline["task_decisions"]
                if item["task_id"] == "variant_to_regulation"
            )
            shuffled_variant["baseline_model_id"] = "shuffled_context"
            with self.assertRaisesRegex(
                SelectionError, "cannot be the champion baseline comparator"
            ):
                build_selection_lock(candidate=candidate, decisions=shuffled_baseline)

            unranked_champion = decision_document(open_champion=True)
            with self.assertRaisesRegex(
                SelectionError, "requires a recursively verified scientific"
            ):
                build_selection_lock(candidate=candidate, decisions=unranked_champion)

    def test_refrozen_run_labels_must_match_the_bound_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            lock = build_selection_lock(
                candidate=candidate, decisions=decision_document()
            )
            payload = lock.to_dict()
            payload["task_decisions"][0]["baseline_runs"] = [
                dict(item)
                for item in payload["task_decisions"][0]["selected_runs"]
            ]
            identity = dict(payload)
            identity.pop("lock_id")
            payload["lock_id"] = canonical_hash(identity)

            forged = root / "forged-run-label-lock"
            forged.mkdir()
            write_json_exclusive(forged / "selection_lock.json", payload)
            freeze_tree(
                forged,
                {"artifact_class": "selection_lock", "lock_id": payload["lock_id"]},
            )
            with self.assertRaisesRegex(
                SelectionError, "run/model/regime/seed binding changed"
            ):
                verify_selection_lock(forged)

    def test_rejects_short_hashes_and_unknown_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            decisions = decision_document()
            decisions["metrics_sha256"] = "abc123"
            with self.assertRaisesRegex(SelectionError, "64-character"):
                build_selection_lock(candidate=candidate, decisions=decisions)

            decisions = decision_document()
            decisions["task_decisions"][0]["post_unblind_repair"] = True
            with self.assertRaisesRegex(SelectionError, "unknown fields"):
                build_selection_lock(candidate=candidate, decisions=decisions)

    def test_frozen_directory_reverifies_and_refrozen_forgery_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            decisions_path = root / "decisions.json"
            write_json_exclusive(decisions_path, decision_document())
            lock_dir = freeze_selection_lock(
                candidate=candidate,
                decisions_path=decisions_path,
                output_root=root / "locks",
            )
            lock = verify_selection_lock(lock_dir)
            self.assertEqual(lock.lock_id, lock_dir.name.removeprefix("selection--"))

            forged = root / "forged-lock"
            forged.mkdir()
            payload = lock.to_dict()
            payload["metrics_sha256"] = digest("forged-metrics")
            write_json_exclusive(forged / "selection_lock.json", payload)
            freeze_tree(
                forged,
                {"artifact_class": "selection_lock", "lock_id": lock.lock_id},
            )
            with self.assertRaisesRegex(SelectionError, "identity hash mismatch"):
                verify_selection_lock(forged)

    def test_mutating_frozen_lock_is_detected_before_use(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, _ = create_candidate(root)
            decisions_path = root / "decisions.json"
            write_json_exclusive(decisions_path, decision_document())
            lock_dir = freeze_selection_lock(
                candidate=candidate,
                decisions_path=decisions_path,
                output_root=root / "locks",
            )
            lock_path = lock_dir / "selection_lock.json"
            payload = json.loads(lock_path.read_text(encoding="utf-8"))
            payload["task_decisions"][0]["selected_model_id"] = SELECTED_MODEL + "-changed"
            os.chmod(lock_path, 0o640)
            lock_path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(SelectionError, "(?:size|checksum) mismatch"):
                verify_selection_lock(lock_dir)


_GUARD_TASK_ID = "cell_state_mapping"
_GUARD_SEEDS = (11, 12, 13, 14, 15)
_GUARD_REGIME = "frozen_features"
_GUARD_ALIGNMENT = digest("guard-shared-ensemble-alignment")
_GUARD_COMPARISON_IDENTITY = {
    "fold_assignment": "development_frozen_v1",
    "evaluator_id": "fixture_development_evaluator_v1",
    "input_budget": "frozen_context_window_v1",
}


def _guard_candidate(model_id: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Build one candidate that ``_eligible_five_seed_candidates`` accepts.

    Every identity field is derived the way the guard rederives it, so the
    record passes cleanly and each test below can invalidate exactly one field.
    In particular the checkpoint, task head, calibration, and fit-state bundle
    all carry the SAME digest (the composite one-fit-output-manifest policy)
    while preprocessing carries a different one (the separate prepare
    manifest); the two candidates share one comparison identity and one
    ensemble alignment.
    """

    configuration = {"model_id": model_id, "adaptation_regime": _GUARD_REGIME}
    configuration_sha256 = canonical_hash(configuration)
    candidate_id = canonical_hash(
        {
            "task_id": _GUARD_TASK_ID,
            "model_id": model_id,
            "adaptation_regime": _GUARD_REGIME,
            "candidate_configuration_sha256": configuration_sha256,
        }
    )
    fit_state_bundle = digest(f"{model_id}-fit-state-bundle")
    runs = [
        {
            "run_id": digest(f"{model_id}-run-{seed}"),
            "seed": seed,
            "task_id": _GUARD_TASK_ID,
            "model_id": model_id,
            "adaptation_regime": _GUARD_REGIME,
        }
        for seed in _GUARD_SEEDS
    ]
    run_ids = sorted(str(run["run_id"]) for run in runs)
    candidate = {
        "task_id": _GUARD_TASK_ID,
        "eligible": True,
        "model_id": model_id,
        "adaptation_regime": _GUARD_REGIME,
        "candidate_id": candidate_id,
        "candidate_configuration": configuration,
        "candidate_configuration_sha256": configuration_sha256,
        "metrics": {"absolute_primary_metric": 0.62},
        "standard_errors": {"absolute_primary_metric": 0.01},
        "standard_error_policy": CANDIDATE_STANDARD_ERROR_POLICY,
        "comparison_identity": dict(_GUARD_COMPARISON_IDENTITY),
        "comparison_identity_sha256": canonical_hash(_GUARD_COMPARISON_IDENTITY),
        "selection_artifact_ensemble_sha256s": {
            "checkpoint_sha256": fit_state_bundle,
            "task_head_sha256": fit_state_bundle,
            "calibration_sha256": fit_state_bundle,
            "preprocessing_sha256": digest(f"{model_id}-prepare-manifest"),
        },
        "fit_state_artifact_policy": FIT_STATE_ARTIFACT_POLICY,
        "fit_state_bundle_ensemble_sha256": fit_state_bundle,
        "metric_basis": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "ensemble_alignment_sha256": _GUARD_ALIGNMENT,
        "ensemble_rows_sha256": digest(f"{model_id}-ensemble-rows"),
        "run_ids": list(run_ids),
        "scheduled_run_ids": list(run_ids),
        "five_seed_universe_complete": True,
        "terminal_dispositions": [],
    }
    return candidate, runs


def _guard_plan() -> dict[str, Any]:
    """A minimal plain-Mapping plan holding a valid two-candidate ledger."""

    first, first_runs = _guard_candidate(SELECTED_MODEL)
    second, second_runs = _guard_candidate(BASELINE_MODEL)
    return {
        "selection_candidate_ledger": {
            "ledger_authority": FINALIST_LEDGER_AUTHORITY,
            "expected_seed_count": FINALIST_SEED_COUNT,
            "task_seed_sets": [
                {"task_id": _GUARD_TASK_ID, "seeds": list(_GUARD_SEEDS)}
            ],
            "model_candidates": [first, second],
            "runs": [*first_runs, *second_runs],
        }
    }


class FinalistCandidateGuardTests(unittest.TestCase):
    """Negative paths through ``_eligible_five_seed_candidates``.

    The function takes a plain ``Mapping``, so these run without the sealed
    candidate fixture.  Each test starts from a plan the function accepts
    (proved by :meth:`test_the_baseline_guard_plan_is_accepted`) and breaks one
    field, which is what keeps the assertions non-vacuous: the only difference
    between the accepted call and the raising call is the field the test names.
    """

    def test_the_baseline_guard_plan_is_accepted(self) -> None:
        candidates = _eligible_five_seed_candidates(
            plan=_guard_plan(), task_id=_GUARD_TASK_ID, seeds=_GUARD_SEEDS
        )
        self.assertIsNotNone(candidates)
        self.assertEqual(
            sorted(item["model_id"] for item in candidates),
            sorted((SELECTED_MODEL, BASELINE_MODEL)),
        )

    def test_a_candidate_not_ranked_on_the_ensemble_is_rejected(self) -> None:
        plan = _guard_plan()
        plan["selection_candidate_ledger"]["model_candidates"][0][
            "metric_basis"
        ] = "per_seed_mean_not_champion_eligible_v1"
        with self.assertRaisesRegex(
            SelectionError, "was not ranked on the five-seed development ensemble"
        ):
            _eligible_five_seed_candidates(
                plan=plan, task_id=_GUARD_TASK_ID, seeds=_GUARD_SEEDS
            )

    def test_finalist_candidates_must_share_one_alignment(self) -> None:
        plan = _guard_plan()
        plan["selection_candidate_ledger"]["model_candidates"][1][
            "ensemble_alignment_sha256"
        ] = digest("guard-divergent-ensemble-alignment")
        with self.assertRaisesRegex(
            SelectionError, "do not share one row/donor alignment"
        ):
            _eligible_five_seed_candidates(
                plan=plan, task_id=_GUARD_TASK_ID, seeds=_GUARD_SEEDS
            )

    def test_independent_fit_state_hashes_are_rejected(self) -> None:
        plan = _guard_plan()
        plan["selection_candidate_ledger"]["model_candidates"][0][
            "selection_artifact_ensemble_sha256s"
        ]["task_head_sha256"] = digest("guard-independent-task-head")
        with self.assertRaisesRegex(
            SelectionError, "claims independent fit-state artifact hashes"
        ):
            _eligible_five_seed_candidates(
                plan=plan, task_id=_GUARD_TASK_ID, seeds=_GUARD_SEEDS
            )

        plan = _guard_plan()
        candidate = plan["selection_candidate_ledger"]["model_candidates"][0]
        candidate["selection_artifact_ensemble_sha256s"]["preprocessing_sha256"] = (
            candidate["fit_state_bundle_ensemble_sha256"]
        )
        with self.assertRaisesRegex(
            SelectionError,
            "preprocessing must come from the separate prepare manifest",
        ):
            _eligible_five_seed_candidates(
                plan=plan, task_id=_GUARD_TASK_ID, seeds=_GUARD_SEEDS
            )

    def test_decision_seeds_must_equal_the_frozen_seed_authority(self) -> None:
        plan = _guard_plan()
        drifted = (*_GUARD_SEEDS[:-1], _GUARD_SEEDS[-1] + 1)
        self.assertNotEqual(drifted, _GUARD_SEEDS)
        with self.assertRaisesRegex(
            SelectionError,
            "decision seeds differ from the frozen task seed authority",
        ):
            _eligible_five_seed_candidates(
                plan=plan, task_id=_GUARD_TASK_ID, seeds=drifted
            )

    def test_a_screening_authority_ledger_yields_no_finalist_candidates(self) -> None:
        plan = _guard_plan()
        ledger = plan["selection_candidate_ledger"]
        ledger["ledger_authority"] = LEDGER_AUTHORITY_SCREENING
        ledger["expected_seed_count"] = 3
        with self.assertRaisesRegex(
            SelectionError, "finalist candidates require a finalist-authority ledger"
        ):
            _eligible_five_seed_candidates(
                plan=plan, task_id=_GUARD_TASK_ID, seeds=_GUARD_SEEDS
            )


class ScoreScaleGuardTests(unittest.TestCase):
    """Pin both edges of the five-seed raw-score comparability guard.

    ``tournament._require_comparable_score_scale`` is pure, so the whole class
    costs microseconds.  Two cases must raise and two must not, which fixes the
    guard at the frozen 4x ratio instead of merely showing that it can fire.
    """

    _FIELD = "raw_score"

    def _check(self, per_seed: Sequence[Sequence[float]]) -> None:
        _require_comparable_score_scale(
            task_id=_GUARD_TASK_ID,
            averaged_fields=(self._FIELD,),
            seeds=_GUARD_SEEDS,
            columns={
                (self._FIELD, seed): values
                for seed, values in zip(_GUARD_SEEDS, per_seed)
            },
        )

    def test_one_seed_five_times_the_others_is_refused(self) -> None:
        with self.assertRaisesRegex(
            TournamentError, r"above the prospectively frozen 4x limit"
        ):
            self._check([[0.0, 1.0]] * 4 + [[0.0, 5.0]])

    def test_a_constant_seed_beside_varying_seeds_is_refused(self) -> None:
        with self.assertRaisesRegex(
            TournamentError, "constant seed alongside varying seeds"
        ):
            self._check([[1.0, 1.0]] + [[0.0, 1.0]] * 4)

    def test_a_ratio_just_under_the_limit_is_allowed(self) -> None:
        self._check([[0.0, 1.0]] * 4 + [[0.0, 3.9]])

    def test_every_seed_constant_is_deliberately_skipped(self) -> None:
        self._check([[2.0, 2.0]] * 5)

if __name__ == "__main__":
    unittest.main()
