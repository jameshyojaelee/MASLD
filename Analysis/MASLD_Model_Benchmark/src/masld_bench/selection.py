"""Freeze complete model-selection decisions before held-back inference."""

from __future__ import annotations

import json
from math import isfinite
from pathlib import Path
import shutil
import tempfile
from typing import Any, Mapping

from .artifacts import (
    ArtifactError,
    canonical_hash,
    freeze_tree,
    publish_directory_noreplace,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)
from .contracts import ContractError, ReleaseState, SelectionLock
from .hashing import HashingError, canonicalize, is_sha256, sha256_file
from .planner import PlanningError, load_frozen_plan


class SelectionError(RuntimeError):
    """Raised when a proposed selection record is incomplete or mutable."""


SELECTION_SCHEMA_VERSION = "masld-bench-selection-lock-v1"
TERMINAL_POLICY = "no_reselection_recalibration_threshold_change_or_repair"
DEFAULT_MULTIPLICITY_PLAN = (
    "confirmatory_holm_fwer_0.05_secondary_bh_by_family"
)
BASELINE_COMPARATOR_POLICY = (
    "strongest_eligible_prespecified_baseline_excluding_selected_model_v1"
)
# Development ranks the frozen five-seed ENSEMBLE endpoint, the same estimand
# held-back inference reports.  Per-seed metrics remain as a stability diagnostic
# and never enter an ordering, so between-seed variance no longer belongs in
# the standard error: the ensemble prediction is deterministic given the fixed
# runs and has no residual seed sampling.
MODEL_SELECTION_POLICY = (
    "highest_five_seed_development_ensemble_primary_metric_"
    "report_one_standard_error_set_v2"
)
CANDIDATE_STANDARD_ERROR_POLICY = (
    "donor_or_block_bootstrap_of_the_five_seed_development_ensemble_endpoint_"
    "between_seed_variance_excluded_seed_stability_reported_separately_v2"
)
SEED_STABILITY_POLICY = (
    "four_of_five_positive_seed_direction_reported_not_ranked_v1"
)
DEVELOPMENT_ENSEMBLE_POLICY_ID = (
    "mean_prediction_five_seed_development_ensemble_v1"
)
# A fit action emits ONE read-only output manifest.  The checkpoint, task head,
# and calibration are three roles of that single composite bundle, not three
# independently produced output files, so the three fields deliberately bind the
# same digest.  Preprocessing comes from the separate prepare manifest and must
# therefore differ.  Naming the policy makes the composite claim explicit and
# machine-checked instead of implied by three identical-looking fields; a
# ledger or selection record that asserts independence is rejected.
FIT_STATE_ARTIFACT_POLICY = (
    "composite_immutable_fit_state_bundle_one_fit_output_manifest_v1"
)
INDEPENDENT_ROLE_ARTIFACT_POLICY = "independent_role_specific_output_manifests_v1"
SUPPORTED_FIT_STATE_ARTIFACT_POLICIES = frozenset({FIT_STATE_ARTIFACT_POLICY})
_COMPOSITE_FIT_STATE_ROLE_FIELDS = (
    "calibration_sha256",
    "checkpoint_sha256",
    "task_head_sha256",
)
_PREPARE_STATE_ROLE_FIELD = "preprocessing_sha256"
UNRANKED_BASELINE_COMPARATOR_POLICY = (
    "prespecified_distinct_baseline_nonchampion_contract_v1"
)
# Mirrored from tournament.LEDGER_AUTHORITY_FINALIST / _FINALIST_SEED_COUNT.
# Duplicated as literals rather than imported because selection.py is imported
# by tournament.py; the contract tests assert the two definitions agree.
FINALIST_LEDGER_AUTHORITY = (
    "finalist_campaign_universe_bound_champion_eligible_v1"
)
FINALIST_SEED_COUNT = 5
NON_SELECTABLE_NEGATIVE_CONTROL_IDS = frozenset(
    {"shuffled_context", "shuffled_outcome"}
)
CONDITIONAL_TRIGGER_INELIGIBLE_MODEL_IDS = frozenset(
    {
        "context_only",
        "sequence_only",
        "trans_only",
        *NON_SELECTABLE_NEGATIVE_CONTROL_IDS,
    }
)
OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS = frozenset(
    CONDITIONAL_TRIGGER_INELIGIBLE_MODEL_IDS
)

_DECISION_DOCUMENT_REQUIRED = frozenset({"metrics_sha256", "task_decisions"})
_DECISION_DOCUMENT_ALLOWED = frozenset(
    {
        *_DECISION_DOCUMENT_REQUIRED,
        "conditional_model_decision",
        "power_decisions",
        "multiplicity_plan",
        "metadata",
        "sealed_results_used",
    }
)

_TASK_DECISION_REQUIRED = frozenset(
    {
        "task_id",
        "selected_model_id",
        "baseline_model_id",
        "selected_adaptation_regime",
        "baseline_adaptation_regime",
        "checkpoint_sha256",
        "preprocessing_sha256",
        "task_head_sha256",
        "seeds",
        "calibration_sha256",
        "thresholds",
        "evaluator_sha256",
        "promotion_gate",
    }
)
_TASK_DECISION_ALLOWED = frozenset(
    {
        *_TASK_DECISION_REQUIRED,
        "open_champion",
        "variant_secondary_comparators",
    }
)
_LOCKED_TASK_DECISION_FIELDS = frozenset(
    {
        *_TASK_DECISION_REQUIRED,
        "open_champion",
        "selected_runs",
        "baseline_runs",
        "candidate_dataset_ids",
        "sealed_dataset_ids",
        "dataset_registry_sha256s",
        "variant_primary_capability",
        "variant_secondary_evaluation",
    }
)
_TASK_HASH_FIELDS = (
    "checkpoint_sha256",
    "preprocessing_sha256",
    "task_head_sha256",
    "calibration_sha256",
    "evaluator_sha256",
)

_VARIANT_TASK_ID = "variant_to_regulation"
_VARIANT_PRIMARY_ENDPOINT_ID = "signed_cell_type_eqtl_effect"
_VARIANT_PRIMARY_EVALUATOR_ID = (
    "variant_ld_block_fisher_z_spearman_gain_v1"
)
_VARIANT_CAPABILITY_BINDING_SCHEMA_VERSION = (
    "masld-bench-locked-variant-primary-capability-v1"
)
_VARIANT_SECONDARY_BINDING_SCHEMA_VERSION = (
    "masld-bench-locked-variant-secondary-evaluation-v1"
)
_VARIANT_SECONDARY_ENDPOINT_ID = "eqtl_retrieval"
_VARIANT_SECONDARY_EVALUATOR_ID = "variant_ld_block_eqtl_retrieval_auprc_v1"
_VARIANT_SECONDARY_SCORE_TRANSFORM_ID = "identity_link_score_v1"
_VARIANT_PRIMARY_RETRIEVAL_TRANSFORM_ID = "absolute_signed_effect_v1"
_VARIANT_MANDATORY_SECONDARY_MODELS = ("abc", "nearest_gene", "re2g")
_VARIANT_SECONDARY_EXPECTED_ROLES = {
    "abc": "link_only",
    "nearest_gene": "baseline",
    "re2g": "link_only",
}
_VARIANT_CAPABILITY_FIELDS = frozenset(
    {
        "model_id",
        "role",
        "native_outputs",
        "allowed_endpoints",
        "primary_eligible",
        "requires_fitted_head",
        "requires_observed_target_context",
        "is_mandatory_baseline",
    }
)


def _strict_mapping(
    value: object,
    *,
    required: frozenset[str],
    allowed: frozenset[str],
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise SelectionError(f"{label} must be an object")
    if any(not isinstance(key, str) for key in value):
        raise SelectionError(f"{label} keys must be strings")
    missing = sorted(required.difference(value))
    unknown = sorted(set(value).difference(allowed))
    if missing:
        raise SelectionError(f"{label} missing fields: {', '.join(missing)}")
    if unknown:
        raise SelectionError(f"{label} has unknown fields: {', '.join(unknown)}")
    return dict(value)


def _nonempty_string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise SelectionError(f"{field_name} must be a non-empty string")
    return value


def _sha256(value: object, field_name: str) -> str:
    if not is_sha256(value):
        raise SelectionError(
            f"{field_name} must be a 64-character lowercase hexadecimal SHA-256"
        )
    return str(value)


def _json_mapping(value: object, field_name: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise SelectionError(f"{field_name} must be an object")
    try:
        normalized = canonicalize(value)
    except (HashingError, TypeError, ValueError) as error:
        raise SelectionError(f"{field_name} is not canonical JSON: {error}") from error
    if not isinstance(normalized, dict):
        raise SelectionError(f"{field_name} must be an object")
    return normalized


def _load_json(path: Path) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SelectionError(f"cannot read selection decision document: {error}") from error
    if not isinstance(value, Mapping):
        raise SelectionError("selection decision document must contain a JSON object")
    return value


def _validated_plan(candidate: Path) -> tuple[dict[str, Any], str, str]:
    if not candidate.is_dir():
        raise SelectionError(f"candidate must be a frozen directory: {candidate}")
    try:
        manifest = verify_frozen_tree(candidate)
        candidate_manifest_sha256 = sha256_file(candidate / "ARTIFACTS.json")
    except (ArtifactError, OSError, ValueError, KeyError) as error:
        raise SelectionError(f"invalid frozen candidate: {error}") from error
    manifest_metadata = manifest.get("metadata")
    if not isinstance(manifest_metadata, Mapping):
        raise SelectionError("frozen candidate lacks artifact metadata")
    artifact_class = manifest_metadata.get("artifact_class")
    if artifact_class == "selection_candidate_ledger":
        try:
            from .tournament import verify_selection_candidate_ledger

            ledger = verify_selection_candidate_ledger(candidate)
        except (ImportError, AttributeError, ValueError, RuntimeError) as error:
            raise SelectionError(f"invalid selection candidate ledger: {error}") from error
        ledger_id = _sha256(ledger.get("ledger_id"), "ledger.ledger_id")
        # Only a ledger bound to an independently reviewed, pre-scoring
        # finalist campaign universe may back a selection record or a best model.  A
        # frozen-screen ledger may produce a shortlist and nothing else.
        if (
            ledger.get("ledger_authority")
            != FINALIST_LEDGER_AUTHORITY
            or ledger.get("lock_construction_allowed") is not True
            or ledger.get("champion_selection_allowed") is not True
            or ledger.get("expected_seed_count") != FINALIST_SEED_COUNT
        ):
            raise SelectionError(
                "a screening-authority selection candidate ledger cannot back "
                "a SelectionLock or a champion"
            )
        plan = {
            "campaign": {"campaign_id": ledger_id},
            "plan_sha256": ledger_id,
            "registry_snapshot_sha256": ledger["registry_snapshot_sha256"],
            "finalist_campaign_universe_binding": {
                "universe_id": ledger["reviewed_universe_id"],
                "universe_kind": ledger["campaign_universe_kind"],
                "campaign_universe_sha256": ledger["campaign_universe_sha256"],
                "manifest_sha256": ledger["campaign_universe_binding"][
                    "manifest_sha256"
                ],
                "document_sha256": ledger["campaign_universe_binding"][
                    "document_sha256"
                ],
                "selection_universe_sha256": ledger[
                    "selection_universe_sha256"
                ],
            },
            "source_lock": {
                "selection_candidate_ledger_id": ledger_id,
                "source_campaign_bindings": ledger["source_campaign_bindings"],
            },
            "resource_firewall": {
                "scientific_receipt_bindings": ledger[
                    "scientific_receipt_bindings"
                ],
                "variant_secondary_receipt_bindings": ledger[
                    "variant_secondary_receipt_bindings"
                ],
            },
            "dataset_locks": ledger["dataset_locks"],
            "model_dispositions": ledger["model_dispositions"],
            "task_dispositions": ledger["task_dispositions"],
            "runs": ledger["runs"],
            "selection_candidate_ledger": ledger,
        }
    elif artifact_class == "candidate_plan":
        try:
            plan = load_frozen_plan(candidate)
        except (PlanningError, OSError, ValueError, KeyError) as error:
            raise SelectionError(f"invalid frozen candidate plan: {error}") from error
    else:
        raise SelectionError(
            "frozen candidate must be a candidate_plan or selection_candidate_ledger"
        )

    _sha256(plan.get("plan_sha256"), "plan.plan_sha256")
    _sha256(
        plan.get("registry_snapshot_sha256"),
        "plan.registry_snapshot_sha256",
    )
    campaign = plan.get("campaign")
    if not isinstance(campaign, Mapping):
        raise SelectionError("candidate plan campaign must be an object")
    _nonempty_string(campaign.get("campaign_id"), "plan.campaign.campaign_id")

    runs = plan.get("runs")
    if not isinstance(runs, list) or not runs:
        raise SelectionError("candidate plan has no schedulable runs")
    observed_run_ids: set[str] = set()
    for index, raw_run in enumerate(runs):
        if not isinstance(raw_run, Mapping):
            raise SelectionError(f"plan.runs[{index}] must be an object")
        run = dict(raw_run)
        run_id = _sha256(run.pop("run_id", None), f"plan.runs[{index}].run_id")
        try:
            actual_run_id = canonical_hash(run)
        except (HashingError, TypeError, ValueError) as error:
            raise SelectionError(
                f"plan.runs[{index}] is not canonical JSON: {error}"
            ) from error
        if actual_run_id != run_id:
            raise SelectionError(
                f"plan.runs[{index}].run_id does not identify its full run payload"
            )
        if run_id in observed_run_ids:
            raise SelectionError(f"candidate plan repeats run_id {run_id}")
        observed_run_ids.add(run_id)
        for field_name in (
            "task_id",
            "model_id",
            "adaptation_regime",
        ):
            _nonempty_string(
                raw_run.get(field_name), f"plan.runs[{index}].{field_name}"
            )
        seed = raw_run.get("seed")
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise SelectionError(f"plan.runs[{index}].seed must be an integer")
        immutable_inputs = raw_run.get("immutable_inputs")
        if not isinstance(immutable_inputs, Mapping):
            raise SelectionError(
                f"plan.runs[{index}].immutable_inputs must be an object"
            )
        for dataset_id, registry_sha256 in immutable_inputs.items():
            _nonempty_string(dataset_id, f"plan.runs[{index}].immutable_inputs key")
            _sha256(
                registry_sha256,
                f"plan.runs[{index}].immutable_inputs[{dataset_id!r}]",
            )
    return plan, candidate_manifest_sha256, str(artifact_class)


def _model_dispositions(plan: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    raw_dispositions = plan.get("model_dispositions")
    if not isinstance(raw_dispositions, list):
        raise SelectionError("candidate plan model_dispositions must be an array")
    result: dict[str, Mapping[str, Any]] = {}
    for index, raw in enumerate(raw_dispositions):
        if not isinstance(raw, Mapping):
            raise SelectionError(f"model_dispositions[{index}] must be an object")
        model_id = _nonempty_string(
            raw.get("model_id"), f"model_dispositions[{index}].model_id"
        )
        if model_id in result:
            raise SelectionError(f"duplicate model disposition: {model_id}")
        result[model_id] = raw
    return result


def _task_baseline_roster(plan: Mapping[str, Any], task_id: str) -> tuple[str, ...]:
    raw_dispositions = plan.get("task_dispositions")
    if not isinstance(raw_dispositions, list):
        raise SelectionError("candidate plan task_dispositions must be an array")
    matches = [
        raw
        for raw in raw_dispositions
        if isinstance(raw, Mapping) and raw.get("task_id") == task_id
    ]
    if len(matches) != 1:
        raise SelectionError(
            f"candidate plan must contain exactly one task disposition for {task_id}"
        )
    raw_roster = matches[0].get("baseline_model_ids")
    if not isinstance(raw_roster, list) or not raw_roster:
        raise SelectionError(
            f"task disposition for {task_id} lacks a non-empty baseline_model_ids roster"
        )
    roster = tuple(
        _nonempty_string(value, f"{task_id}.baseline_model_ids[{index}]")
        for index, value in enumerate(raw_roster)
    )
    if len(set(roster)) != len(roster):
        raise SelectionError(f"{task_id}.baseline_model_ids repeats a model")
    return roster


def _five_seeds(value: object, task_id: str) -> tuple[int, ...]:
    if not isinstance(value, (list, tuple)):
        raise SelectionError(f"finalist {task_id} seeds must be an array")
    seeds = tuple(value)
    if any(isinstance(seed, bool) or not isinstance(seed, int) for seed in seeds):
        raise SelectionError(f"finalist {task_id} seeds must contain integers")
    if len(seeds) != 5 or len(set(seeds)) != 5:
        raise SelectionError(f"finalist {task_id} must freeze five distinct seeds")
    return tuple(sorted(seeds))


def _exact_run(
    run_index: Mapping[tuple[str, str, int, str], tuple[Mapping[str, Any], ...]],
    *,
    task_id: str,
    model_id: str,
    seed: int,
    adaptation_regime: str,
) -> Mapping[str, Any]:
    key = (task_id, model_id, seed, adaptation_regime)
    matches = run_index.get(key, ())
    if len(matches) != 1:
        raise SelectionError(
            "selection must resolve to exactly one frozen run for "
            f"task={task_id}, model={model_id}, seed={seed}, "
            f"regime={adaptation_regime}; observed {len(matches)}"
        )
    return matches[0]


def _eligible_five_seed_candidates(
    *,
    plan: Mapping[str, Any],
    task_id: str,
    seeds: tuple[int, ...],
) -> tuple[dict[str, Any], ...] | None:
    ledger = plan.get("selection_candidate_ledger")
    if not isinstance(ledger, Mapping):
        return None
    raw_candidates = ledger.get("model_candidates")
    raw_runs = ledger.get("runs")
    if not isinstance(raw_candidates, list) or not isinstance(raw_runs, list):
        raise SelectionError(
            "selection candidate ledger lacks model_candidates or runs"
        )
    if (
        ledger.get("ledger_authority") != FINALIST_LEDGER_AUTHORITY
        or ledger.get("expected_seed_count") != FINALIST_SEED_COUNT
    ):
        raise SelectionError(
            f"{task_id} finalist candidates require a finalist-authority ledger"
        )
    raw_task_seed_sets = ledger.get("task_seed_sets")
    if not isinstance(raw_task_seed_sets, list):
        raise SelectionError("selection candidate ledger lacks task_seed_sets")
    matching_seed_sets = [
        item
        for item in raw_task_seed_sets
        if isinstance(item, Mapping) and item.get("task_id") == task_id
    ]
    if len(matching_seed_sets) != 1:
        raise SelectionError(
            f"selection candidate ledger lacks one seed authority for {task_id}"
        )
    authority_seeds = _five_seeds(matching_seed_sets[0].get("seeds"), task_id)
    if seeds != authority_seeds:
        raise SelectionError(
            f"{task_id} decision seeds differ from the frozen task seed authority"
        )
    runs_by_id: dict[str, Mapping[str, Any]] = {}
    for index, raw_run in enumerate(raw_runs):
        if not isinstance(raw_run, Mapping):
            raise SelectionError(f"ledger.runs[{index}] must be an object")
        run_id = _sha256(raw_run.get("run_id"), f"ledger.runs[{index}].run_id")
        if run_id in runs_by_id:
            raise SelectionError("selection candidate ledger repeats run IDs")
        runs_by_id[run_id] = raw_run

    result: list[dict[str, Any]] = []
    for index, candidate in enumerate(raw_candidates):
        if not isinstance(candidate, Mapping):
            raise SelectionError(f"ledger.model_candidates[{index}] must be an object")
        if candidate.get("task_id") != task_id or candidate.get("eligible") is not True:
            continue
        model_id = _nonempty_string(
            candidate.get("model_id"),
            f"ledger.model_candidates[{index}].model_id",
        )
        if model_id in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
            continue
        regime = _nonempty_string(
            candidate.get("adaptation_regime"),
            f"ledger.model_candidates[{index}].adaptation_regime",
        )
        candidate_id = _sha256(
            candidate.get("candidate_id"),
            f"ledger.model_candidates[{index}].candidate_id",
        )
        candidate_configuration = _json_mapping(
            candidate.get("candidate_configuration"),
            f"candidate {candidate_id}.candidate_configuration",
        )
        candidate_configuration_sha256 = _sha256(
            candidate.get("candidate_configuration_sha256"),
            f"candidate {candidate_id}.candidate_configuration_sha256",
        )
        if (
            canonical_hash(candidate_configuration)
            != candidate_configuration_sha256
            or canonical_hash(
                {
                    "task_id": task_id,
                    "model_id": model_id,
                    "adaptation_regime": regime,
                    "candidate_configuration_sha256": (
                        candidate_configuration_sha256
                    ),
                }
            )
            != candidate_id
        ):
            raise SelectionError(
                f"candidate {candidate_id} configuration identity does not rederive"
            )
        metrics = candidate.get("metrics")
        errors = candidate.get("standard_errors")
        if not isinstance(metrics, Mapping) or not isinstance(errors, Mapping):
            raise SelectionError(f"candidate {candidate_id} lacks metrics or errors")
        primary = metrics.get("absolute_primary_metric")
        standard_error = errors.get("absolute_primary_metric")
        for label, value in (
            ("absolute_primary_metric", primary),
            ("absolute_primary_standard_error", standard_error),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise SelectionError(f"candidate {candidate_id} lacks numeric {label}")
        primary = float(primary)
        standard_error = float(standard_error)
        if (
            not isfinite(primary)
            or not isfinite(standard_error)
            or standard_error < 0.0
        ):
            raise SelectionError(
                f"candidate {candidate_id} has an invalid metric or error"
            )
        if candidate.get("standard_error_policy") != CANDIDATE_STANDARD_ERROR_POLICY:
            raise SelectionError(
                f"candidate {candidate_id} uses an unsupported standard-error policy"
            )
        comparison_identity = _json_mapping(
            candidate.get("comparison_identity"),
            f"candidate {candidate_id}.comparison_identity",
        )
        comparison_identity_sha256 = _sha256(
            candidate.get("comparison_identity_sha256"),
            f"candidate {candidate_id}.comparison_identity_sha256",
        )
        if canonical_hash(comparison_identity) != comparison_identity_sha256:
            raise SelectionError(
                f"candidate {candidate_id} comparison identity does not rederive"
            )
        raw_selection_artifacts = candidate.get(
            "selection_artifact_ensemble_sha256s"
        )
        if not isinstance(raw_selection_artifacts, Mapping) or set(
            raw_selection_artifacts
        ) != set(_TASK_HASH_FIELDS).difference({"evaluator_sha256"}):
            raise SelectionError(
                f"candidate {candidate_id} lacks exact ensemble artifact hashes"
            )
        selection_artifacts = {
            field_name: _sha256(
                raw_selection_artifacts[field_name],
                f"candidate {candidate_id}.{field_name}",
            )
            for field_name in sorted(raw_selection_artifacts)
        }
        if (
            candidate.get("fit_state_artifact_policy")
            not in SUPPORTED_FIT_STATE_ARTIFACT_POLICIES
        ):
            raise SelectionError(
                f"candidate {candidate_id} declares an unsupported fit-state "
                "artifact policy"
            )
        fit_state_bundle = _sha256(
            candidate.get("fit_state_bundle_ensemble_sha256"),
            f"candidate {candidate_id}.fit_state_bundle_ensemble_sha256",
        )
        if (
            len(
                {
                    selection_artifacts[field_name]
                    for field_name in _COMPOSITE_FIT_STATE_ROLE_FIELDS
                }
                | {fit_state_bundle}
            )
            != 1
        ):
            raise SelectionError(
                f"candidate {candidate_id} claims independent fit-state "
                f"artifact hashes; the frozen policy is {FIT_STATE_ARTIFACT_POLICY}"
            )
        if selection_artifacts[_PREPARE_STATE_ROLE_FIELD] == fit_state_bundle:
            raise SelectionError(
                f"candidate {candidate_id} preprocessing must come from the "
                "separate prepare manifest"
            )
        if candidate.get("metric_basis") != DEVELOPMENT_ENSEMBLE_POLICY_ID:
            raise SelectionError(
                f"candidate {candidate_id} was not ranked on the five-seed "
                "development ensemble"
            )
        ensemble_alignment_sha256 = _sha256(
            candidate.get("ensemble_alignment_sha256"),
            f"candidate {candidate_id}.ensemble_alignment_sha256",
        )
        ensemble_rows_sha256 = _sha256(
            candidate.get("ensemble_rows_sha256"),
            f"candidate {candidate_id}.ensemble_rows_sha256",
        )
        raw_run_ids = candidate.get("run_ids")
        if not isinstance(raw_run_ids, list) or len(raw_run_ids) != 5:
            raise SelectionError(
                f"eligible candidate {candidate_id} must have exactly five runs"
            )
        run_ids = tuple(
            _sha256(value, f"candidate {candidate_id}.run_ids[{run_index}]")
            for run_index, value in enumerate(raw_run_ids)
        )
        if len(set(run_ids)) != 5 or run_ids != tuple(sorted(run_ids)):
            raise SelectionError(
                f"candidate {candidate_id} run IDs are not unique and sorted"
            )
        if (
            candidate.get("five_seed_universe_complete") is not True
            or candidate.get("scheduled_run_ids") != list(run_ids)
            or candidate.get("terminal_dispositions") != []
        ):
            raise SelectionError(
                f"eligible candidate {candidate_id} lacks complete terminal coverage"
            )
        candidate_runs = [runs_by_id.get(run_id) for run_id in run_ids]
        if any(run is None for run in candidate_runs):
            raise SelectionError(
                f"candidate {candidate_id} names a run absent from its ledger"
            )
        observed_seeds = tuple(
            sorted(int(run["seed"]) for run in candidate_runs if run)
        )
        if observed_seeds != authority_seeds:
            raise SelectionError(
                f"candidate {candidate_id} differs from the frozen task seed authority"
            )
        if any(
            run.get("task_id") != task_id
            or run.get("model_id") != model_id
            or run.get("adaptation_regime") != regime
            for run in candidate_runs
            if run is not None
        ):
            raise SelectionError(
                f"candidate {candidate_id} run identity is inconsistent"
            )
        result.append(
            {
                "candidate_id": candidate_id,
                "model_id": model_id,
                "adaptation_regime": regime,
                "candidate_configuration_sha256": candidate_configuration_sha256,
                "absolute_primary_metric": primary,
                "absolute_primary_standard_error": standard_error,
                "comparison_identity": comparison_identity,
                "comparison_identity_sha256": comparison_identity_sha256,
                "selection_artifact_ensemble_sha256s": selection_artifacts,
                "ensemble_alignment_sha256": ensemble_alignment_sha256,
                "ensemble_rows_sha256": ensemble_rows_sha256,
                "run_ids": list(run_ids),
            }
        )
    result.sort(key=lambda item: item["candidate_id"])
    comparison_hashes = {
        item["comparison_identity_sha256"] for item in result
    }
    if len(comparison_hashes) > 1:
        raise SelectionError(
            f"{task_id} finalist candidates do not share one fold, evaluator, "
            "input budget, and row/unit comparison identity"
        )
    alignment_hashes = {item["ensemble_alignment_sha256"] for item in result}
    if len(alignment_hashes) > 1:
        raise SelectionError(
            f"{task_id} finalist candidates do not share one row/donor "
            "alignment for ensemble ranking"
        )
    return tuple(result)


def _ranked_selected_receipt(
    *,
    plan: Mapping[str, Any],
    task_id: str,
    selected_model_id: str,
    selected_regime: str,
    seeds: tuple[int, ...],
    selected_runs: tuple[Mapping[str, Any], ...] | None,
    open_champion: bool,
) -> dict[str, Any] | None:
    candidates = _eligible_five_seed_candidates(
        plan=plan, task_id=task_id, seeds=seeds
    )
    if candidates is None:
        return None
    dispositions = _model_dispositions(plan)
    eligible = [
        candidate
        for candidate in candidates
        if not open_champion
        or (
            candidate["model_id"] not in OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS
            and dispositions.get(candidate["model_id"], {}).get(
                "champion_eligible"
            )
            is True
        )
    ]
    if not eligible:
        raise SelectionError(
            f"{task_id} has no eligible five-seed candidates for the requested selection"
        )
    best = min(
        eligible,
        key=lambda item: (
            -item["absolute_primary_metric"],
            item["model_id"],
            item["adaptation_regime"],
            item["candidate_id"],
        ),
    )
    threshold = (
        best["absolute_primary_metric"]
        - best["absolute_primary_standard_error"]
    )
    within_one_se = [
        candidate
        for candidate in eligible
        if candidate["absolute_primary_metric"] >= threshold
    ]
    selected = best
    locked_runs = (
        sorted(str(run["run_id"]) for run in selected_runs)
        if selected_runs is not None
        else selected["run_ids"]
    )
    if (
        selected_model_id != selected["model_id"]
        or selected_regime != selected["adaptation_regime"]
        or locked_runs != selected["run_ids"]
    ):
        raise SelectionError(
            f"{task_id} selected model does not match the deterministic "
            "five-seed ensemble winner "
            f"{selected['model_id']}/{selected['adaptation_regime']}"
        )
    return {
        "task_id": task_id,
        "policy": MODEL_SELECTION_POLICY,
        "open_champion_subset": open_champion,
        "eligible_candidate_ids": sorted(
            candidate["candidate_id"] for candidate in eligible
        ),
        "metric_basis": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "seed_stability_policy": SEED_STABILITY_POLICY,
        "best_ensemble_candidate_id": best["candidate_id"],
        "best_ensemble_endpoint": best["absolute_primary_metric"],
        "best_ensemble_endpoint_standard_error": best[
            "absolute_primary_standard_error"
        ],
        "ensemble_alignment_sha256": selected["ensemble_alignment_sha256"],
        "ensemble_rows_sha256": selected["ensemble_rows_sha256"],
        "one_standard_error_threshold": threshold,
        "one_standard_error_candidate_ids": sorted(
            candidate["candidate_id"] for candidate in within_one_se
        ),
        "one_standard_error_policy": CANDIDATE_STANDARD_ERROR_POLICY,
        "comparison_identity": dict(selected["comparison_identity"]),
        "comparison_identity_sha256": selected["comparison_identity_sha256"],
        "selection_artifact_ensemble_sha256s": dict(
            selected["selection_artifact_ensemble_sha256s"]
        ),
        "selected_candidate_id": selected["candidate_id"],
        "selected_model_id": selected["model_id"],
        "selected_adaptation_regime": selected["adaptation_regime"],
        "selected_run_ids": list(selected["run_ids"]),
    }


def _ranked_baseline_receipt(
    *,
    plan: Mapping[str, Any],
    task_id: str,
    selected_model_id: str,
    baseline_model_id: str,
    baseline_regime: str,
    baseline_roster: tuple[str, ...],
    seeds: tuple[int, ...],
    baseline_runs: tuple[Mapping[str, Any], ...] | None,
) -> dict[str, Any] | None:
    """Prove the strongest candidate-excluding baseline from a scientific ledger."""

    candidates = _eligible_five_seed_candidates(
        plan=plan, task_id=task_id, seeds=seeds
    )
    if candidates is None:
        return None
    ranked = [
        candidate
        for candidate in candidates
        if candidate["model_id"] in baseline_roster
        and candidate["model_id"] != selected_model_id
    ]
    if not ranked:
        raise SelectionError(
            f"{task_id} has no eligible five-seed prespecified baseline after "
            "excluding the selected model and negative controls"
        )
    strongest = min(
        ranked,
        key=lambda item: (
            -item["absolute_primary_metric"],
            item["model_id"],
            item["adaptation_regime"],
            item["candidate_id"],
        ),
    )
    locked_runs = (
        sorted(str(run["run_id"]) for run in baseline_runs)
        if baseline_runs is not None
        else strongest["run_ids"]
    )
    if (
        baseline_model_id != strongest["model_id"]
        or baseline_regime != strongest["adaptation_regime"]
        or locked_runs != strongest["run_ids"]
    ):
        raise SelectionError(
            f"{task_id} baseline does not match the {BASELINE_COMPARATOR_POLICY} "
            f"winner {strongest['model_id']}/{strongest['adaptation_regime']}"
        )
    return {
        "task_id": task_id,
        "policy": BASELINE_COMPARATOR_POLICY,
        "candidate_id": strongest["candidate_id"],
        "model_id": strongest["model_id"],
        "adaptation_regime": strongest["adaptation_regime"],
        "absolute_primary_metric": strongest["absolute_primary_metric"],
        "metric_basis": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "comparison_identity_sha256": strongest[
            "comparison_identity_sha256"
        ],
        "selection_artifact_ensemble_sha256s": dict(
            strongest["selection_artifact_ensemble_sha256s"]
        ),
        "run_ids": list(strongest["run_ids"]),
        "excluded_selected_model_id": selected_model_id,
    }


def _validated_primary_variant_capability_record(
    value: object, *, label: str, model_id: str
) -> dict[str, Any]:
    capability = _strict_mapping(
        value,
        required=_VARIANT_CAPABILITY_FIELDS,
        allowed=_VARIANT_CAPABILITY_FIELDS,
        label=label,
    )
    try:
        capability = canonicalize(capability)
    except (HashingError, TypeError, ValueError) as error:
        raise SelectionError(
            f"{label} is not canonical JSON: {error}"
        ) from error
    if not isinstance(capability, dict):
        raise SelectionError(f"{label} must be an object")
    if capability.get("model_id") != model_id:
        raise SelectionError(f"{label} does not identify model {model_id}")
    if capability.get("role") not in {
        "baseline",
        "conditional_only",
        "primary_candidate",
    }:
        raise SelectionError(f"{label} has an invalid primary capability role")
    for field_name in ("native_outputs", "allowed_endpoints"):
        values = capability.get(field_name)
        if (
            not isinstance(values, list)
            or not values
            or any(not isinstance(item, str) or not item for item in values)
            or values != sorted(set(values))
        ):
            raise SelectionError(
                f"{label}.{field_name} must be non-empty, unique, and sorted"
            )
    for field_name in (
        "primary_eligible",
        "requires_fitted_head",
        "requires_observed_target_context",
        "is_mandatory_baseline",
    ):
        if not isinstance(capability.get(field_name), bool):
            raise SelectionError(f"{label}.{field_name} must be boolean")
    if capability["primary_eligible"] is not True:
        raise SelectionError(f"{label} is ineligible for the signed primary endpoint")
    if _VARIANT_PRIMARY_ENDPOINT_ID not in capability["allowed_endpoints"]:
        raise SelectionError(f"{label} does not allow the signed primary endpoint")
    if capability["requires_observed_target_context"] is not False:
        raise SelectionError(f"{label} requires unavailable sealed target context")
    if not set(capability["native_outputs"]).intersection(
        {"gene_expression_delta", "rna_coverage_delta"}
    ):
        raise SelectionError(f"{label} lacks a signed gene-level output")
    return capability


def _variant_capability_for_role(
    *,
    task_id: str,
    role: str,
    model_id: str,
    runs: tuple[Mapping[str, Any], ...],
) -> tuple[dict[str, Any], str, str]:
    observed: list[tuple[dict[str, Any], str, str]] = []
    for index, run in enumerate(runs):
        metadata = run.get("metadata")
        if not isinstance(metadata, Mapping):
            raise SelectionError(
                f"{task_id} {role} run {index} lacks immutable run metadata"
            )
        capability = _validated_primary_variant_capability_record(
            metadata.get("variant_capability"),
            label=f"{task_id}.{role}_runs[{index}].variant_capability",
            model_id=model_id,
        )
        if metadata.get("primary_endpoint_scoring_allowed") is not True:
            raise SelectionError(
                f"{task_id} {role} run forbids primary endpoint scoring"
            )
        capability_sha256 = _sha256(
            metadata.get("variant_capability_sha256"),
            f"{task_id}.{role}_runs[{index}].variant_capability_sha256",
        )
        if capability_sha256 != canonical_hash(capability):
            raise SelectionError(
                f"{task_id} {role} run capability SHA-256 does not match its record"
            )
        registry_sha256 = _sha256(
            metadata.get("variant_capability_registry_sha256"),
            f"{task_id}.{role}_runs[{index}].variant_capability_registry_sha256",
        )
        if metadata.get("primary_evaluator_id") != _VARIANT_PRIMARY_EVALUATOR_ID:
            raise SelectionError(
                f"{task_id} {role} run uses the wrong primary evaluator"
            )
        observed.append((capability, capability_sha256, registry_sha256))
    if any(item != observed[0] for item in observed[1:]):
        raise SelectionError(
            f"{task_id} {role} capability binding differs across frozen seeds"
        )
    return observed[0]


def _locked_variant_primary_capability(
    *,
    task_id: str,
    selected_model_id: str,
    baseline_model_id: str,
    selected_runs: tuple[Mapping[str, Any], ...],
    baseline_runs: tuple[Mapping[str, Any], ...],
) -> dict[str, Any] | None:
    if task_id != _VARIANT_TASK_ID:
        return None
    selected, selected_sha256, selected_registry_sha256 = (
        _variant_capability_for_role(
            task_id=task_id,
            role="selected",
            model_id=selected_model_id,
            runs=selected_runs,
        )
    )
    baseline, baseline_sha256, baseline_registry_sha256 = (
        _variant_capability_for_role(
            task_id=task_id,
            role="baseline",
            model_id=baseline_model_id,
            runs=baseline_runs,
        )
    )
    if selected_registry_sha256 != baseline_registry_sha256:
        raise SelectionError(
            f"{task_id} selected and baseline runs use different capability registries"
        )
    return {
        "schema_version": _VARIANT_CAPABILITY_BINDING_SCHEMA_VERSION,
        "task_id": task_id,
        "primary_endpoint_id": _VARIANT_PRIMARY_ENDPOINT_ID,
        "primary_evaluator_id": _VARIANT_PRIMARY_EVALUATOR_ID,
        "capability_registry_sha256": selected_registry_sha256,
        "selected": {
            "model_id": selected_model_id,
            "capability": selected,
            "capability_sha256": selected_sha256,
        },
        "baseline": {
            "model_id": baseline_model_id,
            "capability": baseline,
            "capability_sha256": baseline_sha256,
        },
    }


def _validated_secondary_variant_capability_record(
    value: object, *, label: str, model_id: str
) -> dict[str, Any]:
    capability = _strict_mapping(
        value,
        required=_VARIANT_CAPABILITY_FIELDS,
        allowed=_VARIANT_CAPABILITY_FIELDS,
        label=label,
    )
    try:
        capability = canonicalize(capability)
    except (HashingError, TypeError, ValueError) as error:
        raise SelectionError(f"{label} is not canonical JSON: {error}") from error
    if not isinstance(capability, dict):
        raise SelectionError(f"{label} must be an object")
    if capability.get("model_id") != model_id:
        raise SelectionError(f"{label} does not identify model {model_id}")
    if model_id not in _VARIANT_SECONDARY_EXPECTED_ROLES:
        raise SelectionError(f"{label} names an unsupported secondary comparator")
    if capability.get("role") != _VARIANT_SECONDARY_EXPECTED_ROLES[model_id]:
        raise SelectionError(f"{label} has the wrong endpoint-specific role")
    if capability.get("native_outputs") != ["enhancer_gene_link_score"]:
        raise SelectionError(f"{label} must expose only an enhancer-gene link score")
    if capability.get("allowed_endpoints") != [
        "enhancer_gene_link",
        _VARIANT_SECONDARY_ENDPOINT_ID,
    ]:
        raise SelectionError(f"{label} has the wrong secondary endpoint scope")
    expected_flags = {
        "primary_eligible": False,
        "requires_fitted_head": False,
        "requires_observed_target_context": False,
        "is_mandatory_baseline": True,
    }
    for field_name, expected in expected_flags.items():
        if capability.get(field_name) is not expected:
            raise SelectionError(
                f"{label}.{field_name} must remain exactly {str(expected).lower()}"
            )
    return capability


def _locked_variant_secondary_evaluation(
    *,
    task_id: str,
    field_provided: bool,
    raw_comparators: object,
    runs_by_id: Mapping[str, Mapping[str, Any]],
    plan: Mapping[str, Any],
    dispositions: Mapping[str, Mapping[str, Any]],
    primary_binding: Mapping[str, Any] | None,
    expected_immutable_inputs: Mapping[str, Any],
) -> dict[str, Any] | None:
    if task_id != _VARIANT_TASK_ID:
        if field_provided:
            raise SelectionError(
                f"{task_id}.variant_secondary_comparators is restricted to the variant task"
            )
        return None
    if not field_provided:
        raise SelectionError(
            f"{task_id}.variant_secondary_comparators must explicitly bind "
            "abc, nearest_gene, and re2g"
        )
    if not isinstance(raw_comparators, (list, tuple)):
        raise SelectionError(
            f"{task_id}.variant_secondary_comparators must be an array"
        )
    parsed_inputs: list[dict[str, str]] = []
    for index, value in enumerate(raw_comparators):
        parsed = _strict_mapping(
            value,
            required=frozenset({"model_id", "run_id"}),
            allowed=frozenset({"model_id", "run_id"}),
            label=f"{task_id}.variant_secondary_comparators[{index}]",
        )
        parsed_inputs.append(
            {
                "model_id": _nonempty_string(
                    parsed["model_id"],
                    f"{task_id}.variant_secondary_comparators[{index}].model_id",
                ),
                "run_id": _sha256(
                    parsed["run_id"],
                    f"{task_id}.variant_secondary_comparators[{index}].run_id",
                ),
            }
        )
    observed_models = tuple(item["model_id"] for item in parsed_inputs)
    if observed_models != _VARIANT_MANDATORY_SECONDARY_MODELS:
        raise SelectionError(
            f"{task_id}.variant_secondary_comparators must contain exactly the "
            "canonical abc, nearest_gene, and re2g roster"
        )
    if len({item["run_id"] for item in parsed_inputs}) != len(parsed_inputs):
        raise SelectionError(
            f"{task_id}.variant_secondary_comparators repeats a run_id"
        )

    ledger = plan.get("selection_candidate_ledger")
    auxiliary_by_run: dict[str, Mapping[str, Any]] | None = None
    scientific_run_ids: set[str] = set()
    if isinstance(ledger, Mapping):
        raw_scientific = ledger.get("scientific_runs")
        if not isinstance(raw_scientific, list):
            raise SelectionError("selection candidate ledger lacks scientific runs")
        scientific_run_ids = {
            str(record.get("run_id"))
            for record in raw_scientific
            if isinstance(record, Mapping)
        }
        if len(scientific_run_ids) != len(raw_scientific):
            raise SelectionError("selection candidate ledger repeats scientific run IDs")
        raw_auxiliary = ledger.get("variant_secondary_runs")
        if not isinstance(raw_auxiliary, list):
            raise SelectionError(
                "selection candidate ledger lacks variant secondary runs"
            )
        auxiliary_by_run = {
            str(record.get("run_id")): record
            for record in raw_auxiliary
            if isinstance(record, Mapping)
        }
        if len(auxiliary_by_run) != len(raw_auxiliary):
            raise SelectionError(
                "selection candidate ledger repeats variant secondary run IDs"
            )

    locked_comparators: list[dict[str, Any]] = []
    registry_sha256: str | None = None
    for index, requested in enumerate(parsed_inputs):
        model_id = requested["model_id"]
        run_id = requested["run_id"]
        disposition = dispositions.get(model_id)
        if disposition is None or disposition.get("disposition") == "blocked":
            raise SelectionError(f"secondary comparator was not admitted: {model_id}")
        run = runs_by_id.get(run_id)
        if run is None:
            raise SelectionError(
                f"secondary comparator run is absent from the candidate: {run_id}"
            )
        if run.get("task_id") != task_id or run.get("model_id") != model_id:
            raise SelectionError(
                f"secondary comparator run {run_id} does not bind {task_id}/{model_id}"
            )
        if dict(run.get("immutable_inputs", {})) != dict(expected_immutable_inputs):
            raise SelectionError(
                f"secondary comparator {model_id} uses different immutable inputs"
            )
        metadata = run.get("metadata")
        if not isinstance(metadata, Mapping):
            raise SelectionError(
                f"secondary comparator run {run_id} lacks immutable capability metadata"
            )
        capability = _validated_secondary_variant_capability_record(
            metadata.get("variant_capability"),
            label=(
                f"{task_id}.variant_secondary_comparators[{index}].variant_capability"
            ),
            model_id=model_id,
        )
        capability_sha256 = _sha256(
            metadata.get("variant_capability_sha256"),
            f"secondary comparator {model_id} capability SHA-256",
        )
        if capability_sha256 != canonical_hash(capability):
            raise SelectionError(
                f"secondary comparator {model_id} capability SHA-256 changed"
            )
        observed_registry_sha256 = _sha256(
            metadata.get("variant_capability_registry_sha256"),
            f"secondary comparator {model_id} capability registry SHA-256",
        )
        if registry_sha256 is None:
            registry_sha256 = observed_registry_sha256
        elif observed_registry_sha256 != registry_sha256:
            raise SelectionError(
                "secondary comparator runs use different capability registries"
            )
        if metadata.get("primary_endpoint_scoring_allowed") is not False:
            raise SelectionError(
                f"secondary comparator {model_id} must forbid primary endpoint scoring"
            )
        if metadata.get("primary_evaluator_id") != _VARIANT_PRIMARY_EVALUATOR_ID:
            raise SelectionError(
                f"secondary comparator {model_id} names the wrong task evaluator"
            )
        if auxiliary_by_run is not None:
            auxiliary = auxiliary_by_run.get(run_id)
            if auxiliary is None:
                raise SelectionError(
                    f"secondary comparator lacks a verified auxiliary receipt: {run_id}"
                )
            if (
                auxiliary.get("task_id") != task_id
                or auxiliary.get("model_id") != model_id
                or auxiliary.get("variant_capability") != capability
                or auxiliary.get("variant_capability_sha256")
                != capability_sha256
                or auxiliary.get("variant_capability_registry_sha256")
                != observed_registry_sha256
                or dict(auxiliary.get("immutable_inputs", {}))
                != dict(expected_immutable_inputs)
                or auxiliary.get("endpoint_id") != _VARIANT_SECONDARY_ENDPOINT_ID
                or auxiliary.get("evaluator_id") != _VARIANT_SECONDARY_EVALUATOR_ID
                or auxiliary.get("score_transform_id")
                != _VARIANT_SECONDARY_SCORE_TRANSFORM_ID
                or auxiliary.get("outcomes_read") is not False
                or auxiliary.get("metrics_computed") is not False
                or auxiliary.get("auxiliary_status") != "succeeded_non_scoring"
            ):
                raise SelectionError(
                    f"secondary comparator auxiliary receipt changed for {run_id}"
                )
            _sha256(
                auxiliary.get("source_secondary_run_receipt_id"),
                f"secondary comparator {model_id} auxiliary receipt ID",
            )
            if run_id in scientific_run_ids:
                raise SelectionError(
                    f"secondary comparator entered primary scientific scoring: {run_id}"
                )
        locked_comparators.append(
            {
                "model_id": model_id,
                "run_id": run_id,
                "capability": capability,
                "capability_sha256": capability_sha256,
                "score_transform_id": _VARIANT_SECONDARY_SCORE_TRANSFORM_ID,
            }
        )
    if registry_sha256 is None:
        raise SelectionError("variant secondary comparator roster is empty")
    if (
        not isinstance(primary_binding, Mapping)
        or primary_binding.get("capability_registry_sha256") != registry_sha256
    ):
        raise SelectionError(
            "primary and secondary variant runs use different capability registries"
        )
    return {
        "schema_version": _VARIANT_SECONDARY_BINDING_SCHEMA_VERSION,
        "task_id": task_id,
        "endpoint_id": _VARIANT_SECONDARY_ENDPOINT_ID,
        "evaluator_id": _VARIANT_SECONDARY_EVALUATOR_ID,
        "capability_registry_sha256": registry_sha256,
        "candidate_transform_id": _VARIANT_PRIMARY_RETRIEVAL_TRANSFORM_ID,
        "primary_baseline_transform_id": _VARIANT_PRIMARY_RETRIEVAL_TRANSFORM_ID,
        "mandatory_model_ids": list(_VARIANT_MANDATORY_SECONDARY_MODELS),
        "comparators": locked_comparators,
    }


def _dataset_bindings(
    *,
    task_id: str,
    selected_runs: tuple[Mapping[str, Any], ...],
    baseline_runs: tuple[Mapping[str, Any], ...],
    plan: Mapping[str, Any],
) -> tuple[tuple[str, ...], tuple[str, ...], dict[str, str]]:
    all_runs = (*selected_runs, *baseline_runs)
    first_inputs = dict(all_runs[0].get("immutable_inputs", {}))
    if any(dict(run.get("immutable_inputs", {})) != first_inputs for run in all_runs[1:]):
        raise SelectionError(
            f"selected and baseline runs for {task_id} do not share immutable inputs"
        )
    dataset_locks = plan.get("dataset_locks")
    if not isinstance(dataset_locks, Mapping):
        raise SelectionError("candidate plan dataset_locks must be an object")
    task_dispositions = plan.get("task_dispositions")
    if not isinstance(task_dispositions, list):
        raise SelectionError("candidate plan task_dispositions must be an array")
    matches = [
        item
        for item in task_dispositions
        if isinstance(item, Mapping) and item.get("task_id") == task_id
    ]
    if len(matches) != 1:
        raise SelectionError(f"candidate plan must bind one TaskSpec for {task_id}")
    task = matches[0]

    def declared(field: str) -> tuple[str, ...]:
        value = task.get(field)
        if not isinstance(value, (list, tuple)):
            raise SelectionError(f"{task_id}.{field} must be an array")
        identifiers = tuple(
            _nonempty_string(item, f"{task_id}.{field}[]") for item in value
        )
        if len(set(identifiers)) != len(identifiers):
            raise SelectionError(f"{task_id}.{field} cannot contain duplicates")
        return identifiers

    training = declared("datasets_train")
    development = declared("datasets_development")
    declared_sealed = declared("datasets_sealed")
    if set(training) & set(development) or (
        set((*training, *development)) & set(declared_sealed)
    ):
        raise SelectionError(f"{task_id} TaskSpec dataset roles must not overlap")
    selection_development: list[str] = []
    prediction_first_stress: list[str] = []
    for dataset_id in development:
        raw_lock = dataset_locks.get(dataset_id)
        if not isinstance(raw_lock, Mapping):
            raise SelectionError(
                f"TaskSpec for {task_id} names dataset absent from dataset_locks: {dataset_id}"
            )
        if raw_lock.get("prediction_first_policy") is None:
            selection_development.append(dataset_id)
        else:
            prediction_first_stress.append(dataset_id)

    expected_run_inputs = set((*training, *selection_development))
    if set(first_inputs).intersection(prediction_first_stress):
        raise SelectionError(
            f"{task_id} selection evidence includes a prediction-first stress dataset"
        )
    if set(first_inputs) != expected_run_inputs:
        raise SelectionError(
            f"{task_id} run inputs must exactly match declared fit and "
            "model-selection development datasets; prediction-first stress "
            "datasets are forbidden"
        )

    sealed: list[str] = []
    for dataset_id in declared_sealed:
        raw_lock = dataset_locks.get(dataset_id)
        if not isinstance(raw_lock, Mapping):
            raise SelectionError(
                f"TaskSpec for {task_id} names dataset absent from dataset_locks: {dataset_id}"
            )
        if (
            raw_lock.get("role") == "withheld_sealed"
            and raw_lock.get("status") == "withheld_sealed"
        ):
            sealed.append(dataset_id)
    candidate_ids = tuple(sorted((*training, *selection_development, *sealed)))
    bindings: dict[str, str] = {}
    for dataset_id in candidate_ids:
        raw_lock = dataset_locks.get(dataset_id)
        if not isinstance(raw_lock, Mapping):
            raise SelectionError(
                f"TaskSpec for {task_id} names dataset absent from dataset_locks: {dataset_id}"
            )
        locked_registry_sha256 = _sha256(
            raw_lock.get("registry_sha256"),
            f"dataset_locks[{dataset_id!r}].registry_sha256",
        )
        bindings[dataset_id] = locked_registry_sha256
        if dataset_id in first_inputs:
            activation_sha256 = _sha256(
                raw_lock.get("activation_sha256"),
                f"dataset_locks[{dataset_id!r}].activation_sha256",
            )
            observed_activation = _sha256(
                first_inputs[dataset_id],
                f"run activation binding {task_id}/{dataset_id}",
            )
            if activation_sha256 != observed_activation:
                raise SelectionError(
                    f"run and plan disagree on activation SHA-256 for dataset {dataset_id}"
                )
    return candidate_ids, tuple(sorted(sealed)), bindings


def _normalize_task_decision(
    raw: Mapping[str, Any],
    *,
    plan: Mapping[str, Any],
    run_index: Mapping[tuple[str, str, int, str], tuple[Mapping[str, Any], ...]],
    runs_by_id: Mapping[str, Mapping[str, Any]],
    known_tasks: set[str],
    dispositions: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    decision = _strict_mapping(
        raw,
        required=_TASK_DECISION_REQUIRED,
        allowed=_TASK_DECISION_ALLOWED,
        label="task_decision",
    )
    task_id = _nonempty_string(decision["task_id"], "task_decision.task_id")
    if task_id not in known_tasks:
        raise SelectionError(f"selection names task absent from campaign runs: {task_id}")
    selected_model_id = _nonempty_string(
        decision["selected_model_id"], f"{task_id}.selected_model_id"
    )
    if selected_model_id in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
        raise SelectionError(
            f"negative control {selected_model_id} cannot be a selected scientific model"
        )
    baseline_model_id = _nonempty_string(
        decision["baseline_model_id"], f"{task_id}.baseline_model_id"
    )
    baseline_roster = _task_baseline_roster(plan, task_id)
    if baseline_model_id in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
        raise SelectionError(
            f"negative control {baseline_model_id} cannot be the champion baseline comparator"
        )
    if baseline_model_id not in baseline_roster:
        raise SelectionError(
            f"{task_id}.baseline_model_id is not a prespecified TaskSpec baseline"
        )
    if selected_model_id == baseline_model_id:
        raise SelectionError(
            f"{task_id} must use the {BASELINE_COMPARATOR_POLICY} policy; "
            "the selected model cannot be its own baseline comparator"
        )
    selected_regime = _nonempty_string(
        decision["selected_adaptation_regime"],
        f"{task_id}.selected_adaptation_regime",
    )
    baseline_regime = _nonempty_string(
        decision["baseline_adaptation_regime"],
        f"{task_id}.baseline_adaptation_regime",
    )
    for role, model_id in (
        ("selected", selected_model_id),
        ("baseline", baseline_model_id),
    ):
        disposition = dispositions.get(model_id)
        if disposition is None or disposition.get("disposition") == "blocked":
            raise SelectionError(f"{role} model was not admitted: {model_id}")

    open_champion = decision.get("open_champion", False)
    if not isinstance(open_champion, bool):
        raise SelectionError(f"{task_id}.open_champion must be boolean")
    if open_champion and not bool(
        dispositions[selected_model_id].get("champion_eligible", False)
    ):
        raise SelectionError(
            f"model {selected_model_id} is ineligible for an open champion claim"
        )
    if open_champion and selected_model_id in OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS:
        raise SelectionError(
            f"control {selected_model_id} cannot receive an open champion claim"
        )
    if open_champion and not isinstance(
        plan.get("selection_candidate_ledger"), Mapping
    ):
        raise SelectionError(
            "an open champion requires a recursively verified scientific selection "
            "candidate ledger with development metrics"
        )

    seeds = _five_seeds(decision["seeds"], task_id)
    ledger_backed = isinstance(plan.get("selection_candidate_ledger"), Mapping)
    if ledger_backed:
        selected_ranking_receipt = _ranked_selected_receipt(
            plan=plan,
            task_id=task_id,
            selected_model_id=selected_model_id,
            selected_regime=selected_regime,
            seeds=seeds,
            selected_runs=None,
            open_champion=open_champion,
        )
        baseline_ranking_receipt = _ranked_baseline_receipt(
            plan=plan,
            task_id=task_id,
            selected_model_id=selected_model_id,
            baseline_model_id=baseline_model_id,
            baseline_regime=baseline_regime,
            baseline_roster=baseline_roster,
            seeds=seeds,
            baseline_runs=None,
        )
        if selected_ranking_receipt is None or baseline_ranking_receipt is None:
            raise SelectionError(
                f"{task_id} ledger does not yield ranked selected and baseline runs"
            )
        selected_runs = tuple(
            sorted(
                (
                    runs_by_id[run_id]
                    for run_id in selected_ranking_receipt["selected_run_ids"]
                ),
                key=lambda run: int(run["seed"]),
            )
        )
        baseline_runs = tuple(
            sorted(
                (
                    runs_by_id[run_id]
                    for run_id in baseline_ranking_receipt["run_ids"]
                ),
                key=lambda run: int(run["seed"]),
            )
        )
    else:
        selected_runs = tuple(
            _exact_run(
                run_index,
                task_id=task_id,
                model_id=selected_model_id,
                seed=seed,
                adaptation_regime=selected_regime,
            )
            for seed in seeds
        )
        baseline_runs = tuple(
            _exact_run(
                run_index,
                task_id=task_id,
                model_id=baseline_model_id,
                seed=seed,
                adaptation_regime=baseline_regime,
            )
            for seed in seeds
        )
        selected_ranking_receipt = None
    if selected_ranking_receipt is not None:
        expected_artifacts = selected_ranking_receipt[
            "selection_artifact_ensemble_sha256s"
        ]
        for field_name in (
            "checkpoint_sha256",
            "preprocessing_sha256",
            "task_head_sha256",
            "calibration_sha256",
        ):
            if decision[field_name] != expected_artifacts[field_name]:
                raise SelectionError(
                    f"{task_id}.{field_name} differs from the five-seed "
                    "scientific artifact ensemble"
                )
    candidate_dataset_ids, sealed_dataset_ids, dataset_bindings = _dataset_bindings(
        task_id=task_id,
        selected_runs=selected_runs,
        baseline_runs=baseline_runs,
        plan=plan,
    )

    for field_name in _TASK_HASH_FIELDS:
        _sha256(decision[field_name], f"{task_id}.{field_name}")
    ledger = plan.get("selection_candidate_ledger")
    if isinstance(ledger, Mapping):
        raw_scientific = ledger.get("scientific_runs")
        if not isinstance(raw_scientific, list):
            raise SelectionError("selection candidate ledger lacks scientific runs")
        science_by_run = {
            str(record.get("run_id")): record
            for record in raw_scientific
            if isinstance(record, Mapping)
        }
        required_runs = (*selected_runs, *baseline_runs)
        if len(science_by_run) != len(raw_scientific):
            raise SelectionError("selection candidate ledger repeats scientific run IDs")
        for run in required_runs:
            scientific = science_by_run.get(str(run["run_id"]))
            if scientific is None:
                raise SelectionError(
                    f"selected run lacks a verified scientific receipt: {run['run_id']}"
                )
            if scientific.get("endpoint_evaluator_sha256") != decision["evaluator_sha256"]:
                raise SelectionError(
                    f"{task_id}.evaluator_sha256 differs from selected scientific receipts"
                )
        receipt_bindings = ledger.get("scientific_receipt_bindings")
        if not isinstance(receipt_bindings, list):
            raise SelectionError(
                "selection candidate ledger lacks scientific receipt bindings"
            )
        for selected_run, baseline_run in zip(
            selected_runs, baseline_runs, strict=True
        ):
            matches = [
                binding
                for binding in receipt_bindings
                if isinstance(binding, Mapping)
                and binding.get("run_id") == selected_run["run_id"]
                and binding.get("baseline_run_id") == baseline_run["run_id"]
            ]
            if len(matches) != 1:
                raise SelectionError(
                    f"{task_id} requires exactly one paired scientific receipt for "
                    f"selected run {selected_run['run_id']} and baseline run "
                    f"{baseline_run['run_id']}"
                )
    thresholds = _json_mapping(decision["thresholds"], f"{task_id}.thresholds")
    promotion_gate = _nonempty_string(
        decision["promotion_gate"], f"{task_id}.promotion_gate"
    )
    variant_primary_capability = _locked_variant_primary_capability(
        task_id=task_id,
        selected_model_id=selected_model_id,
        baseline_model_id=baseline_model_id,
        selected_runs=selected_runs,
        baseline_runs=baseline_runs,
    )
    variant_secondary_evaluation = _locked_variant_secondary_evaluation(
        task_id=task_id,
        field_provided="variant_secondary_comparators" in decision,
        raw_comparators=decision.get("variant_secondary_comparators"),
        runs_by_id=runs_by_id,
        plan=plan,
        dispositions=dispositions,
        primary_binding=variant_primary_capability,
        expected_immutable_inputs=dict(selected_runs[0]["immutable_inputs"]),
    )

    return {
        "task_id": task_id,
        "selected_model_id": selected_model_id,
        "baseline_model_id": baseline_model_id,
        "selected_adaptation_regime": selected_regime,
        "baseline_adaptation_regime": baseline_regime,
        "checkpoint_sha256": decision["checkpoint_sha256"],
        "preprocessing_sha256": decision["preprocessing_sha256"],
        "task_head_sha256": decision["task_head_sha256"],
        "seeds": list(seeds),
        "calibration_sha256": decision["calibration_sha256"],
        "thresholds": thresholds,
        "evaluator_sha256": decision["evaluator_sha256"],
        "promotion_gate": promotion_gate,
        "open_champion": open_champion,
        "selected_runs": [
            {"seed": seed, "run_id": str(run["run_id"])}
            for seed, run in zip(seeds, selected_runs, strict=True)
        ],
        "baseline_runs": [
            {"seed": seed, "run_id": str(run["run_id"])}
            for seed, run in zip(seeds, baseline_runs, strict=True)
        ],
        "candidate_dataset_ids": list(candidate_dataset_ids),
        "sealed_dataset_ids": list(sealed_dataset_ids),
        "dataset_registry_sha256s": dataset_bindings,
        "variant_primary_capability": variant_primary_capability,
        "variant_secondary_evaluation": variant_secondary_evaluation,
    }


def _validate_locked_variant_primary_capability(
    *, task_id: str, decision: Mapping[str, Any]
) -> None:
    raw_binding = decision.get("variant_primary_capability")
    if task_id != _VARIANT_TASK_ID:
        if raw_binding is not None:
            raise SelectionError(
                f"{task_id}.variant_primary_capability must be null"
            )
        return
    binding_fields = frozenset(
        {
            "schema_version",
            "task_id",
            "primary_endpoint_id",
            "primary_evaluator_id",
            "capability_registry_sha256",
            "selected",
            "baseline",
        }
    )
    binding = _strict_mapping(
        raw_binding,
        required=binding_fields,
        allowed=binding_fields,
        label=f"{task_id}.variant_primary_capability",
    )
    expected_identity = {
        "schema_version": _VARIANT_CAPABILITY_BINDING_SCHEMA_VERSION,
        "task_id": _VARIANT_TASK_ID,
        "primary_endpoint_id": _VARIANT_PRIMARY_ENDPOINT_ID,
        "primary_evaluator_id": _VARIANT_PRIMARY_EVALUATOR_ID,
    }
    for field_name, expected in expected_identity.items():
        if binding[field_name] != expected:
            raise SelectionError(
                f"{task_id}.variant_primary_capability.{field_name} is not canonical"
            )
    _sha256(
        binding["capability_registry_sha256"],
        f"{task_id}.variant_primary_capability.capability_registry_sha256",
    )
    role_fields = frozenset({"model_id", "capability", "capability_sha256"})
    for role, model_field in (
        ("selected", "selected_model_id"),
        ("baseline", "baseline_model_id"),
    ):
        role_binding = _strict_mapping(
            binding[role],
            required=role_fields,
            allowed=role_fields,
            label=f"{task_id}.variant_primary_capability.{role}",
        )
        model_id = _nonempty_string(
            role_binding["model_id"],
            f"{task_id}.variant_primary_capability.{role}.model_id",
        )
        if model_id != decision[model_field]:
            raise SelectionError(
                f"{task_id} locked {role} capability names the wrong model"
            )
        capability = _validated_primary_variant_capability_record(
            role_binding["capability"],
            label=f"{task_id}.variant_primary_capability.{role}.capability",
            model_id=model_id,
        )
        capability_sha256 = _sha256(
            role_binding["capability_sha256"],
            f"{task_id}.variant_primary_capability.{role}.capability_sha256",
        )
        if capability_sha256 != canonical_hash(capability):
            raise SelectionError(
                f"{task_id} locked {role} capability SHA-256 changed"
            )


def _validate_locked_variant_secondary_evaluation(
    *, task_id: str, decision: Mapping[str, Any]
) -> None:
    raw_binding = decision.get("variant_secondary_evaluation")
    if task_id != _VARIANT_TASK_ID:
        if raw_binding is not None:
            raise SelectionError(
                f"{task_id}.variant_secondary_evaluation must be null"
            )
        return
    binding_fields = frozenset(
        {
            "schema_version",
            "task_id",
            "endpoint_id",
            "evaluator_id",
            "capability_registry_sha256",
            "candidate_transform_id",
            "primary_baseline_transform_id",
            "mandatory_model_ids",
            "comparators",
        }
    )
    binding = _strict_mapping(
        raw_binding,
        required=binding_fields,
        allowed=binding_fields,
        label=f"{task_id}.variant_secondary_evaluation",
    )
    expected_identity = {
        "schema_version": _VARIANT_SECONDARY_BINDING_SCHEMA_VERSION,
        "task_id": _VARIANT_TASK_ID,
        "endpoint_id": _VARIANT_SECONDARY_ENDPOINT_ID,
        "evaluator_id": _VARIANT_SECONDARY_EVALUATOR_ID,
        "candidate_transform_id": _VARIANT_PRIMARY_RETRIEVAL_TRANSFORM_ID,
        "primary_baseline_transform_id": _VARIANT_PRIMARY_RETRIEVAL_TRANSFORM_ID,
    }
    for field_name, expected in expected_identity.items():
        if binding[field_name] != expected:
            raise SelectionError(
                f"{task_id}.variant_secondary_evaluation.{field_name} is not canonical"
            )
    registry_sha256 = _sha256(
        binding["capability_registry_sha256"],
        f"{task_id}.variant_secondary_evaluation.capability_registry_sha256",
    )
    primary = decision.get("variant_primary_capability")
    if (
        not isinstance(primary, Mapping)
        or primary.get("capability_registry_sha256") != registry_sha256
    ):
        raise SelectionError(
            "locked primary and secondary variant capability registries differ"
        )
    # SelectionLock normalizes every nested array in a task decision to a
    # tuple (contracts._as_metadata), so comparing against list(...) here could
    # never succeed once the decision had been through the requirements.  Compare
    # type-agnostically against the frozen roster.
    if tuple(binding["mandatory_model_ids"]) != _VARIANT_MANDATORY_SECONDARY_MODELS:
        raise SelectionError(
            f"{task_id}.variant_secondary_evaluation mandatory roster changed"
        )
    raw_comparators = binding["comparators"]
    if not isinstance(raw_comparators, (list, tuple)):
        raise SelectionError(
            f"{task_id}.variant_secondary_evaluation.comparators must be an array"
        )
    comparator_fields = frozenset(
        {
            "model_id",
            "run_id",
            "capability",
            "capability_sha256",
            "score_transform_id",
        }
    )
    observed_models: list[str] = []
    observed_runs: list[str] = []
    for index, value in enumerate(raw_comparators):
        comparator = _strict_mapping(
            value,
            required=comparator_fields,
            allowed=comparator_fields,
            label=f"{task_id}.variant_secondary_evaluation.comparators[{index}]",
        )
        model_id = _nonempty_string(
            comparator["model_id"],
            f"{task_id}.variant_secondary_evaluation.comparators[{index}].model_id",
        )
        observed_models.append(model_id)
        observed_runs.append(
            _sha256(
                comparator["run_id"],
                f"{task_id}.variant_secondary_evaluation.comparators[{index}].run_id",
            )
        )
        if comparator["score_transform_id"] != _VARIANT_SECONDARY_SCORE_TRANSFORM_ID:
            raise SelectionError(
                f"secondary comparator {model_id} score transform changed"
            )
        capability = _validated_secondary_variant_capability_record(
            comparator["capability"],
            label=f"locked secondary comparator {model_id} capability",
            model_id=model_id,
        )
        capability_sha256 = _sha256(
            comparator["capability_sha256"],
            f"locked secondary comparator {model_id} capability SHA-256",
        )
        if capability_sha256 != canonical_hash(capability):
            raise SelectionError(
                f"locked secondary comparator {model_id} capability SHA-256 changed"
            )
    if tuple(observed_models) != _VARIANT_MANDATORY_SECONDARY_MODELS:
        raise SelectionError(
            f"{task_id}.variant_secondary_evaluation comparator roster changed"
        )
    if len(set(observed_runs)) != len(observed_runs):
        raise SelectionError(
            f"{task_id}.variant_secondary_evaluation repeats a comparator run_id"
        )
    primary_runs = {
        str(record["run_id"])
        for role in ("selected_runs", "baseline_runs")
        for record in decision[role]
    }
    if primary_runs.intersection(observed_runs):
        raise SelectionError(
            "variant secondary comparator runs overlap signed primary runs"
        )


def _validate_locked_task_decision(raw: Mapping[str, Any]) -> dict[str, Any]:
    decision = _strict_mapping(
        raw,
        required=_LOCKED_TASK_DECISION_FIELDS,
        allowed=_LOCKED_TASK_DECISION_FIELDS,
        label="SelectionLock.task_decision",
    )
    task_id = _nonempty_string(decision["task_id"], "task_decision.task_id")
    for field_name in (
        "selected_model_id",
        "baseline_model_id",
        "selected_adaptation_regime",
        "baseline_adaptation_regime",
        "promotion_gate",
    ):
        _nonempty_string(decision[field_name], f"{task_id}.{field_name}")
    if decision["selected_model_id"] in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
        raise SelectionError(
            f"negative control {decision['selected_model_id']} cannot be selected"
        )
    if decision["selected_model_id"] == decision["baseline_model_id"]:
        raise SelectionError(
            f"{task_id} selected and baseline model IDs must be distinct under the "
            f"{BASELINE_COMPARATOR_POLICY} policy"
        )
    for field_name in _TASK_HASH_FIELDS:
        _sha256(decision[field_name], f"{task_id}.{field_name}")
    seeds = _five_seeds(decision["seeds"], task_id)
    if list(seeds) != list(decision["seeds"]):
        raise SelectionError(f"{task_id}.seeds must be in canonical sorted order")
    if not isinstance(decision["open_champion"], bool):
        raise SelectionError(f"{task_id}.open_champion must be boolean")
    if (
        decision["open_champion"]
        and decision["selected_model_id"] in OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS
    ):
        raise SelectionError(
            f"control {decision['selected_model_id']} cannot receive an open champion claim"
        )
    _json_mapping(decision["thresholds"], f"{task_id}.thresholds")

    for role in ("selected", "baseline"):
        records = decision[f"{role}_runs"]
        if not isinstance(records, (list, tuple)) or len(records) != len(seeds):
            raise SelectionError(
                f"{task_id}.{role}_runs must bind every frozen seed exactly once"
            )
        observed: list[dict[str, Any]] = []
        for index, record in enumerate(records):
            parsed = _strict_mapping(
                record,
                required=frozenset({"seed", "run_id"}),
                allowed=frozenset({"seed", "run_id"}),
                label=f"{task_id}.{role}_runs[{index}]",
            )
            seed = parsed["seed"]
            if isinstance(seed, bool) or not isinstance(seed, int):
                raise SelectionError(
                    f"{task_id}.{role}_runs[{index}].seed must be an integer"
                )
            observed.append(
                {
                    "seed": seed,
                    "run_id": _sha256(
                        parsed["run_id"], f"{task_id}.{role}_runs[{index}].run_id"
                    ),
                }
            )
        if [record["seed"] for record in observed] != list(seeds):
            raise SelectionError(
                f"{task_id}.{role}_runs are not aligned to the locked seeds"
            )
        if len({record["run_id"] for record in observed}) != len(observed):
            raise SelectionError(f"{task_id}.{role}_runs repeat a run_id")

    candidate_dataset_ids = decision["candidate_dataset_ids"]
    sealed_dataset_ids = decision["sealed_dataset_ids"]
    if not isinstance(candidate_dataset_ids, (list, tuple)):
        raise SelectionError(f"{task_id}.candidate_dataset_ids must be an array")
    if not isinstance(sealed_dataset_ids, (list, tuple)):
        raise SelectionError(f"{task_id}.sealed_dataset_ids must be an array")
    for field_name, values in (
        ("candidate_dataset_ids", candidate_dataset_ids),
        ("sealed_dataset_ids", sealed_dataset_ids),
    ):
        if any(not isinstance(item, str) or not item for item in values):
            raise SelectionError(f"{task_id}.{field_name} must contain identifiers")
        if len(values) != len(set(values)) or list(values) != sorted(values):
            raise SelectionError(
                f"{task_id}.{field_name} must be unique and canonically sorted"
            )
    if not set(sealed_dataset_ids).issubset(candidate_dataset_ids):
        raise SelectionError(
            f"{task_id}.sealed_dataset_ids must be a subset of candidate_dataset_ids"
        )
    bindings = decision["dataset_registry_sha256s"]
    if not isinstance(bindings, Mapping) or set(bindings) != set(candidate_dataset_ids):
        raise SelectionError(
            f"{task_id}.dataset_registry_sha256s must cover candidate datasets exactly"
        )
    for dataset_id, digest in bindings.items():
        _sha256(digest, f"{task_id}.dataset_registry_sha256s[{dataset_id!r}]")
    _validate_locked_variant_primary_capability(task_id=task_id, decision=decision)
    _validate_locked_variant_secondary_evaluation(
        task_id=task_id, decision=decision
    )
    return decision


def _validate_lock_structure(lock: SelectionLock) -> None:
    if lock.schema_version != SELECTION_SCHEMA_VERSION:
        raise SelectionError(
            f"unsupported selection-lock schema_version: {lock.schema_version}"
        )
    if lock.terminal_policy != TERMINAL_POLICY:
        raise SelectionError("selection lock terminal policy is not fail-closed")
    if not lock.locked or lock.outcomes_unlocked:
        raise SelectionError("selection lock must precede every sealed outcome unlock")
    if lock.release_state is not ReleaseState.CANDIDATE:
        raise SelectionError("new selection locks must have release_state=candidate")
    if list(lock.candidate_run_ids) != sorted(lock.candidate_run_ids):
        raise SelectionError("candidate_run_ids must be canonically sorted")
    if list(lock.selected_run_ids) != sorted(lock.selected_run_ids):
        raise SelectionError("selected_run_ids must be canonically sorted")
    if not set(lock.selected_run_ids).issubset(lock.candidate_run_ids):
        raise SelectionError("selected task runs must be candidate_run_ids")
    decisions = [
        _validate_locked_task_decision(decision) for decision in lock.task_decisions
    ]
    task_ids = [str(decision["task_id"]) for decision in decisions]
    if task_ids != sorted(task_ids):
        raise SelectionError("task_decisions must be canonically sorted by task_id")
    selected_ids = sorted(
        str(record["run_id"])
        for decision in decisions
        for record in decision["selected_runs"]
    )
    if selected_ids != list(lock.selected_run_ids):
        raise SelectionError(
            "aggregate selected_run_ids do not equal task-decision selected runs"
        )
    baseline_ids = {
        str(record["run_id"])
        for decision in decisions
        for record in decision["baseline_runs"]
    }
    if not baseline_ids.issubset(lock.candidate_run_ids):
        raise SelectionError("baseline task runs must be candidate_run_ids")
    secondary_ids = {
        str(record["run_id"])
        for decision in decisions
        if isinstance(decision.get("variant_secondary_evaluation"), Mapping)
        for record in decision["variant_secondary_evaluation"]["comparators"]
    }
    if not secondary_ids.issubset(lock.candidate_run_ids):
        raise SelectionError("secondary comparator runs must be candidate_run_ids")


def _canonical_metadata(value: object, label: str) -> object:
    """Compare requirements-normalized metadata against freshly built structures.

    ``contracts._as_metadata`` freezes every nested array into a tuple, so a
    plain ``==`` between selection record metadata and a rebuilt list is always False no
    matter what it holds.  ``canonicalize`` renders tuples and lists alike as
    JSON arrays, so this compares structure and values, not container type.
    """

    try:
        return canonicalize(value)
    except (HashingError, TypeError, ValueError) as error:
        raise SelectionError(f"{label} is not canonical JSON: {error}") from error


def _validate_lock_against_candidate(
    lock: SelectionLock,
    *,
    candidate_path: Path,
    plan: Mapping[str, Any],
    candidate_manifest_sha256: str,
    artifact_class: str,
) -> None:
    """Recursively bind fixed model labels to the exact candidate run records."""

    expected_binding = {
        "path": candidate_path.as_posix(),
        "manifest_sha256": candidate_manifest_sha256,
        "artifact_class": artifact_class,
        "source_identity_sha256": str(plan["plan_sha256"]),
    }
    metadata = lock.metadata
    if not isinstance(metadata, Mapping):
        raise SelectionError("selection lock metadata must be an object")
    if _canonical_metadata(
        metadata.get("candidate_binding"), "selection lock candidate_binding"
    ) != _canonical_metadata(expected_binding, "expected candidate_binding"):
        raise SelectionError("selection lock candidate binding changed")
    expected_universe = (
        dict(plan["finalist_campaign_universe_binding"])
        if artifact_class == "selection_candidate_ledger"
        else {}
    )
    if _canonical_metadata(
        metadata.get("finalist_campaign_universe_binding"),
        "selection lock finalist_campaign_universe_binding",
    ) != _canonical_metadata(
        expected_universe, "expected finalist_campaign_universe_binding"
    ):
        raise SelectionError(
            "selection lock finalist campaign universe binding changed"
        )
    if metadata.get("fit_state_artifact_policy") != FIT_STATE_ARTIFACT_POLICY:
        raise SelectionError("selection lock fit-state artifact policy changed")
    if (
        lock.campaign_id != str(plan["campaign"]["campaign_id"])
        or lock.plan_sha256 != str(plan["plan_sha256"])
        or lock.candidate_manifest_sha256 != candidate_manifest_sha256
        or lock.registry_sha256 != str(plan["registry_snapshot_sha256"])
    ):
        raise SelectionError("selection lock identities differ from its candidate")
    if (
        artifact_class == "selection_candidate_ledger"
        and lock.metrics_sha256 != str(plan["plan_sha256"])
    ):
        raise SelectionError(
            "ledger-backed selection metrics do not identify the verified ledger"
        )

    raw_runs = plan.get("runs")
    if not isinstance(raw_runs, list):
        raise SelectionError("selection candidate runs must be an array")
    runs_by_id: dict[str, Mapping[str, Any]] = {}
    for index, raw_run in enumerate(raw_runs):
        if not isinstance(raw_run, Mapping):
            raise SelectionError(f"candidate runs[{index}] must be an object")
        run_id = _sha256(raw_run.get("run_id"), f"candidate runs[{index}].run_id")
        if run_id in runs_by_id:
            raise SelectionError("selection candidate repeats run IDs")
        runs_by_id[run_id] = raw_run
    if set(lock.candidate_run_ids) != set(runs_by_id):
        raise SelectionError("selection lock candidate_run_ids changed from its candidate")

    run_index_lists: dict[
        tuple[str, str, int, str], list[Mapping[str, Any]]
    ] = {}
    for run in runs_by_id.values():
        key = (
            str(run["task_id"]),
            str(run["model_id"]),
            int(run["seed"]),
            str(run["adaptation_regime"]),
        )
        run_index_lists.setdefault(key, []).append(run)
    run_index = {key: tuple(value) for key, value in run_index_lists.items()}
    known_tasks = {str(run["task_id"]) for run in runs_by_id.values()}

    dispositions = _model_dispositions(plan)
    expected_model_receipts: list[dict[str, Any]] = []
    expected_receipts: list[dict[str, Any]] = []
    observed_selected_run_ids: set[str] = set()
    for decision in lock.task_decisions:
        task_id = str(decision["task_id"])
        selected_model_id = str(decision["selected_model_id"])
        baseline_model_id = str(decision["baseline_model_id"])
        baseline_roster = _task_baseline_roster(plan, task_id)
        if selected_model_id in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
            raise SelectionError(
                f"negative control {selected_model_id} cannot be selected"
            )
        if baseline_model_id in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
            raise SelectionError(
                f"negative control {baseline_model_id} cannot be a baseline comparator"
            )
        if baseline_model_id not in baseline_roster:
            raise SelectionError(
                f"{task_id} locked baseline is absent from its TaskSpec roster"
            )
        if selected_model_id == baseline_model_id:
            raise SelectionError(
                f"{task_id} selected and baseline model IDs must be distinct"
            )
        if decision["open_champion"]:
            disposition = dispositions.get(selected_model_id)
            if (
                artifact_class != "selection_candidate_ledger"
                or disposition is None
                or disposition.get("champion_eligible") is not True
                or selected_model_id in OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS
            ):
                raise SelectionError(
                    f"{task_id} open champion lacks an eligible scientific ledger"
                )

        role_runs: dict[str, tuple[Mapping[str, Any], ...]] = {}
        for role, model_field, regime_field in (
            ("selected", "selected_model_id", "selected_adaptation_regime"),
            ("baseline", "baseline_model_id", "baseline_adaptation_regime"),
        ):
            resolved: list[Mapping[str, Any]] = []
            for locked_run in decision[f"{role}_runs"]:
                run_id = str(locked_run["run_id"])
                run = runs_by_id.get(run_id)
                if run is None:
                    raise SelectionError(
                        f"{task_id} {role} run is absent from its candidate: {run_id}"
                    )
                if (
                    run.get("task_id") != task_id
                    or run.get("model_id") != decision[model_field]
                    or run.get("adaptation_regime") != decision[regime_field]
                    or run.get("seed") != locked_run["seed"]
                ):
                    raise SelectionError(
                        f"{task_id} {role} run/model/regime/seed binding changed"
                    )
                resolved.append(run)
                if role == "selected":
                    observed_selected_run_ids.add(run_id)
            role_runs[role] = tuple(resolved)
        if set(str(run["run_id"]) for run in role_runs["selected"]) & set(
            str(run["run_id"]) for run in role_runs["baseline"]
        ):
            raise SelectionError(
                f"{task_id} selected and baseline run sets must be disjoint"
            )

        raw_decision = {
            field_name: decision[field_name]
            for field_name in _TASK_DECISION_REQUIRED
        }
        raw_decision["open_champion"] = decision["open_champion"]
        if task_id == _VARIANT_TASK_ID:
            secondary = decision["variant_secondary_evaluation"]
            if not isinstance(secondary, Mapping):
                raise SelectionError(
                    "variant lock lacks its secondary evaluation binding"
                )
            raw_decision["variant_secondary_comparators"] = [
                {
                    "model_id": comparator["model_id"],
                    "run_id": comparator["run_id"],
                }
                for comparator in secondary["comparators"]
            ]
        expected_decision = _normalize_task_decision(
            raw_decision,
            plan=plan,
            run_index=run_index,
            runs_by_id=runs_by_id,
            known_tasks=known_tasks,
            dispositions=dispositions,
        )
        try:
            locked_decision = canonicalize(decision)
        except (HashingError, TypeError, ValueError) as error:
            raise SelectionError(
                f"{task_id} locked task decision is not canonical JSON: {error}"
            ) from error
        if expected_decision != locked_decision:
            raise SelectionError(
                f"{task_id} task decision no longer rederives from its candidate"
            )

        task_records = [
            item
            for item in plan.get("task_dispositions", ())
            if isinstance(item, Mapping) and item.get("task_id") == task_id
        ]
        if (
            artifact_class == "selection_candidate_ledger"
            and (
                len(task_records) != 1
                or task_records[0].get("promotion_gate_id")
                != decision["promotion_gate"]
            )
        ):
            raise SelectionError(
                f"{task_id}.promotion_gate differs from its frozen TaskSpec"
            )

        selected_receipt = _ranked_selected_receipt(
            plan=plan,
            task_id=task_id,
            selected_model_id=selected_model_id,
            selected_regime=str(decision["selected_adaptation_regime"]),
            seeds=tuple(int(seed) for seed in decision["seeds"]),
            selected_runs=role_runs["selected"],
            open_champion=bool(decision["open_champion"]),
        )
        if selected_receipt is not None:
            expected_model_receipts.append(selected_receipt)
        receipt = _ranked_baseline_receipt(
            plan=plan,
            task_id=task_id,
            selected_model_id=selected_model_id,
            baseline_model_id=baseline_model_id,
            baseline_regime=str(decision["baseline_adaptation_regime"]),
            baseline_roster=baseline_roster,
            seeds=tuple(int(seed) for seed in decision["seeds"]),
            baseline_runs=role_runs["baseline"],
        )
        if receipt is not None:
            expected_receipts.append(receipt)

    if observed_selected_run_ids != set(lock.selected_run_ids):
        raise SelectionError("selection lock selected_run_ids changed from task decisions")
    expected_policy = (
        BASELINE_COMPARATOR_POLICY
        if artifact_class == "selection_candidate_ledger"
        else UNRANKED_BASELINE_COMPARATOR_POLICY
    )
    if metadata.get("baseline_comparator_policy") != expected_policy:
        raise SelectionError("selection lock baseline comparator policy changed")
    expected_model_policy = (
        MODEL_SELECTION_POLICY
        if artifact_class == "selection_candidate_ledger"
        else "explicit_nonchampion_contract_selection_v1"
    )
    if metadata.get("model_selection_policy") != expected_model_policy:
        raise SelectionError("selection lock model-selection policy changed")
    expected_model_receipts.sort(key=lambda item: item["task_id"])
    # SelectionLock normalizes every nested array in metadata to a tuple
    # (contracts._as_metadata), so a direct == against freshly built lists can
    # never succeed.  canonicalize() renders both sides in one canonical JSON
    # shape, which compares structure and values rather than container type.
    if _canonical_metadata(
        metadata.get("model_selection_receipts"),
        "selection lock model_selection_receipts",
    ) != _canonical_metadata(
        expected_model_receipts, "expected model_selection_receipts"
    ):
        raise SelectionError("selection lock model-selection receipts changed")
    expected_receipts.sort(key=lambda item: item["task_id"])
    if _canonical_metadata(
        metadata.get("baseline_comparator_receipts"),
        "selection lock baseline_comparator_receipts",
    ) != _canonical_metadata(
        expected_receipts, "expected baseline_comparator_receipts"
    ):
        raise SelectionError("selection lock baseline ranking receipts changed")


def build_selection_lock(
    *, candidate: str | Path, decisions: Mapping[str, Any]
) -> SelectionLock:
    """Build a strict aggregate record from a frozen campaign and development metrics."""

    decision_document = _strict_mapping(
        decisions,
        required=_DECISION_DOCUMENT_REQUIRED,
        allowed=_DECISION_DOCUMENT_ALLOWED,
        label="selection decision document",
    )
    sealed_results_used = decision_document.get("sealed_results_used", False)
    if not isinstance(sealed_results_used, bool):
        raise SelectionError("sealed_results_used must be boolean")
    if sealed_results_used:
        raise SelectionError("selection decisions may not use project-sealed results")

    try:
        candidate_path = reject_symlink_components(
            Path(candidate), label="selection candidate"
        ).resolve()
    except (ArtifactError, OSError) as error:
        raise SelectionError(f"invalid selection candidate: {error}") from error
    plan, candidate_manifest_sha256, artifact_class = _validated_plan(
        candidate_path
    )
    ledger = plan.get("selection_candidate_ledger")
    if isinstance(ledger, Mapping) and decision_document["metrics_sha256"] != ledger.get(
        "ledger_id"
    ):
        raise SelectionError(
            "ledger-backed selection metrics_sha256 must equal the verified ledger_id"
        )
    runs = tuple(plan["runs"])
    known_tasks = {str(run["task_id"]) for run in runs}
    run_index_lists: dict[
        tuple[str, str, int, str], list[Mapping[str, Any]]
    ] = {}
    for run in runs:
        key = (
            str(run["task_id"]),
            str(run["model_id"]),
            int(run["seed"]),
            str(run["adaptation_regime"]),
        )
        run_index_lists.setdefault(key, []).append(run)
    run_index = {key: tuple(value) for key, value in run_index_lists.items()}
    runs_by_id = {str(run["run_id"]): run for run in runs}
    dispositions = _model_dispositions(plan)

    raw_task_decisions = decision_document["task_decisions"]
    if not isinstance(raw_task_decisions, list) or not raw_task_decisions:
        raise SelectionError("at least one task_decision is required")
    normalized: list[dict[str, Any]] = []
    seen_tasks: set[str] = set()
    for raw in raw_task_decisions:
        if not isinstance(raw, Mapping):
            raise SelectionError("each task_decision must be an object")
        decision = _normalize_task_decision(
            raw,
            plan=plan,
            run_index=run_index,
            runs_by_id=runs_by_id,
            known_tasks=known_tasks,
            dispositions=dispositions,
        )
        task_id = str(decision["task_id"])
        if task_id in seen_tasks:
            raise SelectionError(f"duplicate task decision: {task_id}")
        seen_tasks.add(task_id)
        normalized.append(decision)
    missing_tasks = sorted(known_tasks.difference(seen_tasks))
    extra_tasks = sorted(seen_tasks.difference(known_tasks))
    if missing_tasks or extra_tasks:
        details: list[str] = []
        if missing_tasks:
            details.append("missing=" + ",".join(missing_tasks))
        if extra_tasks:
            details.append("extra=" + ",".join(extra_tasks))
        raise SelectionError(
            "selection lock must cover every and only schedulable task: "
            + "; ".join(details)
        )
    normalized.sort(key=lambda item: str(item["task_id"]))

    conditional = _json_mapping(
        decision_document.get(
            "conditional_model_decision", {"status": "not_evaluated"}
        ),
        "conditional_model_decision",
    )
    power_decisions = _json_mapping(
        decision_document.get("power_decisions", {}), "power_decisions"
    )
    multiplicity_plan = _nonempty_string(
        decision_document.get("multiplicity_plan", DEFAULT_MULTIPLICITY_PLAN),
        "multiplicity_plan",
    )
    decision_metadata = _json_mapping(
        decision_document.get("metadata", {}), "metadata"
    )
    model_selection_receipts: list[dict[str, Any]] = []
    baseline_comparator_receipts: list[dict[str, Any]] = []
    if artifact_class == "selection_candidate_ledger":
        for decision in normalized:
            task_id = str(decision["task_id"])
            selected_receipt = _ranked_selected_receipt(
                plan=plan,
                task_id=task_id,
                selected_model_id=str(decision["selected_model_id"]),
                selected_regime=str(decision["selected_adaptation_regime"]),
                seeds=tuple(int(seed) for seed in decision["seeds"]),
                selected_runs=tuple(
                    runs_by_id[str(record["run_id"])]
                    for record in decision["selected_runs"]
                ),
                open_champion=bool(decision["open_champion"]),
            )
            if selected_receipt is None:
                raise SelectionError(
                    f"{task_id} lacks a ranked model-selection receipt"
                )
            model_selection_receipts.append(selected_receipt)
            receipt = _ranked_baseline_receipt(
                plan=plan,
                task_id=task_id,
                selected_model_id=str(decision["selected_model_id"]),
                baseline_model_id=str(decision["baseline_model_id"]),
                baseline_regime=str(decision["baseline_adaptation_regime"]),
                baseline_roster=_task_baseline_roster(plan, task_id),
                seeds=tuple(int(seed) for seed in decision["seeds"]),
                baseline_runs=tuple(
                    runs_by_id[str(record["run_id"])]
                    for record in decision["baseline_runs"]
                ),
            )
            if receipt is None:
                raise SelectionError(
                    f"{task_id} lacks a ranked baseline receipt"
                )
            baseline_comparator_receipts.append(receipt)
        model_selection_receipts.sort(key=lambda item: item["task_id"])
        baseline_comparator_receipts.sort(key=lambda item: item["task_id"])
    try:
        metadata = {
            "source_lock_sha256": canonical_hash(plan.get("source_lock")),
            "resource_firewall_sha256": canonical_hash(plan.get("resource_firewall")),
            "candidate_binding": {
                "path": candidate_path.as_posix(),
                "manifest_sha256": candidate_manifest_sha256,
                "artifact_class": artifact_class,
                "source_identity_sha256": str(plan["plan_sha256"]),
            },
            "baseline_comparator_policy": (
                BASELINE_COMPARATOR_POLICY
                if artifact_class == "selection_candidate_ledger"
                else UNRANKED_BASELINE_COMPARATOR_POLICY
            ),
            "model_selection_policy": (
                MODEL_SELECTION_POLICY
                if artifact_class == "selection_candidate_ledger"
                else "explicit_nonchampion_contract_selection_v1"
            ),
            "model_selection_receipts": model_selection_receipts,
            "baseline_comparator_receipts": baseline_comparator_receipts,
            "decision_metadata": decision_metadata,
            "fit_state_artifact_policy": FIT_STATE_ARTIFACT_POLICY,
            "finalist_campaign_universe_binding": (
                dict(plan["finalist_campaign_universe_binding"])
                if artifact_class == "selection_candidate_ledger"
                else {}
            ),
        }
    except (HashingError, TypeError, ValueError) as error:
        raise SelectionError(f"candidate source locks are not canonical JSON: {error}") from error
    candidate_run_ids = sorted(str(run["run_id"]) for run in runs)
    selected_run_ids = sorted(
        str(record["run_id"])
        for decision in normalized
        for record in decision["selected_runs"]
    )
    identity = {
        "schema_version": SELECTION_SCHEMA_VERSION,
        "campaign_id": str(plan["campaign"]["campaign_id"]),
        "plan_sha256": str(plan["plan_sha256"]),
        "candidate_manifest_sha256": candidate_manifest_sha256,
        "registry_sha256": str(plan["registry_snapshot_sha256"]),
        "metrics_sha256": _sha256(
            decision_document["metrics_sha256"], "metrics_sha256"
        ),
        "locked": True,
        "outcomes_unlocked": False,
        "release_state": ReleaseState.CANDIDATE.value,
        "candidate_run_ids": candidate_run_ids,
        "selected_run_ids": selected_run_ids,
        "task_decisions": normalized,
        "conditional_model_decision": conditional,
        "power_decisions": power_decisions,
        "multiplicity_plan": multiplicity_plan,
        "terminal_policy": TERMINAL_POLICY,
        "metadata": metadata,
    }
    payload = {"lock_id": canonical_hash(identity), **identity}
    try:
        lock = SelectionLock.from_dict(payload)
    except ContractError as error:
        raise SelectionError(f"invalid aggregate SelectionLock: {error}") from error
    _validate_lock_structure(lock)
    _validate_lock_against_candidate(
        lock,
        candidate_path=candidate_path,
        plan=plan,
        candidate_manifest_sha256=candidate_manifest_sha256,
        artifact_class=artifact_class,
    )
    return lock


def freeze_selection_lock(
    *,
    candidate: str | Path,
    decisions_path: str | Path,
    output_root: str | Path,
) -> Path:
    """Write a SelectionLock as a read-only, independently verifiable directory."""

    try:
        decisions_source = reject_symlink_components(
            Path(decisions_path), label="selection decisions"
        )
    except (ArtifactError, OSError) as error:
        raise SelectionError(f"invalid selection decisions path: {error}") from error
    decisions = _load_json(decisions_source)
    lock = build_selection_lock(candidate=candidate, decisions=decisions)
    try:
        root = reject_symlink_components(
            Path(output_root), label="selection output root"
        ).resolve()
    except (ArtifactError, OSError) as error:
        raise SelectionError(f"invalid selection output root: {error}") from error
    root.mkdir(parents=True, exist_ok=True)
    target = root / f"selection--{lock.lock_id}"
    if target.exists():
        raise SelectionError(f"selection lock already exists and is immutable: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".selection.", dir=root))
    try:
        write_json_exclusive(
            staging / "selection_lock.json", lock.to_dict(), mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "selection_lock",
                "lock_id": lock.lock_id,
            },
        )
        try:
            publish_directory_noreplace(staging, target)
        except (ArtifactError, OSError) as error:
            raise SelectionError(
                f"could not publish immutable selection lock {target}: {error}"
            ) from error
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verified = verify_selection_lock(target)
    if verified.lock_id != lock.lock_id:
        raise SelectionError("published selection lock failed identity verification")
    return target


def verify_selection_lock(path: str | Path) -> SelectionLock:
    """Verify the frozen tree, strict schema, and canonical selection record identity."""

    try:
        root = reject_symlink_components(
            Path(path), label="selection lock"
        ).resolve()
    except (ArtifactError, OSError) as error:
        raise SelectionError(f"invalid selection lock path: {error}") from error
    if not root.is_dir():
        raise SelectionError(f"selection lock must be a frozen directory: {root}")
    try:
        manifest = verify_frozen_tree(root)
        lock = SelectionLock.load_json(root / "selection_lock.json")
    except (ArtifactError, ContractError, OSError, ValueError, KeyError, TypeError) as error:
        raise SelectionError(f"invalid frozen selection lock: {error}") from error
    metadata = manifest.get("metadata")
    if not isinstance(metadata, Mapping):
        raise SelectionError("selection-lock artifact manifest metadata is missing")
    if metadata.get("artifact_class") != "selection_lock":
        raise SelectionError("artifact directory is not classified as a selection lock")
    if metadata.get("lock_id") != lock.lock_id:
        raise SelectionError("artifact manifest and SelectionLock disagree on lock_id")
    identity = lock.to_dict()
    claimed = identity.pop("lock_id")
    if canonical_hash(identity) != claimed:
        raise SelectionError("selection lock identity hash mismatch")
    _validate_lock_structure(lock)
    candidate_binding = lock.metadata.get("candidate_binding")
    if not isinstance(candidate_binding, Mapping):
        raise SelectionError("selection lock lacks a candidate binding")
    try:
        candidate_path = reject_symlink_components(
            Path(str(candidate_binding.get("path", ""))),
            label="selection lock candidate binding",
        ).resolve()
        plan, candidate_manifest_sha256, artifact_class = _validated_plan(
            candidate_path
        )
    except (ArtifactError, OSError, ValueError, TypeError) as error:
        raise SelectionError(
            f"selection lock candidate binding is invalid: {error}"
        ) from error
    _validate_lock_against_candidate(
        lock,
        candidate_path=candidate_path,
        plan=plan,
        candidate_manifest_sha256=candidate_manifest_sha256,
        artifact_class=artifact_class,
    )
    try:
        lock.require_locked(for_external_scoring=True)
    except ContractError as error:
        raise SelectionError(str(error)) from error
    return lock


__all__ = [
    "BASELINE_COMPARATOR_POLICY",
    "CANDIDATE_STANDARD_ERROR_POLICY",
    "CONDITIONAL_TRIGGER_INELIGIBLE_MODEL_IDS",
    "DEFAULT_MULTIPLICITY_PLAN",
    "DEVELOPMENT_ENSEMBLE_POLICY_ID",
    "FINALIST_LEDGER_AUTHORITY",
    "FINALIST_SEED_COUNT",
    "FIT_STATE_ARTIFACT_POLICY",
    "INDEPENDENT_ROLE_ARTIFACT_POLICY",
    "MODEL_SELECTION_POLICY",
    "NON_SELECTABLE_NEGATIVE_CONTROL_IDS",
    "OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS",
    "SEED_STABILITY_POLICY",
    "SELECTION_SCHEMA_VERSION",
    "SUPPORTED_FIT_STATE_ARTIFACT_POLICIES",
    "SelectionError",
    "TERMINAL_POLICY",
    "UNRANKED_BASELINE_COMPARATOR_POLICY",
    "build_selection_lock",
    "freeze_selection_lock",
    "verify_selection_lock",
]
