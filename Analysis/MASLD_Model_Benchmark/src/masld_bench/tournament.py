"""Deterministic finalist selection and fail-closed best-model adoption.

The aggregate SelectionLock freezes five seeded selected and baseline runs per
task.  A separate hashed terminal authorization binds those exact runs to one
held-back outcome bundle after predictions and prospective power are frozen.
Every opened bundle is terminal whether its check passes or fails.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
import csv
from hashlib import sha256
import json
from math import fsum, inf, isfinite, sqrt
from pathlib import Path
import re
import shutil
import tempfile
import tomllib
from typing import Any

from .artifacts import (
    ArtifactError,
    freeze_tree,
    publish_directory_noreplace,
    reject_symlink_components,
    sha256_file,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)
from .conditional_model import ComplementarityEvidence, ConditionalModelError
from .contracts import (
    ArtifactRef,
    ContractError,
    MissingState,
    PredictionBundle,
    ReleaseState,
    RunSpec,
    SelectionLock,
)
from .evaluators.endpoint_power import (
    EndpointBootstrapResult,
    EndpointPowerError,
    recompute_endpoint_bootstrap,
)
from .evaluators.sealed_metrics import (
    SealedMetricError,
    recompute_sealed_metrics,
)
from .evaluators.stats import holm_correction
from .firewall import (
    FirewallError,
    OutcomeConsumption,
    PowerDecision,
    PredictionCommit,
    verify_development_power_evidence,
    verify_outcome_consumption,
    verify_power_decision,
    verify_prediction_commit,
    verify_sealed_outcome_bundle,
)
from .hashing import HashingError, canonical_sha256, canonicalize
from .planner import PlanningError, load_frozen_plan
from .selection import (
    CANDIDATE_STANDARD_ERROR_POLICY,
    CONDITIONAL_TRIGGER_INELIGIBLE_MODEL_IDS,
    FIT_STATE_ARTIFACT_POLICY,
    NON_SELECTABLE_NEGATIVE_CONTROL_IDS,
    OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS,
    SUPPORTED_FIT_STATE_ARTIFACT_POLICIES,
    SelectionError,
    verify_selection_lock,
)


_NO_DEFAULT = object()
_MISSING = object()
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_PROMOTION_CONFIG_SHA256 = (
    "396729908f39ad5e6d32ff419ecaf73d110c568814cbe1c254bfc53518f5a25a"
)
_PROMOTION_CONFIG = (
    Path(__file__).resolve().parents[2]
    / "config"
    / "evaluation"
    / "promotion_gates.toml"
)

CELL_TASK = "cell_state_mapping"
VARIANT_TASK = "variant_to_regulation"
RNA_ATAC_TASK = "rna_conditioned_atac"
BULK_TASK = "bulk_state_transfer"
GRAPH_TASK = "typed_evidence_graph"
PERTURBATION_TASK = "perturbation_transfer"
UNIFIED_TASK = "unified_model"
CONDITIONAL_CONTEXT_TASK = "conditional_context"

_VARIANT_SECONDARY_MODELS = ("abc", "nearest_gene", "re2g")
_VARIANT_SECONDARY_ROLES = {
    "abc": "link_only",
    "nearest_gene": "baseline",
    "re2g": "link_only",
}
_VARIANT_PRIMARY_EVALUATOR_ID = "variant_ld_block_fisher_z_spearman_gain_v1"
_VARIANT_SECONDARY_ENDPOINT_ID = "eqtl_retrieval"
_VARIANT_SECONDARY_EVALUATOR_ID = "variant_ld_block_eqtl_retrieval_auprc_v1"
_VARIANT_SECONDARY_SCORE_TRANSFORM_ID = "identity_link_score_v1"

EXPECTED_TASK_IDS = frozenset(
    {
        CELL_TASK,
        VARIANT_TASK,
        RNA_ATAC_TASK,
        BULK_TASK,
        GRAPH_TASK,
        PERTURBATION_TASK,
        UNIFIED_TASK,
    }
)

SEALED_CONFIRMATORY = "sealed_confirmatory"
CONDITIONAL_SEALED = "conditional_sealed"
DEVELOPMENT_ONLY = "development_only"
EXPLORATORY_ONLY = "exploratory_only"
BLOCKED_MISSING_SEALED_SOURCE = "blocked_missing_sealed_source"
BLOCKED_MISSING_EXTERNAL_FAMILY = "blocked_missing_external_family"
PROMOTABLE_MODES = frozenset({SEALED_CONFIRMATORY, CONDITIONAL_SEALED})

SELECTION_LOCK_SCHEMA_VERSION = "masld-bench-selection-lock-v1"
DEVELOPMENT_SHORTLIST_SCHEMA_VERSION = "masld-bench-development-shortlist-v3"
# Which wave may produce a shortlist, the exact seeds that wave must have run,
# and what that shortlist is then allowed to authorize.  Mirrors
# config/campaigns/v1_frozen_screen.toml and v1_specialist_screen.toml; the
# contract tests assert the two agree.  A one-seed smoke ledger appears in
# neither entry and therefore cannot shortlist anything.
SHORTLIST_SOURCE_WAVE_CONTRACTS: Mapping[str, Mapping[str, Any]] = {
    "frozen_screen": {
        "expected_seeds": (1103, 2909, 4721),
        "authorizes_wave": "full_specialist_screen",
    },
    "full_specialist_screen": {
        "expected_seeds": (1103, 2909, 4721, 6673, 8111),
        "authorizes_wave": "adaptation",
    },
}
SCIENTIFIC_RUN_RECEIPT_SCHEMA_VERSION = "masld-bench-scientific-run-receipt-v3"
VARIANT_SECONDARY_RUN_RECEIPT_SCHEMA_VERSION = (
    "masld-bench-variant-secondary-run-receipt-v1"
)
DEVELOPMENT_OUTCOME_BUNDLE_SCHEMA_VERSION = (
    "masld-bench-development-outcome-bundle-v1"
)
SELECTION_CANDIDATE_LEDGER_SCHEMA_VERSION = (
    "masld-bench-selection-candidate-ledger-v5"
)
FINALIST_CAMPAIGN_UNIVERSE_SCHEMA_VERSION = (
    "masld-bench-finalist-campaign-universe-v2"
)
SCREENING_CAMPAIGN_UNIVERSE_SCHEMA_VERSION = (
    "masld-bench-screening-campaign-universe-v1"
)
# A ledger's authority decides what it may support.  A screening ledger may
# produce a shortlist and nothing else; only a finalist ledger bound to a
# pre-scoring campaign universe may back a SelectionLock or a best model.
LEDGER_AUTHORITY_SCREENING = "frozen_screen_shortlist_only_no_lock_v1"
LEDGER_AUTHORITY_FINALIST = (
    "finalist_campaign_universe_bound_champion_eligible_v1"
)
UNIVERSE_KIND_SCREENING = "screening_frozen_screen"
UNIVERSE_KIND_FINALIST = "finalist"
_SCREENING_UNIVERSE_WAVES = frozenset({"frozen_screen"})
_FINALIST_UNIVERSE_WAVES = frozenset(
    {"full_specialist_screen", "adaptation", "conditional_model"}
)
_SCREENING_SEED_COUNT = 3
_FINALIST_SEED_COUNT = 5
_CAMPAIGN_UNIVERSE_POLICIES: Mapping[str, tuple[frozenset[str], int, str]] = {
    UNIVERSE_KIND_SCREENING: (
        _SCREENING_UNIVERSE_WAVES,
        _SCREENING_SEED_COUNT,
        SCREENING_CAMPAIGN_UNIVERSE_SCHEMA_VERSION,
    ),
    UNIVERSE_KIND_FINALIST: (
        _FINALIST_UNIVERSE_WAVES,
        _FINALIST_SEED_COUNT,
        FINALIST_CAMPAIGN_UNIVERSE_SCHEMA_VERSION,
    ),
}
_LEDGER_AUTHORITY_POLICIES: Mapping[str, tuple[str, int]] = {
    LEDGER_AUTHORITY_SCREENING: (UNIVERSE_KIND_SCREENING, _SCREENING_SEED_COUNT),
    LEDGER_AUTHORITY_FINALIST: (UNIVERSE_KIND_FINALIST, _FINALIST_SEED_COUNT),
}
# Development ranking and held-back inference must estimate the SAME quantity.
# Held-back inference averages the five fixed seeds' predictions and evaluates the
# task-native endpoint once (evaluators/sealed_metrics.py).  Macro-F1, Fisher-z
# Spearman, average precision, and relative profile deviance are all nonlinear
# in the prediction, so mean(metric(seed)) != metric(mean(seed)).  Development
# therefore ranks the five-seed ensemble itself; per-seed metrics survive only
# as a stability diagnostic and never enter an ordering.
DEVELOPMENT_ENSEMBLE_POLICY_ID = (
    "mean_prediction_five_seed_development_ensemble_v1"
)
_DEVELOPMENT_ENSEMBLE_AGGREGATION = (
    "task_native_endpoint_on_mean_of_exactly_five_locked_seed_predictions"
)
_ENSEMBLE_ALIGNMENT_SCHEMA_VERSION = (
    "masld-bench-development-ensemble-alignment-v1"
)
_FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION = (
    "masld-bench-finalist-endpoint-distribution-v2"
)
_ENSEMBLE_SEED_COUNT = 5
_ENSEMBLE_PROBABILITY_TOLERANCE = 1e-6
# Averaging RAW predicted scores is not scale invariant: a seed that emits
# scores on a wider scale dominates the ensemble for every rank-based endpoint
# (Spearman, Fisher-z Spearman, average precision).  Held-back inference averages
# raw scores, so rank-averaging here would break parity; instead require the
# five seeds to be on a comparable scale and fail closed when they are not.
# The ratio is prospectively frozen and is not tuned on any outcome.
ENSEMBLE_SCORE_SCALE_POLICY = (
    "five_seed_raw_score_dispersion_ratio_at_most_4x_fail_closed_v1"
)
_ENSEMBLE_SCORE_SCALE_RATIO_MAX = 4.0
_CAMPAIGN_UNIVERSE_DOCUMENT = "campaign_universe.json"
_CAMPAIGN_UNIVERSE_ARTIFACT_CLASS = "campaign_universe"
_UNIVERSE_IDENTITY_KEYS = (
    "adaptation_regime",
    "candidate_configuration_sha256",
    "fold",
    "model_id",
    "run_id",
    "seed",
    "task_id",
)
FINALIST_DEVELOPMENT_METRIC_BUNDLE_SCHEMA_VERSION = (
    "masld-bench-finalist-development-metric-bundle-v1"
)
DEVELOPMENT_RESIDUAL_BUNDLE_SCHEMA_VERSION = (
    "masld-bench-development-residual-bundle-v1"
)
CONDITIONAL_DECISION_SCHEMA_VERSION = "masld-bench-conditional-decision-v4"
GATE_DECISION_SCHEMA_VERSION = "masld-bench-champion-gate-decision-v1"
SEALED_METRIC_BUNDLE_SCHEMA_VERSION = "masld-bench-sealed-metric-bundle-v2"
SEALED_EVALUATOR_SOURCE_BUNDLE_SCHEMA_VERSION = (
    "masld-bench-sealed-evaluator-source-bundle-v1"
)
CONFIRMATORY_MULTIPLICITY_BUNDLE_SCHEMA_VERSION = (
    "masld-bench-confirmatory-multiplicity-bundle-v1"
)
TERMINAL_AUTHORIZATION_SCHEMA_VERSION = "masld-bench-terminal-authorization-v1"
LOCKED_MULTIPLICITY_PLAN = (
    "confirmatory_holm_fwer_0.05_secondary_bh_by_family"
)
LOCKED_TERMINAL_POLICY = "no_reselection_recalibration_threshold_change_or_repair"

_LOCKED_TASK_DECISION_FIELDS = frozenset(
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

_TASK_ENDPOINT_EVALUATORS: Mapping[str, str] = {
    CELL_TASK: "cell_donor_balanced_macro_f1_v1",
    VARIANT_TASK: "variant_ld_block_fisher_z_spearman_gain_v1",
    RNA_ATAC_TASK: "rna_atac_two_way_deviance_reduction_v1",
    BULK_TASK: "bulk_paired_spearman_gain_v1",
    GRAPH_TASK: "graph_ld_block_auprc_gain_v1",
}
_PREDICTION_TABLE_FIELDS: Mapping[str, tuple[str, ...]] = {
    CELL_TASK: ("row_hash", "unit_hash", "predicted_class"),
    VARIANT_TASK: (
        "row_hash",
        "unit_hash",
        "block_hash",
        "stratum",
        "predicted",
    ),
    RNA_ATAC_TASK: (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "predicted",
    ),
    BULK_TASK: ("row_hash", "unit_hash", "predicted"),
    GRAPH_TASK: ("row_hash", "unit_hash", "block_hash", "predicted"),
}
_OUTCOME_TABLE_FIELDS: Mapping[str, tuple[str, ...]] = {
    CELL_TASK: ("row_hash", "unit_hash", "observed_class"),
    VARIANT_TASK: (
        "row_hash",
        "unit_hash",
        "block_hash",
        "stratum",
        "observed",
    ),
    RNA_ATAC_TASK: (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "observed",
    ),
    BULK_TASK: ("row_hash", "unit_hash", "observed"),
    GRAPH_TASK: ("row_hash", "unit_hash", "block_hash", "observed_binary"),
}
_JOINED_ENDPOINT_FIELDS: Mapping[str, tuple[str, ...]] = {
    CELL_TASK: (
        "row_hash",
        "unit_hash",
        "observed_class",
        "candidate_class",
        "baseline_class",
    ),
    VARIANT_TASK: (
        "row_hash",
        "unit_hash",
        "block_hash",
        "stratum",
        "observed",
        "candidate",
        "baseline",
    ),
    RNA_ATAC_TASK: (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "observed",
        "candidate",
        "baseline",
    ),
    BULK_TASK: ("row_hash", "unit_hash", "observed", "candidate", "baseline"),
    GRAPH_TASK: (
        "row_hash",
        "unit_hash",
        "block_hash",
        "observed_binary",
        "candidate",
        "baseline",
    ),
}
# Numeric prediction columns averaged across seeds to form the development
# ensemble.  The cell task is empty because its ensemble is built from class
# probabilities and then argmaxed, exactly as held-back inference does.
_ENSEMBLE_AVERAGED_FIELDS: Mapping[str, tuple[str, ...]] = {
    CELL_TASK: (),
    VARIANT_TASK: ("candidate", "baseline"),
    RNA_ATAC_TASK: ("candidate", "baseline"),
    BULK_TASK: ("candidate", "baseline"),
    GRAPH_TASK: ("candidate", "baseline"),
}
_SEALED_ENSEMBLE_FIELDS: Mapping[str, tuple[str, ...]] = {
    CELL_TASK: (
        "row_hash",
        "unit_hash",
        "observed_class",
        "candidate_class",
        "baseline_class",
        "candidate_probabilities_json",
        "baseline_probabilities_json",
    ),
    VARIANT_TASK: (
        "row_hash",
        "unit_hash",
        "block_hash",
        "stratum",
        "observed",
        "observed_binary",
        "summary_statistics_complete",
        "candidate",
        "baseline",
        "retrieval_scores_json",
    ),
}
_HASH_ROW_FIELDS = frozenset({"row_hash", "unit_hash", "donor_hash", "block_hash"})
_SELECTION_BOOTSTRAP_SEED = 20260821
_SELECTION_BOOTSTRAP_RESAMPLES = 1_000
_POWER_BOOTSTRAP_SEED = 20260821
_POWER_BOOTSTRAP_RESAMPLES = 10_000
_ROSTER_AUTHORITY_SCHEMA_VERSION = "masld-bench-evaluator-roster-v1"
_SEALED_METRIC_RESULT_CACHE: dict[str, Any] = {}


def _sealed_evaluator_source_bundle() -> dict[str, Any]:
    """Hash the local source closure that constructs held-back endpoint metrics."""

    package_root = Path(__file__).resolve().parent
    source_paths = (
        ("src/masld_bench/tournament.py", package_root / "tournament.py"),
        (
            "src/masld_bench/evaluators/sealed_metrics.py",
            package_root / "evaluators" / "sealed_metrics.py",
        ),
        (
            "src/masld_bench/evaluators/metrics.py",
            package_root / "evaluators" / "metrics.py",
        ),
        ("src/masld_bench/hashing.py", package_root / "hashing.py"),
    )
    records: list[dict[str, Any]] = []
    for relative_path, source_path in source_paths:
        if not source_path.is_file() or source_path.is_symlink():
            raise TournamentError(
                f"sealed evaluator source is missing or symlinked: {relative_path}"
            )
        records.append(
            {
                "path": relative_path,
                "sha256": sha256_file(source_path),
                "size_bytes": source_path.stat().st_size,
            }
        )
    return {
        "schema_version": SEALED_EVALUATOR_SOURCE_BUNDLE_SCHEMA_VERSION,
        "sources": records,
        "bundle_sha256": canonical_sha256(records),
    }


def _load_gate_policy() -> dict[str, Any]:
    try:
        payload = _PROMOTION_CONFIG.read_bytes()
    except OSError as exc:
        raise RuntimeError(f"cannot read promotion-gate registry: {exc}") from exc
    observed_sha256 = sha256(payload).hexdigest()
    if observed_sha256 != _PROMOTION_CONFIG_SHA256:
        raise RuntimeError(
            "promotion-gate registry hash mismatch: "
            f"expected {_PROMOTION_CONFIG_SHA256}, observed {observed_sha256}"
        )
    try:
        raw = tomllib.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise RuntimeError(f"invalid promotion-gate registry: {exc}") from exc
    if raw.get("schema_version") != "masld-bench-promotion-gates-v1":
        raise RuntimeError("unexpected promotion-gate schema_version")
    tasks = raw.get("tasks")
    if not isinstance(tasks, dict) or set(tasks) != EXPECTED_TASK_IDS:
        observed = sorted(tasks) if isinstance(tasks, dict) else []
        raise RuntimeError(
            "promotion-gate task IDs diverge from the executable policy: "
            f"expected={sorted(EXPECTED_TASK_IDS)!r}, observed={observed!r}"
        )
    if raw.get("confirmatory_method") != "Holm":
        raise RuntimeError("confirmatory method must remain Holm")
    if float(raw.get("confirmatory_fwer", -1.0)) != 0.05:
        raise RuntimeError("confirmatory FWER must remain 0.05")
    for task_id, policy in tasks.items():
        if not isinstance(policy.get("promotion_gate_id"), str):
            raise RuntimeError(f"{task_id} lacks a promotion_gate_id")
        has_power_metric = isinstance(policy.get("power_metric"), str)
        has_power_method = isinstance(policy.get("power_method_id"), str)
        if has_power_metric != has_power_method:
            raise RuntimeError(
                f"{task_id} must bind power_metric and power_method_id together"
            )
    return raw


PROMOTION_GATE_POLICY = _load_gate_policy()
CHAMPION_GATE_THRESHOLDS: dict[str, Mapping[str, Any]] = dict(
    PROMOTION_GATE_POLICY["tasks"]
)


class TournamentError(ValueError):
    """Base exception for an invalid or disallowed tournament operation."""


class SelectionLockError(TournamentError):
    """Raised when a strict aggregate SelectionLock cannot be verified."""


class TerminalFailureError(TournamentError):
    """Raised when an opened bundle or terminally failed model is reused."""


def _safe_resolve(path: str | Path, label: str) -> Path:
    try:
        return reject_symlink_components(Path(path), label=label).resolve()
    except (ArtifactError, OSError) as exc:
        raise TournamentError(f"invalid {label}: {exc}") from exc


@dataclass(frozen=True, slots=True)
class Objective:
    name: str
    maximize: bool = True
    tolerance: float = 0.0


@dataclass(frozen=True, slots=True)
class GateCheck:
    name: str
    passed: bool
    observed: Any
    threshold: Any
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class DevelopmentGateDecision:
    decision_sha256: str
    task_id: str
    selected_run_ids: tuple[str, ...]
    selection_lock_id: str
    passed: bool
    created_before_outcome_unblind: bool
    source_family_id: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class ComplementarityDecision:
    decision_id: str
    triggered: bool
    evidence_receipt_sha256: str
    complementarity_evidence_sha256: str
    pairwise_correlations_sha256: str
    shortlist_ids: tuple[str, ...]
    passing_tracks: tuple[str, ...]
    source_families: tuple[str, ...]
    selected_candidate_ids: tuple[str, ...]
    complementary_pairs: tuple[tuple[str, str], ...]
    created_before_outcome_unblind: bool
    reasons: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision_id": self.decision_id,
            "triggered": self.triggered,
            "evidence_receipt_sha256": self.evidence_receipt_sha256,
            "complementarity_evidence_sha256": self.complementarity_evidence_sha256,
            "pairwise_correlations_sha256": self.pairwise_correlations_sha256,
            "shortlist_ids": list(self.shortlist_ids),
            "passing_tracks": list(self.passing_tracks),
            "source_families": list(self.source_families),
            "selected_candidate_ids": list(self.selected_candidate_ids),
            "complementary_pairs": [list(pair) for pair in self.complementary_pairs],
            "created_before_outcome_unblind": self.created_before_outcome_unblind,
            "reasons": list(self.reasons),
        }


def _read(record: Any, name: str, default: Any = _NO_DEFAULT) -> Any:
    if isinstance(record, Mapping):
        if name in record:
            return record[name]
    elif hasattr(record, name):
        return getattr(record, name)
    if default is _NO_DEFAULT:
        raise KeyError(name)
    return default


def _read_first(record: Any, names: Sequence[str], default: Any = _MISSING) -> Any:
    for name in names:
        value = _read(record, name, _MISSING)
        if value is not _MISSING:
            return value
    if default is _MISSING:
        raise KeyError("|".join(names))
    return default


def _nonempty_identifier(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise TournamentError(f"{name} must be a non-empty string")
    return value.strip()


def _sha256_identifier(value: Any, name: str) -> str:
    identifier = _nonempty_identifier(value, name)
    if not _SHA256.fullmatch(identifier):
        raise TournamentError(
            f"{name} must contain 64 lowercase hexadecimal characters"
        )
    return identifier


def _finite_number(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise TournamentError(f"{name} must be numeric, not boolean")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise TournamentError(f"{name} must be numeric") from exc
    if not isfinite(number):
        raise TournamentError(f"{name} must be finite")
    return number


def _unique_strings(
    value: Any, name: str, *, allow_empty: bool = False, sorted_required: bool = False
) -> tuple[str, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TournamentError(f"{name} must be an array of strings")
    result = tuple(_nonempty_identifier(item, f"{name}[]") for item in value)
    if not allow_empty and not result:
        raise TournamentError(f"{name} cannot be empty")
    if len(set(result)) != len(result):
        raise TournamentError(f"{name} cannot contain duplicates")
    if sorted_required and result != tuple(sorted(result)):
        raise TournamentError(f"{name} must be canonically sorted")
    return result


def _sha256_sequence(value: Any, name: str) -> tuple[str, ...]:
    values = _unique_strings(value, name)
    return tuple(
        _sha256_identifier(item, f"{name}[{index}]")
        for index, item in enumerate(values)
    )


def _candidate_id(candidate: Any) -> str:
    return _nonempty_identifier(
        _read_first(candidate, ("candidate_id", "run_id", "id"), None),
        "candidate_id",
    )


def _candidate_family(candidate: Any) -> str:
    return _nonempty_identifier(
        _read_first(candidate, ("family", "model_family"), None), "family"
    )


def _candidate_model_id(candidate: Any) -> str:
    """Return a model identity, falling back only for legacy generic records."""

    return _nonempty_identifier(
        _read(candidate, "model_id", _candidate_id(candidate)), "model_id"
    )


def _candidate_eligible(candidate: Any) -> bool:
    value = _read(candidate, "eligible", _MISSING)
    if not isinstance(value, bool):
        raise TournamentError(
            f"eligible must be an explicit boolean for {_candidate_id(candidate)}"
        )
    return value


def _candidate_metric(candidate: Any, name: str) -> float:
    metrics = _read(candidate, "metrics", {})
    if not isinstance(metrics, Mapping):
        raise TournamentError(f"metrics must be a mapping for {_candidate_id(candidate)}")
    value = metrics.get(name, _MISSING)
    if value is _MISSING:
        value = _read(candidate, name, _MISSING)
    if value is _MISSING:
        raise TournamentError(
            f"candidate {_candidate_id(candidate)!r} is missing metric {name!r}"
        )
    return _finite_number(value, f"{_candidate_id(candidate)}.{name}")


def _candidate_standard_error(candidate: Any, metric: str) -> float:
    errors = _read(candidate, "standard_errors", {})
    if not isinstance(errors, Mapping):
        raise TournamentError(
            f"standard_errors must be a mapping for {_candidate_id(candidate)}"
        )
    value = errors.get(metric, _MISSING)
    if value is _MISSING:
        value = _read(candidate, f"{metric}_se", _MISSING)
    if value is _MISSING:
        raise TournamentError(
            f"candidate {_candidate_id(candidate)!r} lacks a standard error for {metric!r}"
        )
    result = _finite_number(value, f"{_candidate_id(candidate)}.{metric}.standard_error")
    if result < 0.0:
        raise TournamentError("standard errors cannot be negative")
    return result


def _candidate_complexity(candidate: Any) -> float:
    value = _read(candidate, "complexity", _MISSING)
    if value is _MISSING:
        metrics = _read(candidate, "metrics", {})
        value = metrics.get("complexity", metrics.get("compute_cost", inf))
    if value == inf:
        return inf
    result = _finite_number(value, f"{_candidate_id(candidate)}.complexity")
    if result < 0.0:
        raise TournamentError("complexity cannot be negative")
    return result


def _validated_candidates(candidates: Iterable[Any]) -> list[Any]:
    records = list(candidates)
    seen: set[str] = set()
    for candidate in records:
        identifier = _candidate_id(candidate)
        if identifier in seen:
            raise TournamentError(f"duplicate candidate_id: {identifier!r}")
        seen.add(identifier)
    return records


def _normalized_finite_mapping(value: Any, label: str) -> dict[str, float]:
    if not isinstance(value, Mapping) or not value:
        raise TournamentError(f"{label} must be a non-empty mapping")
    result = {
        _nonempty_identifier(str(key), f"{label} key"): _finite_number(
            item, f"{label}[{key}]"
        )
        for key, item in value.items()
    }
    return dict(sorted(result.items()))


def _candidate_plan_source(path: str | Path) -> tuple[Path, dict[str, Any], dict[str, str]]:
    root = reject_symlink_components(Path(path), label="candidate campaign")
    if not root.is_absolute() or root.is_symlink() or not root.is_dir():
        raise TournamentError("candidate campaign must be an absolute non-symlink directory")
    root = root.resolve()
    try:
        manifest = verify_frozen_tree(root)
        plan = load_frozen_plan(root)
    except (ArtifactError, PlanningError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid frozen candidate campaign: {exc}") from exc
    metadata = manifest.get("metadata")
    if not isinstance(metadata, Mapping) or metadata.get("artifact_class") != "candidate_plan":
        raise TournamentError("scientific receipt source is not a candidate_plan")
    plan_sha256 = _sha256_identifier(plan.get("plan_sha256"), "plan.plan_sha256")
    binding = {
        "path": root.as_posix(),
        "manifest_sha256": sha256_file(root / "ARTIFACTS.json"),
        "document_sha256": sha256_file(root / "plan.json"),
        "plan_sha256": plan_sha256,
    }
    return root, plan, binding


def _run_from_plan(plan: Mapping[str, Any], run_id: str) -> Mapping[str, Any]:
    runs = plan.get("runs")
    if not isinstance(runs, list):
        raise TournamentError("candidate plan runs must be an array")
    matches = [run for run in runs if isinstance(run, Mapping) and run.get("run_id") == run_id]
    if len(matches) != 1:
        raise TournamentError(f"run_id must identify one candidate-plan run: {run_id}")
    raw = dict(matches[0])
    claimed = _sha256_identifier(raw.pop("run_id", None), "run.run_id")
    if canonical_sha256(raw) != claimed:
        raise TournamentError("candidate-plan run identity hash mismatch")
    return matches[0]


def _model_disposition_for_plan(
    plan: Mapping[str, Any], model_id: str
) -> Mapping[str, Any]:
    records = plan.get("model_dispositions")
    if not isinstance(records, list):
        raise TournamentError("candidate plan model_dispositions must be an array")
    matches = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("model_id") == model_id
    ]
    if len(matches) != 1:
        raise TournamentError(f"candidate plan lacks one disposition for {model_id}")
    _nonempty_identifier(matches[0].get("family_id"), f"{model_id}.family_id")
    return matches[0]


def _run_dataset_registry_bindings(
    plan: Mapping[str, Any], run: Mapping[str, Any]
) -> dict[str, str]:
    dataset_ids = _unique_strings(run.get("dataset_ids"), "run.dataset_ids")
    locks = plan.get("dataset_locks")
    if not isinstance(locks, Mapping) or not set(dataset_ids).issubset(locks):
        raise TournamentError(
            "candidate plan dataset_locks do not cover the scientific run"
        )
    immutable_inputs = run.get("immutable_inputs")
    if not isinstance(immutable_inputs, Mapping) or set(immutable_inputs) != set(
        dataset_ids
    ):
        raise TournamentError(
            "scientific run immutable_inputs must exactly cover its datasets"
        )
    registry: dict[str, str] = {}
    for dataset_id in dataset_ids:
        lock = locks[dataset_id]
        if not isinstance(lock, Mapping):
            raise TournamentError(f"dataset lock {dataset_id!r} must be an object")
        registry[dataset_id] = _sha256_identifier(
            lock.get("registry_sha256"),
            f"dataset_locks[{dataset_id}].registry_sha256",
        )
        activation_sha256 = _sha256_identifier(
            lock.get("activation_sha256"),
            f"dataset_locks[{dataset_id}].activation_sha256",
        )
        if immutable_inputs[dataset_id] != activation_sha256:
            raise TournamentError(
                f"scientific run changed activation binding for {dataset_id}"
            )
    return dict(sorted(registry.items()))


def _verify_execution_attempt(path: str | Path) -> tuple[Path, dict[str, Any], dict[str, str]]:
    root = _safe_resolve(path, "run execution attempt")
    try:
        from .campaign import verify_run_execution_attempt

        verify_run_execution_attempt(root, require_succeeded=True)
        binding, _ = _frozen_document_binding(
            root,
            filename="run_execution_receipt.json",
            artifact_class="run_execution_attempt",
        )
        payload = json.loads(
            (root / "run_execution_receipt.json").read_text(encoding="utf-8")
        )
    except (ImportError, AttributeError) as exc:
        raise TournamentError(
            "campaign verifier does not expose verify_run_execution_attempt"
        ) from exc
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TournamentError(f"invalid successful run execution attempt: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise TournamentError("run execution receipt must be a JSON object")
    if (
        payload.get("schema_version") != "masld-bench-run-execution-receipt-v1"
        or payload.get("status") != "succeeded"
    ):
        raise TournamentError("scientific receipt requires a succeeded execution receipt")
    return root, dict(payload), binding


def _terminal_execution_attempt_binding(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    root = _safe_resolve(path, "terminal run execution attempt")
    try:
        from .campaign import verify_run_execution_attempt

        payload = verify_run_execution_attempt(root, require_succeeded=False)
        binding, _ = _frozen_document_binding(
            root,
            filename="run_execution_receipt.json",
            artifact_class="run_execution_attempt",
        )
    except (ImportError, AttributeError) as exc:
        raise TournamentError(
            "campaign verifier does not expose verify_run_execution_attempt"
        ) from exc
    except (ArtifactError, OSError, ValueError, RuntimeError) as exc:
        raise TournamentError(f"invalid terminal execution attempt: {exc}") from exc
    status = str(payload.get("status", ""))
    if status not in {
        "invalid_output",
        "failed_software",
        "failed_resource",
        "blocked_upstream",
        "blocked_license",
        "blocked_data",
        "skipped_gate",
    }:
        raise TournamentError(
            "terminal selection disposition must be a frozen non-success state"
        )
    return (
        {
            **binding,
            "run_id": _sha256_identifier(payload.get("run_id"), "terminal run_id"),
            "plan_sha256": _sha256_identifier(
                payload.get("plan_sha256"), "terminal plan_sha256"
            ),
            "status": status,
            "attempt": int(payload.get("attempt", 0)),
        },
        dict(payload),
    )


def _table_schema_sha256(fields: Sequence[str]) -> str:
    return canonical_sha256({"format": "tsv", "fields": list(fields)})


def _read_exact_tsv(path: Path, fields: Sequence[str], label: str) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != tuple(fields):
                raise TournamentError(
                    f"{label} header must be exactly {' '.join(fields)}"
                )
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise TournamentError(f"cannot read {label}: {exc}") from exc
    if not rows or any(
        None in row
        or any(value is None or any(char in value for char in "\t\r\n") for value in row.values())
        for row in rows
    ):
        raise TournamentError(f"{label} must contain non-empty exact TSV rows")
    row_ids = [row[fields[0]] for row in rows]
    if any(_SHA256.fullmatch(value) is None for value in row_ids):
        raise TournamentError(f"{label} row IDs must be lowercase SHA-256 values")
    if len(set(row_ids)) != len(row_ids) or row_ids != sorted(row_ids):
        raise TournamentError(f"{label} row IDs must be unique and canonically sorted")
    for field in set(fields) & _HASH_ROW_FIELDS:
        if any(_SHA256.fullmatch(row[field]) is None for row in rows):
            raise TournamentError(f"{label}.{field} must contain lowercase SHA-256 values")
    return rows


def _tsv_text(fields: Sequence[str], rows: Sequence[Mapping[str, str]]) -> str:
    lines = ["\t".join(fields)]
    lines.extend("\t".join(str(row[field]) for field in fields) for row in rows)
    return "\n".join(lines) + "\n"


def freeze_development_outcome_bundle(
    *,
    task_id: str,
    dataset_ids: Sequence[str],
    dataset_registry_sha256s: Mapping[str, str],
    split_id: str,
    row_id_field: str,
    unit_id_field: str,
    unit_id_namespace: str,
    biological_unit: str,
    outcome_table_path: str | Path,
    output_root: str | Path,
) -> Path:
    """Freeze known development outcomes under the same join identity as predictions."""

    if task_id not in _TASK_ENDPOINT_EVALUATORS:
        raise TournamentError(f"task {task_id!r} has no registered scientific evaluator")
    fields = _OUTCOME_TABLE_FIELDS[task_id]
    if (row_id_field, unit_id_field) != fields[:2]:
        raise TournamentError("development outcome row/unit fields differ from task schema")
    datasets = _unique_strings(dataset_ids, "dataset_ids")
    if not isinstance(dataset_registry_sha256s, Mapping) or set(dataset_registry_sha256s) != set(datasets):
        raise TournamentError("dataset_registry_sha256s must cover outcome datasets exactly")
    dataset_bindings = {
        dataset_id: _sha256_identifier(
            dataset_registry_sha256s[dataset_id], f"dataset {dataset_id}"
        )
        for dataset_id in datasets
    }
    split_id = _nonempty_identifier(split_id, "split_id")
    unit_id_namespace = _nonempty_identifier(unit_id_namespace, "unit_id_namespace")
    biological_unit = _nonempty_identifier(biological_unit, "biological_unit")
    rows = _read_exact_tsv(
        _safe_resolve(outcome_table_path, "development outcome table"),
        fields,
        "development outcome table",
    )
    join_key = canonical_sha256(
        {
            "task_id": task_id,
            "dataset_ids": list(datasets),
            "split_id": split_id,
            "row_id_field": row_id_field,
            "unit_id_field": unit_id_field,
            "unit_id_namespace": unit_id_namespace,
            "biological_unit": biological_unit,
        }
    )
    output = _safe_resolve(output_root, "development outcome output root")
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".development-outcomes.", dir=output))
    try:
        table_path = staging / "outcomes.tsv"
        row_path = staging / "row_ids.tsv"
        write_text_exclusive(table_path, _tsv_text(fields, rows), mode=0o440)
        row_fields = (row_id_field, unit_id_field)
        write_text_exclusive(
            row_path,
            _tsv_text(row_fields, rows),
            mode=0o440,
        )
        table_ref = ArtifactRef.from_path(
            table_path,
            relative_to=staging,
            media_type="text/tab-separated-values",
            role=f"standardized_development_outcomes:{task_id}",
        )
        row_ref = ArtifactRef.from_path(
            row_path,
            relative_to=staging,
            media_type="text/tab-separated-values",
            role=f"development_outcome_row_ids:{task_id}",
        )
        identity = {
            "schema_version": DEVELOPMENT_OUTCOME_BUNDLE_SCHEMA_VERSION,
            "task_id": task_id,
            "dataset_ids": list(datasets),
            "dataset_registry_sha256s": dataset_bindings,
            "split_id": split_id,
            "row_id_field": row_id_field,
            "unit_id_field": unit_id_field,
            "unit_id_namespace": unit_id_namespace,
            "biological_unit": biological_unit,
            "table_schema_sha256": _table_schema_sha256(fields),
            "source_join_key_sha256": join_key,
            "n_outcomes": len(rows),
            "standardized_table": table_ref.to_dict(),
            "row_ids": row_ref.to_dict(),
            "created_before_outcome_unblind": True,
        }
        payload = {"outcome_bundle_id": canonical_sha256(identity), **identity}
        write_json_exclusive(
            staging / "development_outcome_bundle.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "development_outcome_bundle",
                "outcome_bundle_id": payload["outcome_bundle_id"],
            },
        )
        target = output / f"development-outcomes--{payload['outcome_bundle_id']}"
        if target.exists():
            raise TournamentError(f"development outcome bundle already exists: {target}")
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_development_outcome_bundle(target)
    return target


def verify_development_outcome_bundle(path: str | Path) -> dict[str, Any]:
    root = _safe_resolve(path, "development outcome bundle")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "development_outcome_bundle.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid frozen development outcome bundle: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != DEVELOPMENT_OUTCOME_BUNDLE_SCHEMA_VERSION:
        raise TournamentError("unsupported development outcome bundle schema")
    expected_fields = {
        "outcome_bundle_id",
        "schema_version",
        "task_id",
        "dataset_ids",
        "dataset_registry_sha256s",
        "split_id",
        "row_id_field",
        "unit_id_field",
        "unit_id_namespace",
        "biological_unit",
        "table_schema_sha256",
        "source_join_key_sha256",
        "n_outcomes",
        "standardized_table",
        "row_ids",
        "created_before_outcome_unblind",
    }
    if set(payload) != expected_fields:
        raise TournamentError("development outcome bundle has the wrong exact schema")
    task_id = str(payload.get("task_id"))
    if task_id not in _TASK_ENDPOINT_EVALUATORS:
        raise TournamentError("development outcome task has no registered evaluator")
    fields = _OUTCOME_TABLE_FIELDS[task_id]
    try:
        table_ref = ArtifactRef.from_dict(payload.get("standardized_table"))
        row_ref = ArtifactRef.from_dict(payload.get("row_ids"))
    except (ContractError, TypeError) as exc:
        raise TournamentError(f"invalid development outcome ArtifactRef: {exc}") from exc
    if (
        table_ref.role != f"standardized_development_outcomes:{task_id}"
        or row_ref.role != f"development_outcome_row_ids:{task_id}"
    ):
        raise TournamentError("development outcome artifact roles are invalid")
    table_path = table_ref.validate(root, require_relative=True)
    row_path = row_ref.validate(root, require_relative=True)
    rows = _read_exact_tsv(table_path, fields, "development outcome table")
    row_fields = (str(payload.get("row_id_field")), str(payload.get("unit_id_field")))
    inventory = _read_exact_tsv(row_path, row_fields, "development outcome row inventory")
    if [tuple(row[field] for field in row_fields) for row in rows] != [
        tuple(row[field] for field in row_fields) for row in inventory
    ]:
        raise TournamentError("development outcome row inventory differs from its table")
    datasets = _unique_strings(payload.get("dataset_ids"), "dataset_ids")
    bindings = payload.get("dataset_registry_sha256s")
    if not isinstance(bindings, Mapping) or set(bindings) != set(datasets):
        raise TournamentError("development outcome dataset bindings are incomplete")
    for dataset_id in datasets:
        _sha256_identifier(
            bindings[dataset_id],
            f"dataset_registry_sha256s[{dataset_id}]",
        )
    if (
        (payload.get("row_id_field"), payload.get("unit_id_field")) != fields[:2]
        or not isinstance(payload.get("n_outcomes"), int)
        or isinstance(payload.get("n_outcomes"), bool)
    ):
        raise TournamentError("development outcome row/unit/count contract changed")
    for field in ("split_id", "unit_id_namespace", "biological_unit"):
        _nonempty_identifier(payload.get(field), field)
    identity = dict(payload)
    claimed = _sha256_identifier(identity.pop("outcome_bundle_id", None), "outcome_bundle_id")
    expected_join = canonical_sha256(
        {
            "task_id": task_id,
            "dataset_ids": list(datasets),
            "split_id": payload.get("split_id"),
            "row_id_field": payload.get("row_id_field"),
            "unit_id_field": payload.get("unit_id_field"),
            "unit_id_namespace": payload.get("unit_id_namespace"),
            "biological_unit": payload.get("biological_unit"),
        }
    )
    if (
        payload.get("table_schema_sha256") != _table_schema_sha256(fields)
        or payload.get("source_join_key_sha256") != expected_join
        or payload.get("n_outcomes") != len(rows)
        or payload.get("created_before_outcome_unblind") is not True
        or canonical_sha256(identity) != claimed
    ):
        raise TournamentError("development outcome bundle identity does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": "development_outcome_bundle",
        "outcome_bundle_id": claimed,
    }:
        raise TournamentError("development outcome manifest metadata mismatch")
    return dict(payload)


def _prediction_source(
    bundle_path: str | Path,
    *,
    run: Mapping[str, Any],
    class_roster: Sequence[str] | None = None,
) -> tuple[PredictionBundle, list[dict[str, str]], dict[str, Any]]:
    path = _safe_resolve(bundle_path, "development PredictionBundle")
    try:
        bundle = PredictionBundle.load_json(path)
        bundle.validate_artifacts(path.parent)
    except (OSError, ValueError, ContractError) as exc:
        raise TournamentError(f"invalid frozen development PredictionBundle: {exc}") from exc
    if bundle.run_id != run["run_id"] or bundle.task_id != run["task_id"] or bundle.model_id != run["model_id"]:
        raise TournamentError("PredictionBundle does not bind its successful execution run")
    if (
        tuple(bundle.dataset_ids) != tuple(run["dataset_ids"])
        or bundle.split_id != run["split_id"]
    ):
        raise TournamentError("PredictionBundle changed the frozen run dataset/split")
    fields = _PREDICTION_TABLE_FIELDS.get(bundle.task_id)
    if fields is None:
        raise TournamentError(f"task {bundle.task_id!r} has no registered scientific evaluator")
    if bundle.table_schema_sha256 != _table_schema_sha256(fields):
        raise TournamentError("PredictionBundle table schema is not task-native")
    table_path = bundle.standardized_table.validate(path.parent, require_relative=True)
    row_path = bundle.row_ids.validate(path.parent, require_relative=True)
    rows = _read_exact_tsv(table_path, fields, f"PredictionBundle {bundle.bundle_id}")
    row_fields = (bundle.row_id_field, bundle.unit_id_field)
    inventory = _read_exact_tsv(
        row_path, row_fields, f"PredictionBundle {bundle.bundle_id} row inventory"
    )
    if [tuple(row[field] for field in row_fields) for row in rows] != [
        tuple(row[field] for field in row_fields) for row in inventory
    ]:
        raise TournamentError("PredictionBundle row inventory differs from its table")
    if bundle.n_predictions != len(rows) or bundle.missing_state is not MissingState.OBSERVED:
        raise TournamentError("PredictionBundle is incomplete for scientific scoring")
    binding = {
        "path": path.as_posix(),
        "document_sha256": sha256_file(path),
        "bundle_id": bundle.bundle_id,
        "run_id": bundle.run_id,
        "standardized_table_sha256": bundle.standardized_table.sha256,
        "row_ids_sha256": bundle.row_ids.sha256,
    }
    if bundle.task_id == CELL_TASK:
        # Cell-state ensembling averages class probabilities and then argmaxes,
        # exactly as held-back inference does.  A development bundle that ships
        # only hard labels cannot be ensemble-ranked, so require the output file
        # here rather than discovering the gap at selection time.
        if class_roster is None:
            raise TournamentError(
                "cell PredictionBundles require the frozen class roster"
            )
        probability_rows = _cell_probability_rows(
            bundle=bundle,
            bundle_path=path,
            prediction_rows=rows,
            class_roster=class_roster,
        )
        hard_by_row = {
            str(row["row_hash"]): str(row["predicted_class"]) for row in rows
        }
        for row in probability_rows:
            values = {
                class_id: _finite_number(
                    row["probabilities"][class_id],
                    f"{bundle.bundle_id}.probability::{class_id}",
                )
                for class_id in class_roster
            }
            if any(not 0.0 <= value <= 1.0 for value in values.values()):
                raise TournamentError(
                    "development class probabilities must lie in [0, 1]"
                )
            if abs(fsum(values.values()) - 1.0) > _ENSEMBLE_PROBABILITY_TOLERANCE:
                raise TournamentError(
                    "development class probabilities must sum to one"
                )
            winner = min(class_roster, key=lambda item: (-values[item], item))
            if hard_by_row.get(str(row["row_hash"])) != winner:
                raise TournamentError(
                    "development cell PredictionBundle hard class differs from "
                    "its probability argmax"
                )
        probability_path = _one_bundle_artifact(
            bundle, root=path.parent, role=f"class_probabilities:{CELL_TASK}"
        )
        binding["class_probabilities_sha256"] = sha256_file(probability_path)
    return bundle, rows, binding


def _variant_secondary_run_capability(
    run: Mapping[str, Any],
) -> tuple[dict[str, Any], str, str]:
    if run.get("task_id") != VARIANT_TASK:
        raise TournamentError("secondary run receipt is restricted to variant_to_regulation")
    model_id = _nonempty_identifier(run.get("model_id"), "secondary run model_id")
    if model_id not in _VARIANT_SECONDARY_MODELS:
        raise TournamentError(
            f"secondary run receipt does not admit model {model_id!r}"
        )
    metadata = run.get("metadata")
    if not isinstance(metadata, Mapping):
        raise TournamentError("secondary variant run lacks immutable capability metadata")
    if metadata.get("primary_endpoint_scoring_allowed") is not False:
        raise TournamentError(
            "secondary variant run must forbid primary endpoint scoring"
        )
    if metadata.get("primary_evaluator_id") != _VARIANT_PRIMARY_EVALUATOR_ID:
        raise TournamentError("secondary variant run names the wrong task evaluator")
    capability = metadata.get("variant_capability")
    capability_fields = {
        "model_id",
        "role",
        "native_outputs",
        "allowed_endpoints",
        "primary_eligible",
        "requires_fitted_head",
        "requires_observed_target_context",
        "is_mandatory_baseline",
    }
    if not isinstance(capability, Mapping) or set(capability) != capability_fields:
        raise TournamentError("secondary variant capability has an invalid exact schema")
    normalized = dict(capability)
    expected_values = {
        "model_id": model_id,
        "role": _VARIANT_SECONDARY_ROLES[model_id],
        "native_outputs": ["enhancer_gene_link_score"],
        "allowed_endpoints": ["enhancer_gene_link", _VARIANT_SECONDARY_ENDPOINT_ID],
        "primary_eligible": False,
        "requires_fitted_head": False,
        "requires_observed_target_context": False,
        "is_mandatory_baseline": True,
    }
    if normalized != expected_values:
        raise TournamentError(
            f"secondary variant capability for {model_id} is not canonical"
        )
    capability_sha256 = _sha256_identifier(
        metadata.get("variant_capability_sha256"),
        f"secondary variant capability SHA-256 for {model_id}",
    )
    if capability_sha256 != canonical_sha256(normalized):
        raise TournamentError(
            f"secondary variant capability SHA-256 changed for {model_id}"
        )
    registry_sha256 = _sha256_identifier(
        metadata.get("variant_capability_registry_sha256"),
        f"secondary variant capability registry SHA-256 for {model_id}",
    )
    return normalized, capability_sha256, registry_sha256


def _variant_secondary_run_sources(
    *,
    candidate_campaign_dir: str | Path,
    execution_attempt_dir: str | Path,
    prediction_bundle_path: str | Path,
) -> dict[str, Any]:
    _, plan, campaign_binding = _candidate_plan_source(candidate_campaign_dir)
    _, execution, execution_binding = _verify_execution_attempt(execution_attempt_dir)
    run = _run_from_plan(
        plan,
        _sha256_identifier(execution.get("run_id"), "secondary execution run_id"),
    )
    for field in (
        "plan_sha256",
        "candidate_manifest_sha256",
        "task_id",
        "model_id",
        "seed",
        "fold",
    ):
        if field == "plan_sha256":
            expected = campaign_binding["plan_sha256"]
        elif field == "candidate_manifest_sha256":
            expected = campaign_binding["manifest_sha256"]
        else:
            expected = run[field]
        if execution.get(field) != expected:
            raise TournamentError(
                f"secondary execution receipt changed frozen {field}"
            )
    capability, capability_sha256, capability_registry_sha256 = (
        _variant_secondary_run_capability(run)
    )
    model_id = str(run["model_id"])
    disposition = _model_disposition_for_plan(plan, model_id)
    if disposition.get("disposition") == "blocked":
        raise TournamentError(f"secondary variant model is blocked: {model_id}")
    dataset_registry_sha256s = _run_dataset_registry_bindings(plan, run)
    immutable_inputs = run.get("immutable_inputs")
    if not isinstance(immutable_inputs, Mapping):
        raise TournamentError("secondary variant run immutable_inputs must be an object")
    bundle, rows, prediction_binding = _prediction_source(
        prediction_bundle_path, run=run
    )
    expected_prediction_metadata = {
        "variant_endpoint_id": _VARIANT_SECONDARY_ENDPOINT_ID,
        "variant_evaluator_id": _VARIANT_SECONDARY_EVALUATOR_ID,
        "variant_score_transform_id": _VARIANT_SECONDARY_SCORE_TRANSFORM_ID,
        "primary_endpoint_scoring_allowed": False,
        "variant_capability_sha256": capability_sha256,
    }
    if dict(bundle.metadata) != expected_prediction_metadata:
        raise TournamentError(
            "secondary variant PredictionBundle endpoint metadata differs from its run"
        )
    for index, row in enumerate(rows):
        _finite_number(row.get("predicted"), f"secondary prediction row {index}")
    return {
        "campaign_binding": campaign_binding,
        "execution_binding": execution_binding,
        "prediction_binding": prediction_binding,
        "run": run,
        "family_id": _nonempty_identifier(disposition.get("family_id"), "family_id"),
        "dataset_registry_sha256s": dataset_registry_sha256s,
        "immutable_inputs": dict(sorted(immutable_inputs.items())),
        "capability": capability,
        "capability_sha256": capability_sha256,
        "capability_registry_sha256": capability_registry_sha256,
        "bundle": bundle,
    }


def _ensemble_alignment_sha256(
    task_id: str, rows: Sequence[Mapping[str, str]]
) -> str:
    """Hash the order-sensitive, prediction-free identity of a joined table.

    ``unit_set_sha256`` and ``row_set_sha256`` are SET hashes, so they prove the
    five seeds cover the same donors and rows but say nothing about per-row
    correspondence.  Averaging predictions requires per-row correspondence, so
    this hash pins the row order together with every unit, block, stratum, and
    observed value.
    """

    fields = _JOINED_ENDPOINT_FIELDS.get(task_id)
    if fields is None:
        raise TournamentError(f"task {task_id!r} has no joined endpoint schema")
    prediction_fields = {
        "candidate",
        "baseline",
        "candidate_class",
        "baseline_class",
    }
    metadata_fields = tuple(
        field for field in fields if field not in prediction_fields
    )
    row_ids = [str(row["row_hash"]) for row in rows]
    if not row_ids or len(set(row_ids)) != len(row_ids) or row_ids != sorted(row_ids):
        raise TournamentError(
            "ensemble alignment requires unique, canonically sorted rows"
        )
    return canonical_sha256(
        {
            "schema_version": _ENSEMBLE_ALIGNMENT_SCHEMA_VERSION,
            "task_id": task_id,
            "fields": list(metadata_fields),
            "rows": [[str(row[field]) for field in metadata_fields] for row in rows],
        }
    )


def _cell_probability_index(
    *,
    receipt: Mapping[str, Any],
    binding_key: str,
    class_roster: Sequence[str],
) -> dict[str, dict[str, float]]:
    """Load one frozen development PredictionBundle's class probabilities."""

    binding = receipt.get(binding_key)
    if not isinstance(binding, Mapping) or not isinstance(binding.get("path"), str):
        raise TournamentError(
            f"scientific receipt lacks a {binding_key} for ensemble ranking"
        )
    bundle_path = _safe_resolve(binding["path"], "development PredictionBundle")
    try:
        bundle = PredictionBundle.load_json(bundle_path)
        bundle.validate_artifacts(bundle_path.parent)
    except (ContractError, OSError, ValueError) as exc:
        raise TournamentError(
            f"invalid development PredictionBundle for ensemble ranking: {exc}"
        ) from exc
    if str(bundle.task_id) != CELL_TASK:
        raise TournamentError("cell ensemble requires a cell PredictionBundle")
    prediction_rows = _read_exact_tsv(
        bundle.standardized_table.validate(
            bundle_path.parent, require_relative=True
        ),
        _PREDICTION_TABLE_FIELDS[CELL_TASK],
        f"development PredictionBundle {bundle.bundle_id}",
    )
    probability_rows = _cell_probability_rows(
        bundle=bundle,
        bundle_path=bundle_path,
        prediction_rows=prediction_rows,
        class_roster=class_roster,
    )
    hard_by_row = {
        str(row["row_hash"]): str(row["predicted_class"])
        for row in prediction_rows
    }
    index: dict[str, dict[str, float]] = {}
    for row in probability_rows:
        row_hash = str(row["row_hash"])
        values: dict[str, float] = {}
        for class_id in class_roster:
            value = _finite_number(
                row["probabilities"][class_id],
                f"{bundle.bundle_id}.probability::{class_id}",
            )
            if not 0.0 <= value <= 1.0:
                raise TournamentError(
                    "development class probabilities must lie in [0, 1]"
                )
            values[class_id] = value
        if abs(fsum(values.values()) - 1.0) > _ENSEMBLE_PROBABILITY_TOLERANCE:
            raise TournamentError(
                "development class probabilities must sum to one"
            )
        winner = min(class_roster, key=lambda item: (-values[item], item))
        if hard_by_row.get(row_hash) != winner:
            raise TournamentError(
                "development cell PredictionBundle hard class differs from its "
                "probability argmax"
            )
        index[row_hash] = values
    return index


def _sample_dispersion(values: Sequence[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = fsum(values) / len(values)
    return sqrt(
        fsum((value - mean) ** 2 for value in values) / (len(values) - 1)
    )


def _require_comparable_score_scale(
    *,
    task_id: str,
    averaged_fields: Sequence[str],
    seeds: Sequence[int],
    columns: Mapping[tuple[str, int], Sequence[float]],
) -> None:
    """Refuse to average raw scores the five seeds put on different scales."""

    for field in averaged_fields:
        dispersions = {
            seed: _sample_dispersion(columns[(field, seed)]) for seed in seeds
        }
        smallest = min(dispersions.values())
        largest = max(dispersions.values())
        if largest == 0.0:
            # Every seed predicts a constant.  Degenerate but self-consistent;
            # the endpoint evaluator reports it, this check does not.
            continue
        if smallest <= 0.0:
            raise TournamentError(
                f"{task_id}.{field} has a constant seed alongside varying "
                f"seeds; raw-score averaging is not comparable "
                f"({ENSEMBLE_SCORE_SCALE_POLICY})"
            )
        if largest / smallest > _ENSEMBLE_SCORE_SCALE_RATIO_MAX:
            raise TournamentError(
                f"{task_id}.{field} seed score dispersions differ by "
                f"{largest / smallest:.3g}x, above the prospectively frozen "
                f"{_ENSEMBLE_SCORE_SCALE_RATIO_MAX:g}x limit "
                f"({ENSEMBLE_SCORE_SCALE_POLICY})"
            )


def _development_ensemble_rows(
    *,
    task_id: str,
    per_seed: Sequence[tuple[int, Path, Mapping[str, Any]]],
    endpoint_parameters: Mapping[str, Any],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Build the exact five-seed development ensemble joined table.

    Returns the ensemble rows and an alignment/stability record.  The ensemble
    is the deployed output file, so this is the quantity development must rank.
    """

    fields = _JOINED_ENDPOINT_FIELDS.get(task_id)
    averaged_fields = _ENSEMBLE_AVERAGED_FIELDS.get(task_id)
    if fields is None or averaged_fields is None:
        raise TournamentError(
            f"task {task_id!r} has no development ensemble contract"
        )
    ordered = sorted(per_seed, key=lambda item: int(item[0]))
    seeds = tuple(int(item[0]) for item in ordered)
    if len(seeds) != _ENSEMBLE_SEED_COUNT or len(set(seeds)) != _ENSEMBLE_SEED_COUNT:
        raise TournamentError(
            "development ensemble requires exactly five distinct seeds"
        )
    tables: dict[int, list[dict[str, str]]] = {}
    alignments: set[str] = set()
    for seed, receipt_root, receipt in ordered:
        joined_ref = ArtifactRef.from_dict(receipt["joined_endpoint_table"])
        joined_path = joined_ref.validate(receipt_root, require_relative=True)
        rows = _read_exact_tsv(
            joined_path, fields, f"development joined endpoint seed {seed}"
        )
        alignments.add(_ensemble_alignment_sha256(task_id, rows))
        tables[seed] = rows
    if len(alignments) != 1:
        raise TournamentError(
            "development seed joined tables are not row-aligned"
        )
    ensemble_alignment_sha256 = next(iter(alignments))
    reference = tables[seeds[0]]
    ensemble_rows: list[dict[str, str]] = []

    if task_id == CELL_TASK:
        roster = tuple(
            _unique_strings(
                endpoint_parameters.get("class_roster"),
                "class_roster",
                sorted_required=True,
            )
        )
        candidate_probabilities: dict[int, dict[str, dict[str, float]]] = {}
        baseline_probabilities: dict[int, dict[str, dict[str, float]]] = {}
        for seed, _receipt_root, receipt in ordered:
            candidate_probabilities[seed] = _cell_probability_index(
                receipt=receipt,
                binding_key="candidate_prediction_binding",
                class_roster=roster,
            )
            baseline_probabilities[seed] = _cell_probability_index(
                receipt=receipt,
                binding_key="baseline_prediction_binding",
                class_roster=roster,
            )
        for row in reference:
            row_hash = str(row["row_hash"])
            record = {field: str(row[field]) for field in fields}
            for label, source in (
                ("candidate_class", candidate_probabilities),
                ("baseline_class", baseline_probabilities),
            ):
                mean_probabilities: dict[str, float] = {}
                for class_id in roster:
                    values = []
                    for seed in seeds:
                        seed_index = source[seed]
                        if row_hash not in seed_index:
                            raise TournamentError(
                                "development ensemble seed lacks a joined row"
                            )
                        values.append(seed_index[row_hash][class_id])
                    mean_probabilities[class_id] = fsum(values) / len(values)
                record[label] = min(
                    roster,
                    key=lambda item: (-mean_probabilities[item], item),
                )
            ensemble_rows.append(record)
    else:
        columns: dict[tuple[str, int], list[float]] = {}
        for field in averaged_fields:
            for seed in seeds:
                columns[(field, seed)] = [
                    _finite_number(row[field], f"{task_id}.{field}")
                    for row in tables[seed]
                ]
        _require_comparable_score_scale(
            task_id=task_id,
            averaged_fields=averaged_fields,
            seeds=seeds,
            columns=columns,
        )
        for index, row in enumerate(reference):
            record = {field: str(row[field]) for field in fields}
            for field in averaged_fields:
                values = []
                for seed in seeds:
                    seed_row = tables[seed][index]
                    if str(seed_row["row_hash"]) != str(row["row_hash"]):
                        raise TournamentError(
                            "development ensemble rows drifted across seeds"
                        )
                    values.append(columns[(field, seed)][index])
                record[field] = repr(fsum(values) / len(values))
            ensemble_rows.append(record)

    alignment = {
        "development_ensemble_policy_id": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "score_scale_policy": (
            ENSEMBLE_SCORE_SCALE_POLICY if averaged_fields else "not_applicable"
        ),
        "ensemble_alignment_sha256": ensemble_alignment_sha256,
        "ensemble_rows_sha256": canonical_sha256(
            {
                "task_id": task_id,
                "fields": list(fields),
                "rows": [
                    [row[field] for field in fields] for row in ensemble_rows
                ],
            }
        ),
        "seeds": list(seeds),
    }
    return ensemble_rows, alignment


def _join_development_rows(
    *,
    task_id: str,
    candidate_rows: Sequence[Mapping[str, str]],
    baseline_rows: Sequence[Mapping[str, str]],
    outcome_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    indices = [
        {row["row_hash"]: row for row in rows}
        for rows in (candidate_rows, baseline_rows, outcome_rows)
    ]
    if not (set(indices[0]) == set(indices[1]) == set(indices[2])):
        raise TournamentError("candidate, baseline, and outcome row-ID sets differ")
    joined: list[dict[str, str]] = []
    metadata_fields = [
        field
        for field in _OUTCOME_TABLE_FIELDS[task_id]
        if field not in {"observed", "observed_binary", "observed_class"}
    ]
    for row_hash in sorted(indices[0]):
        candidate, baseline, outcome = (index[row_hash] for index in indices)
        if any(
            candidate.get(field) != outcome.get(field)
            or baseline.get(field) != outcome.get(field)
            for field in metadata_fields
        ):
            raise TournamentError(f"development join metadata differs for row {row_hash}")
        if task_id == CELL_TASK:
            joined.append(
                {
                    "row_hash": row_hash,
                    "unit_hash": outcome["unit_hash"],
                    "observed_class": outcome["observed_class"],
                    "candidate_class": candidate["predicted_class"],
                    "baseline_class": baseline["predicted_class"],
                }
            )
        else:
            joined.append(
                {
                    **{field: outcome[field] for field in metadata_fields},
                    **{
                        field: outcome[field]
                        for field in ("observed", "observed_binary")
                        if field in outcome
                    },
                    "candidate": candidate["predicted"],
                    "baseline": baseline["predicted"],
                }
            )
    return joined


def _validated_roster_authority(
    *,
    plan: Mapping[str, Any],
    task_id: str,
    evaluator_id: str,
    roster_field: str,
    roster: Sequence[str],
    raw: Any,
) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TournamentError(f"{task_id}.roster_authority must be an ArtifactRef")
    try:
        authority = ArtifactRef.from_dict(raw)
    except ContractError as exc:
        raise TournamentError(f"invalid {task_id} roster authority: {exc}") from exc
    expected_role = f"task_evaluator_roster:{task_id}"
    if authority.role != expected_role:
        raise TournamentError(
            f"{task_id} roster authority role must be {expected_role!r}"
        )
    configured = Path(authority.path)
    config_root = None if configured.is_absolute() else plan.get("config_root")
    if config_root is not None and not isinstance(config_root, str):
        raise TournamentError("candidate plan config_root must be a string")
    try:
        source = authority.validate(config_root)
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (
        ContractError,
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
    ) as exc:
        raise TournamentError(f"invalid {task_id} roster authority: {exc}") from exc
    expected = {
        "schema_version": _ROSTER_AUTHORITY_SCHEMA_VERSION,
        "task_id": task_id,
        "primary_evaluator_id": evaluator_id,
        "roster_field": roster_field,
        "roster": list(roster),
    }
    if not isinstance(payload, Mapping) or dict(payload) != expected:
        raise TournamentError(
            f"{task_id} roster authority content differs from its frozen TaskSpec"
        )
    return authority.to_dict()


def _endpoint_contract_from_plan(
    plan: Mapping[str, Any], task_id: str
) -> tuple[str, dict[str, Any], str, dict[str, Any] | None]:
    records = [
        item
        for item in plan.get("task_dispositions", ())
        if isinstance(item, Mapping) and item.get("task_id") == task_id
    ]
    if len(records) != 1:
        raise TournamentError(f"candidate plan lacks one TaskSpec binding for {task_id}")
    record = records[0]
    evaluator_id = _nonempty_identifier(
        record.get("primary_evaluator_id"), f"{task_id}.primary_evaluator_id"
    )
    if evaluator_id != _TASK_ENDPOINT_EVALUATORS.get(task_id):
        raise TournamentError(f"{task_id} TaskSpec names the wrong primary evaluator")
    parameters = record.get("evaluator_parameters")
    if not isinstance(parameters, Mapping):
        raise TournamentError(f"{task_id} TaskSpec evaluator_parameters must be an object")
    expected_evaluator_contract_sha256 = canonical_sha256(
        {
            "primary_evaluator_id": evaluator_id,
            "evaluator_parameters": dict(parameters),
        }
    )
    if record.get("evaluator_contract_sha256") != expected_evaluator_contract_sha256:
        raise TournamentError(f"{task_id} TaskSpec evaluator contract hash changed")
    roster_field = (
        "class_roster"
        if task_id == CELL_TASK
        else "strata"
        if task_id in {VARIANT_TASK, RNA_ATAC_TASK}
        else None
    )
    expected_fields = (
        {roster_field, "roster_authority"}
        if roster_field is not None
        else set()
    )
    if set(parameters) != expected_fields:
        raise TournamentError(
            f"{task_id} TaskSpec evaluator parameters do not match the registered schema"
        )
    normalized: dict[str, Any] = {}
    roster_authority: dict[str, Any] | None = None
    if roster_field is not None:
        normalized[roster_field] = list(
            _unique_strings(
                parameters[roster_field],
                f"{task_id}.{roster_field}",
                sorted_required=True,
            )
        )
        roster_authority = _validated_roster_authority(
            plan=plan,
            task_id=task_id,
            evaluator_id=evaluator_id,
            roster_field=roster_field,
            roster=normalized[roster_field],
            raw=parameters["roster_authority"],
        )
    contract_sha256 = _sha256_identifier(
        record.get("contract_sha256"), f"{task_id}.contract_sha256"
    )
    return evaluator_id, normalized, contract_sha256, roster_authority


def _sample_standard_error(values: Sequence[float]) -> float:
    if len(values) < 2:
        raise TournamentError("standard-error estimation requires at least two resamples")
    mean = fsum(values) / len(values)
    return sqrt(fsum((value - mean) ** 2 for value in values) / (len(values) - 1))


def _endpoint_summary(endpoint: EndpointBootstrapResult) -> dict[str, Any]:
    return {
        "evaluator_id": endpoint.evaluator_id,
        "observed_effect": endpoint.observed_effect,
        "candidate_primary_metric": endpoint.candidate_primary_metric,
        "baseline_primary_metric": endpoint.baseline_primary_metric,
        "candidate_primary_standard_error": _sample_standard_error(
            endpoint.candidate_primary_bootstrap
        ),
        "baseline_primary_standard_error": _sample_standard_error(
            endpoint.baseline_primary_bootstrap
        ),
        "paired_effect_standard_error": _sample_standard_error(
            endpoint.paired_difference_bootstrap
        ),
        "candidate_primary_bootstrap_sha256": canonical_sha256(
            endpoint.candidate_primary_bootstrap
        ),
        "baseline_primary_bootstrap_sha256": canonical_sha256(
            endpoint.baseline_primary_bootstrap
        ),
        "paired_difference_bootstrap_sha256": canonical_sha256(
            endpoint.paired_difference_bootstrap
        ),
        "n_units": endpoint.n_units,
        "independent_unit_counts": dict(endpoint.independent_unit_counts),
        "unit_set_sha256": endpoint.unit_set_sha256,
        "n_rows": endpoint.n_rows,
        "row_set_sha256": endpoint.row_set_sha256,
        "strata": list(endpoint.strata),
        "seed": endpoint.seed,
        "n_resamples": endpoint.n_resamples,
    }


def _scientific_receipt_identity(
    *,
    candidate_execution_binding: Mapping[str, str],
    baseline_execution_binding: Mapping[str, str],
    candidate_campaign_binding: Mapping[str, str],
    baseline_campaign_binding: Mapping[str, str],
    candidate_prediction_binding: Mapping[str, Any],
    baseline_prediction_binding: Mapping[str, Any],
    outcome_binding: Mapping[str, Any],
    run: Mapping[str, Any],
    baseline_run: Mapping[str, Any],
    family_id: str,
    joined_table: ArtifactRef,
    endpoint_parameters: Mapping[str, Any],
    endpoint_summary: Mapping[str, Any],
    evaluator_sha256: str,
    task_contract_sha256: str,
    endpoint_roster_authority: Mapping[str, Any] | None,
) -> dict[str, Any]:
    metrics = {
        "absolute_primary_metric": endpoint_summary["candidate_primary_metric"],
        "candidate_minus_fixed_baseline_effect": endpoint_summary["observed_effect"],
    }
    standard_errors = {
        "absolute_primary_metric": endpoint_summary["candidate_primary_standard_error"],
        "candidate_minus_fixed_baseline_effect": endpoint_summary[
            "paired_effect_standard_error"
        ],
    }
    return {
        "schema_version": SCIENTIFIC_RUN_RECEIPT_SCHEMA_VERSION,
        "candidate_execution_receipt_binding": dict(candidate_execution_binding),
        "baseline_execution_receipt_binding": dict(baseline_execution_binding),
        "candidate_campaign_binding": dict(candidate_campaign_binding),
        "baseline_campaign_binding": dict(baseline_campaign_binding),
        "candidate_prediction_binding": dict(candidate_prediction_binding),
        "baseline_prediction_binding": dict(baseline_prediction_binding),
        "development_outcome_binding": dict(outcome_binding),
        "run_id": run["run_id"],
        "baseline_run_id": baseline_run["run_id"],
        "task_id": run["task_id"],
        "model_id": run["model_id"],
        "baseline_model_id": baseline_run["model_id"],
        "family_id": family_id,
        "adaptation_regime": run["adaptation_regime"],
        "seed": run["seed"],
        "fold": run["fold"],
        "endpoint_evaluator_id": endpoint_summary["evaluator_id"],
        "endpoint_evaluator_sha256": evaluator_sha256,
        "task_contract_sha256": task_contract_sha256,
        "endpoint_parameters": dict(endpoint_parameters),
        "endpoint_roster_authority": (
            dict(endpoint_roster_authority)
            if endpoint_roster_authority is not None
            else None
        ),
        "joined_endpoint_table": joined_table.to_dict(),
        "endpoint_summary": dict(endpoint_summary),
        "metrics": metrics,
        "standard_errors": standard_errors,
        "compute_cost_status": "unavailable_no_comparable_measured_basis",
        "scientific_status": "succeeded",
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
    }


def _scientific_sources(
    *,
    candidate_campaign_dir: str | Path,
    candidate_execution_attempt_dir: str | Path,
    candidate_prediction_bundle_path: str | Path,
    baseline_campaign_dir: str | Path,
    baseline_execution_attempt_dir: str | Path,
    baseline_prediction_bundle_path: str | Path,
    development_outcome_bundle_dir: str | Path,
) -> dict[str, Any]:
    _, candidate_plan, candidate_campaign = _candidate_plan_source(candidate_campaign_dir)
    _, baseline_plan, baseline_campaign = _candidate_plan_source(baseline_campaign_dir)
    _, candidate_execution, candidate_execution_binding = _verify_execution_attempt(
        candidate_execution_attempt_dir
    )
    _, baseline_execution, baseline_execution_binding = _verify_execution_attempt(
        baseline_execution_attempt_dir
    )
    run = _run_from_plan(
        candidate_plan,
        _sha256_identifier(candidate_execution.get("run_id"), "candidate execution run_id"),
    )
    baseline_run = _run_from_plan(
        baseline_plan,
        _sha256_identifier(baseline_execution.get("run_id"), "baseline execution run_id"),
    )
    for execution, campaign, frozen_run, label in (
        (candidate_execution, candidate_campaign, run, "candidate"),
        (baseline_execution, baseline_campaign, baseline_run, "baseline"),
    ):
        for field in ("plan_sha256", "candidate_manifest_sha256", "task_id", "model_id", "seed", "fold"):
            if field == "plan_sha256":
                expected = campaign["plan_sha256"]
            elif field == "candidate_manifest_sha256":
                expected = campaign["manifest_sha256"]
            else:
                expected = frozen_run[field]
            if execution.get(field) != expected:
                raise TournamentError(f"{label} execution receipt changed frozen {field}")
    if (
        run["task_id"] != baseline_run["task_id"]
        or run["seed"] != baseline_run["seed"]
        or run["fold"] != baseline_run["fold"]
    ):
        raise TournamentError("candidate and fixed baseline runs must align by task/seed/fold")
    if run["model_id"] == baseline_run["model_id"]:
        raise TournamentError(
            "scientific receipt candidate and baseline models must be distinct"
        )
    task_id = str(run["task_id"])
    if task_id not in _TASK_ENDPOINT_EVALUATORS:
        raise TournamentError(f"task {task_id!r} has no registered scientific evaluator")
    candidate_task_records = [
        item
        for item in candidate_plan.get("task_dispositions", ())
        if isinstance(item, Mapping) and item.get("task_id") == task_id
    ]
    baseline_task_records = [
        item
        for item in baseline_plan.get("task_dispositions", ())
        if isinstance(item, Mapping) and item.get("task_id") == task_id
    ]
    if len(candidate_task_records) != 1 or len(baseline_task_records) != 1:
        raise TournamentError(
            f"scientific receipt requires one TaskSpec disposition for {task_id}"
        )
    candidate_roster = candidate_task_records[0].get("baseline_model_ids")
    baseline_roster = baseline_task_records[0].get("baseline_model_ids")
    if (
        not isinstance(candidate_roster, list)
        or not isinstance(baseline_roster, list)
        or baseline_run["model_id"] not in candidate_roster
        or baseline_run["model_id"] not in baseline_roster
    ):
        raise TournamentError(
            "scientific receipt baseline is absent from the frozen TaskSpec roster"
        )
    if task_id == VARIANT_TASK:
        for label, frozen_run in (("candidate", run), ("baseline", baseline_run)):
            metadata = frozen_run.get("metadata")
            if (
                not isinstance(metadata, Mapping)
                or metadata.get("primary_endpoint_scoring_allowed") is not True
            ):
                raise TournamentError(
                    "variant primary ScientificRunReceipt forbids "
                    f"{label} run {frozen_run['run_id']}: "
                    "primary_endpoint_scoring_allowed must be true"
                )
    candidate_endpoint_contract = _endpoint_contract_from_plan(candidate_plan, task_id)
    baseline_endpoint_contract = _endpoint_contract_from_plan(baseline_plan, task_id)
    task_class_roster = (
        candidate_endpoint_contract[1].get("class_roster")
        if task_id == CELL_TASK
        else None
    )
    candidate_bundle, candidate_rows, candidate_prediction = _prediction_source(
        candidate_prediction_bundle_path, run=run, class_roster=task_class_roster
    )
    baseline_bundle, baseline_rows, baseline_prediction = _prediction_source(
        baseline_prediction_bundle_path,
        run=baseline_run,
        class_roster=task_class_roster,
    )
    outcome_root = _safe_resolve(
        development_outcome_bundle_dir, "development outcome bundle"
    )
    outcome = verify_development_outcome_bundle(outcome_root)
    outcome_binding = {
        "path": outcome_root.as_posix(),
        "manifest_sha256": sha256_file(outcome_root / "ARTIFACTS.json"),
        "document_sha256": sha256_file(outcome_root / "development_outcome_bundle.json"),
        "outcome_bundle_id": outcome["outcome_bundle_id"],
    }
    join_identity = (
        candidate_bundle.source_join_key_sha256,
        baseline_bundle.source_join_key_sha256,
        outcome["source_join_key_sha256"],
    )
    if len(set(join_identity)) != 1:
        raise TournamentError("candidate, baseline, and outcome source-join identities differ")
    if candidate_bundle.dataset_ids != baseline_bundle.dataset_ids or list(candidate_bundle.dataset_ids) != outcome["dataset_ids"]:
        raise TournamentError("candidate, baseline, and outcome dataset order differs")
    candidate_dataset_bindings = _run_dataset_registry_bindings(candidate_plan, run)
    baseline_dataset_bindings = _run_dataset_registry_bindings(
        baseline_plan, baseline_run
    )
    if (
        candidate_dataset_bindings != baseline_dataset_bindings
        or candidate_dataset_bindings
        != dict(sorted(outcome["dataset_registry_sha256s"].items()))
    ):
        raise TournamentError(
            "candidate, baseline, and outcome dataset registry bindings differ"
        )
    outcome_rows = _read_exact_tsv(
        outcome_root / outcome["standardized_table"]["path"],
        _OUTCOME_TABLE_FIELDS[task_id],
        "development outcome table",
    )
    joined_rows = _join_development_rows(
        task_id=task_id,
        candidate_rows=candidate_rows,
        baseline_rows=baseline_rows,
        outcome_rows=outcome_rows,
    )
    if candidate_endpoint_contract != baseline_endpoint_contract:
        raise TournamentError(
            "candidate and baseline campaigns use different frozen TaskSpec evaluators"
        )
    (
        evaluator_id,
        endpoint_parameters,
        task_contract_sha256,
        endpoint_roster_authority,
    ) = candidate_endpoint_contract
    disposition = _model_disposition_for_plan(candidate_plan, str(run["model_id"]))
    return {
        "candidate_execution_binding": candidate_execution_binding,
        "baseline_execution_binding": baseline_execution_binding,
        "candidate_campaign_binding": candidate_campaign,
        "baseline_campaign_binding": baseline_campaign,
        "candidate_prediction_binding": candidate_prediction,
        "baseline_prediction_binding": baseline_prediction,
        "outcome_binding": outcome_binding,
        "run": run,
        "baseline_run": baseline_run,
        "family_id": _nonempty_identifier(disposition["family_id"], "family_id"),
        "task_id": task_id,
        "joined_rows": joined_rows,
        "endpoint_evaluator_id": evaluator_id,
        "endpoint_parameters": endpoint_parameters,
        "task_contract_sha256": task_contract_sha256,
        "endpoint_roster_authority": endpoint_roster_authority,
    }


def freeze_scientific_run_receipt(
    *,
    candidate_campaign_dir: str | Path,
    candidate_execution_attempt_dir: str | Path,
    candidate_prediction_bundle_path: str | Path,
    baseline_campaign_dir: str | Path,
    baseline_execution_attempt_dir: str | Path,
    baseline_prediction_bundle_path: str | Path,
    development_outcome_bundle_dir: str | Path,
    output_root: str | Path,
) -> Path:
    """Join verified predictions/outcomes and freeze task-native development metrics."""

    sources = _scientific_sources(
        candidate_campaign_dir=candidate_campaign_dir,
        candidate_execution_attempt_dir=candidate_execution_attempt_dir,
        candidate_prediction_bundle_path=candidate_prediction_bundle_path,
        baseline_campaign_dir=baseline_campaign_dir,
        baseline_execution_attempt_dir=baseline_execution_attempt_dir,
        baseline_prediction_bundle_path=baseline_prediction_bundle_path,
        development_outcome_bundle_dir=development_outcome_bundle_dir,
    )
    output = _safe_resolve(output_root, "ScientificRunReceipt output root")
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".scientific-run.", dir=output))
    try:
        joined_path = staging / "joined_endpoint.tsv"
        joined_fields = _JOINED_ENDPOINT_FIELDS[sources["task_id"]]
        write_text_exclusive(
            joined_path,
            _tsv_text(joined_fields, sources["joined_rows"]),
            mode=0o440,
        )
        endpoint = recompute_endpoint_bootstrap(
            evaluator_id=sources["endpoint_evaluator_id"],
            table_path=joined_path,
            n_resamples=_SELECTION_BOOTSTRAP_RESAMPLES,
            seed=_SELECTION_BOOTSTRAP_SEED,
            parameters=sources["endpoint_parameters"],
        )
        joined_ref = ArtifactRef.from_path(
            joined_path,
            relative_to=staging,
            media_type="text/tab-separated-values",
            role=f"joined_development_endpoint:{sources['task_id']}",
        )
        identity = _scientific_receipt_identity(
            candidate_execution_binding=sources["candidate_execution_binding"],
            baseline_execution_binding=sources["baseline_execution_binding"],
            candidate_campaign_binding=sources["candidate_campaign_binding"],
            baseline_campaign_binding=sources["baseline_campaign_binding"],
            candidate_prediction_binding=sources["candidate_prediction_binding"],
            baseline_prediction_binding=sources["baseline_prediction_binding"],
            outcome_binding=sources["outcome_binding"],
            run=sources["run"],
            baseline_run=sources["baseline_run"],
            family_id=sources["family_id"],
            joined_table=joined_ref,
            endpoint_parameters=sources["endpoint_parameters"],
            endpoint_summary=_endpoint_summary(endpoint),
            evaluator_sha256=sha256_file(Path(recompute_endpoint_bootstrap.__code__.co_filename)),
            task_contract_sha256=sources["task_contract_sha256"],
            endpoint_roster_authority=sources["endpoint_roster_authority"],
        )
        payload = {"scientific_receipt_id": canonical_sha256(identity), **identity}
        write_json_exclusive(
            staging / "scientific_run_receipt.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "scientific_run_receipt",
                "scientific_receipt_id": payload["scientific_receipt_id"],
                "run_id": sources["run"]["run_id"],
            },
        )
        target = output / f"scientific-run--{payload['scientific_receipt_id']}"
        if target.exists():
            raise TournamentError(f"scientific run receipt already exists: {target}")
        publish_directory_noreplace(staging, target)
    except (EndpointPowerError, ContractError, OSError, ValueError) as exc:
        shutil.rmtree(staging, ignore_errors=True)
        raise TournamentError(f"scientific endpoint evaluation failed: {exc}") from exc
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_scientific_run_receipt(target)
    return target


def verify_scientific_run_receipt(path: str | Path) -> dict[str, Any]:
    root = _safe_resolve(path, "ScientificRunReceipt")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "scientific_run_receipt.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid frozen ScientificRunReceipt: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != SCIENTIFIC_RUN_RECEIPT_SCHEMA_VERSION:
        raise TournamentError("unsupported ScientificRunReceipt schema")
    try:
        sources = _scientific_sources(
            candidate_campaign_dir=payload["candidate_campaign_binding"]["path"],
            candidate_execution_attempt_dir=payload[
                "candidate_execution_receipt_binding"
            ]["path"],
            candidate_prediction_bundle_path=payload[
                "candidate_prediction_binding"
            ]["path"],
            baseline_campaign_dir=payload["baseline_campaign_binding"]["path"],
            baseline_execution_attempt_dir=payload[
                "baseline_execution_receipt_binding"
            ]["path"],
            baseline_prediction_bundle_path=payload[
                "baseline_prediction_binding"
            ]["path"],
            development_outcome_bundle_dir=payload[
                "development_outcome_binding"
            ]["path"],
        )
        joined_ref = ArtifactRef.from_dict(payload.get("joined_endpoint_table"))
    except (KeyError, TypeError, ContractError) as exc:
        raise TournamentError(f"ScientificRunReceipt has invalid bindings: {exc}") from exc
    if (
        joined_ref.media_type != "text/tab-separated-values"
        or joined_ref.role
        != f"joined_development_endpoint:{sources['task_id']}"
    ):
        raise TournamentError("ScientificRunReceipt joined endpoint ArtifactRef is invalid")
    joined_path = joined_ref.validate(root, require_relative=True)
    expected_rows = _read_exact_tsv(
        joined_path,
        _JOINED_ENDPOINT_FIELDS[sources["task_id"]],
        "joined scientific endpoint",
    )
    if expected_rows != sources["joined_rows"]:
        raise TournamentError("ScientificRunReceipt joined table does not rederive")
    try:
        endpoint = recompute_endpoint_bootstrap(
            evaluator_id=sources["endpoint_evaluator_id"],
            table_path=joined_path,
            n_resamples=_SELECTION_BOOTSTRAP_RESAMPLES,
            seed=_SELECTION_BOOTSTRAP_SEED,
            parameters=sources["endpoint_parameters"],
        )
    except EndpointPowerError as exc:
        raise TournamentError(f"ScientificRunReceipt endpoint does not rederive: {exc}") from exc
    identity = _scientific_receipt_identity(
        candidate_execution_binding=sources["candidate_execution_binding"],
        baseline_execution_binding=sources["baseline_execution_binding"],
        candidate_campaign_binding=sources["candidate_campaign_binding"],
        baseline_campaign_binding=sources["baseline_campaign_binding"],
        candidate_prediction_binding=sources["candidate_prediction_binding"],
        baseline_prediction_binding=sources["baseline_prediction_binding"],
        outcome_binding=sources["outcome_binding"],
        run=sources["run"],
        baseline_run=sources["baseline_run"],
        family_id=sources["family_id"],
        joined_table=joined_ref,
        endpoint_parameters=sources["endpoint_parameters"],
        endpoint_summary=_endpoint_summary(endpoint),
        evaluator_sha256=sha256_file(Path(recompute_endpoint_bootstrap.__code__.co_filename)),
        task_contract_sha256=sources["task_contract_sha256"],
        endpoint_roster_authority=sources["endpoint_roster_authority"],
    )
    claimed = _sha256_identifier(payload.get("scientific_receipt_id"), "scientific_receipt_id")
    if dict(payload) != {"scientific_receipt_id": claimed, **identity} or canonical_sha256(identity) != claimed:
        raise TournamentError("ScientificRunReceipt identity does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": "scientific_run_receipt",
        "scientific_receipt_id": claimed,
        "run_id": sources["run"]["run_id"],
    }:
        raise TournamentError("ScientificRunReceipt manifest metadata mismatch")
    return dict(payload)


def _variant_secondary_run_receipt_identity(
    sources: Mapping[str, Any],
) -> dict[str, Any]:
    run = sources["run"]
    bundle = sources["bundle"]
    return {
        "schema_version": VARIANT_SECONDARY_RUN_RECEIPT_SCHEMA_VERSION,
        "candidate_campaign_binding": dict(sources["campaign_binding"]),
        "execution_receipt_binding": dict(sources["execution_binding"]),
        "prediction_binding": dict(sources["prediction_binding"]),
        "run_id": run["run_id"],
        "task_id": VARIANT_TASK,
        "model_id": run["model_id"],
        "family_id": sources["family_id"],
        "adaptation_regime": run["adaptation_regime"],
        "seed": run["seed"],
        "fold": run["fold"],
        "dataset_ids": list(run["dataset_ids"]),
        "dataset_registry_sha256s": dict(sources["dataset_registry_sha256s"]),
        "immutable_inputs": dict(sources["immutable_inputs"]),
        "variant_capability": dict(sources["capability"]),
        "variant_capability_sha256": sources["capability_sha256"],
        "variant_capability_registry_sha256": sources[
            "capability_registry_sha256"
        ],
        "endpoint_id": _VARIANT_SECONDARY_ENDPOINT_ID,
        "evaluator_id": _VARIANT_SECONDARY_EVALUATOR_ID,
        "score_transform_id": _VARIANT_SECONDARY_SCORE_TRANSFORM_ID,
        "prediction_bundle_id": bundle.bundle_id,
        "prediction_table_schema_sha256": bundle.table_schema_sha256,
        "source_join_key_sha256": bundle.source_join_key_sha256,
        "n_predictions": bundle.n_predictions,
        "primary_endpoint_scoring_allowed": False,
        "outcomes_read": False,
        "metrics_computed": False,
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
        "auxiliary_status": "succeeded_non_scoring",
    }


def freeze_variant_secondary_run_receipt(
    *,
    candidate_campaign_dir: str | Path,
    execution_attempt_dir: str | Path,
    prediction_bundle_path: str | Path,
    output_root: str | Path,
) -> Path:
    """Freeze a recursively verified link-score run without reading outcomes."""

    sources = _variant_secondary_run_sources(
        candidate_campaign_dir=candidate_campaign_dir,
        execution_attempt_dir=execution_attempt_dir,
        prediction_bundle_path=prediction_bundle_path,
    )
    identity = _variant_secondary_run_receipt_identity(sources)
    payload = {
        "secondary_run_receipt_id": canonical_sha256(identity),
        **identity,
    }
    output = _safe_resolve(output_root, "variant secondary receipt output root")
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"variant-secondary-run--{payload['secondary_run_receipt_id']}"
    if target.exists():
        raise TournamentError(f"variant secondary run receipt already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".variant-secondary-run.", dir=output))
    try:
        write_json_exclusive(
            staging / "variant_secondary_run_receipt.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "variant_secondary_run_receipt",
                "secondary_run_receipt_id": payload["secondary_run_receipt_id"],
                "run_id": identity["run_id"],
            },
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_variant_secondary_run_receipt(target)
    return target


def verify_variant_secondary_run_receipt(path: str | Path) -> dict[str, Any]:
    """Recursively verify one read-only non-scoring variant comparator receipt."""

    root = _safe_resolve(path, "variant secondary run receipt")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "variant_secondary_run_receipt.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid variant secondary run receipt: {exc}") from exc
    if (
        not isinstance(payload, Mapping)
        or payload.get("schema_version")
        != VARIANT_SECONDARY_RUN_RECEIPT_SCHEMA_VERSION
    ):
        raise TournamentError("unsupported variant secondary run receipt schema")
    expected_fields = {
        "secondary_run_receipt_id",
        "schema_version",
        "candidate_campaign_binding",
        "execution_receipt_binding",
        "prediction_binding",
        "run_id",
        "task_id",
        "model_id",
        "family_id",
        "adaptation_regime",
        "seed",
        "fold",
        "dataset_ids",
        "dataset_registry_sha256s",
        "immutable_inputs",
        "variant_capability",
        "variant_capability_sha256",
        "variant_capability_registry_sha256",
        "endpoint_id",
        "evaluator_id",
        "score_transform_id",
        "prediction_bundle_id",
        "prediction_table_schema_sha256",
        "source_join_key_sha256",
        "n_predictions",
        "primary_endpoint_scoring_allowed",
        "outcomes_read",
        "metrics_computed",
        "created_before_outcome_unblind",
        "sealed_results_used",
        "auxiliary_status",
    }
    if set(payload) != expected_fields:
        raise TournamentError("variant secondary run receipt has the wrong exact schema")
    campaign = payload.get("candidate_campaign_binding")
    execution = payload.get("execution_receipt_binding")
    prediction = payload.get("prediction_binding")
    if not all(isinstance(item, Mapping) for item in (campaign, execution, prediction)):
        raise TournamentError("variant secondary run receipt bindings are invalid")
    try:
        sources = _variant_secondary_run_sources(
            candidate_campaign_dir=campaign["path"],
            execution_attempt_dir=execution["path"],
            prediction_bundle_path=prediction["path"],
        )
    except KeyError as exc:
        raise TournamentError(
            f"variant secondary run receipt binding is incomplete: {exc}"
        ) from exc
    identity = _variant_secondary_run_receipt_identity(sources)
    claimed = _sha256_identifier(
        payload.get("secondary_run_receipt_id"), "secondary_run_receipt_id"
    )
    if dict(payload) != {"secondary_run_receipt_id": claimed, **identity}:
        raise TournamentError("variant secondary run receipt does not rederive")
    if canonical_sha256(identity) != claimed:
        raise TournamentError("variant secondary run receipt identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "variant_secondary_run_receipt",
        "secondary_run_receipt_id": claimed,
        "run_id": identity["run_id"],
    }:
        raise TournamentError("variant secondary run receipt manifest metadata mismatch")
    return dict(payload)


def _variant_secondary_receipt_binding(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    root = _safe_resolve(path, "variant secondary run receipt")
    payload = verify_variant_secondary_run_receipt(root)
    binding, _ = _frozen_document_binding(
        root,
        filename="variant_secondary_run_receipt.json",
        artifact_class="variant_secondary_run_receipt",
    )
    return (
        {
            **binding,
            "secondary_run_receipt_id": payload["secondary_run_receipt_id"],
            "run_id": payload["run_id"],
            "model_id": payload["model_id"],
        },
        payload,
    )


def _receipt_binding(path: str | Path) -> tuple[dict[str, Any], dict[str, Any]]:
    root = _safe_resolve(path, "ScientificRunReceipt")
    payload = verify_scientific_run_receipt(root)
    binding, _ = _frozen_document_binding(
        root,
        filename="scientific_run_receipt.json",
        artifact_class="scientific_run_receipt",
    )
    return (
        {
            **binding,
            "scientific_receipt_id": payload["scientific_receipt_id"],
            "run_id": payload["run_id"],
            "baseline_run_id": payload["baseline_run_id"],
        },
        payload,
    )


def _campaign_from_receipt_binding(
    binding: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    root, plan, observed = _candidate_plan_source(binding.get("path"))
    expected = {
        "path": root.as_posix(),
        "manifest_sha256": _sha256_identifier(
            binding.get("manifest_sha256"), "campaign manifest_sha256"
        ),
        "document_sha256": _sha256_identifier(
            binding.get("document_sha256"), "campaign document_sha256"
        ),
        "plan_sha256": _sha256_identifier(
            binding.get("plan_sha256"), "campaign plan_sha256"
        ),
    }
    if observed != expected:
        raise TournamentError("ScientificRunReceipt campaign binding changed")
    return plan, observed


def _merge_exact_record(
    records: dict[str, dict[str, Any]],
    *,
    identifier: str,
    value: Mapping[str, Any],
    label: str,
) -> None:
    normalized = dict(value)
    prior = records.get(identifier)
    if prior is not None and prior != normalized:
        raise TournamentError(f"conflicting {label} across frozen campaigns: {identifier}")
    records[identifier] = normalized


def _aggregate_standard_error(
    values: Sequence[float], within_run_errors: Sequence[float]
) -> float:
    if len(values) != len(within_run_errors) or not values:
        raise TournamentError("model metric aggregation requires aligned observations")
    between_component = 0.0
    if len(values) > 1:
        mean = fsum(values) / len(values)
        between_component = (
            fsum((value - mean) ** 2 for value in values)
            / (len(values) - 1)
        )
    within_component = fsum(
        error * error for error in within_run_errors
    ) / len(within_run_errors)
    return sqrt(between_component + within_component)


def _scientific_comparison_identity(
    *,
    run: Mapping[str, Any],
    baseline_run: Mapping[str, Any],
    receipt: Mapping[str, Any],
) -> dict[str, Any]:
    endpoint = receipt["endpoint_summary"]
    identity = {
        "task_id": run["task_id"],
        # The fixed baseline is part of the comparison.  For
        # rna_conditioned_atac the primary metric is a relative deviance
        # reduction against this baseline, so two candidates scored against
        # different baselines are not comparable quantities at all.
        "fixed_baseline_model_id": receipt["baseline_model_id"],
        "fixed_baseline_candidate_configuration_sha256": canonical_sha256(
            _candidate_configuration_identity(baseline_run)
        ),
        "split_id": run["split_id"],
        "fold": run["fold"],
        "dataset_ids": list(run["dataset_ids"]),
        "immutable_inputs": dict(run["immutable_inputs"]),
        "endpoint_evaluator_id": receipt["endpoint_evaluator_id"],
        "endpoint_evaluator_sha256": receipt["endpoint_evaluator_sha256"],
        "endpoint_parameters_sha256": canonical_sha256(
            receipt["endpoint_parameters"]
        ),
        "task_contract_sha256": receipt["task_contract_sha256"],
        "development_outcome_bundle_id": receipt["development_outcome_binding"][
            "outcome_bundle_id"
        ],
        "development_outcome_manifest_sha256": receipt[
            "development_outcome_binding"
        ]["manifest_sha256"],
        "n_units": endpoint["n_units"],
        "unit_set_sha256": endpoint["unit_set_sha256"],
        "n_rows": endpoint["n_rows"],
        "row_set_sha256": endpoint["row_set_sha256"],
        "strata": list(endpoint["strata"]),
        "bootstrap_seed": endpoint["seed"],
        "bootstrap_replicates": endpoint["n_resamples"],
    }
    return identity


def _candidate_configuration_identity(run: Mapping[str, Any]) -> dict[str, Any]:
    raw = dict(run)
    claimed_run_id = raw.pop("run_id", None)
    try:
        spec = RunSpec.from_dict(raw)
    except ContractError as exc:
        raise TournamentError(
            f"candidate run does not satisfy the complete RunSpec: {exc}"
        ) from exc
    if claimed_run_id is not None and spec.run_id != claimed_run_id:
        raise TournamentError("candidate run identity does not rederive")
    identity = spec.identity_payload
    identity.pop("seed")
    return identity


def _selection_artifact_hashes(
    execution_binding: Mapping[str, Any],
) -> dict[str, str]:
    _, execution, observed_binding = _verify_execution_attempt(
        execution_binding.get("path")
    )
    if any(
        execution_binding.get(field) != observed_binding.get(field)
        for field in ("path", "manifest_sha256", "document_sha256")
    ):
        raise TournamentError("scientific execution binding changed")
    raw_actions = execution.get("action_receipts")
    if not isinstance(raw_actions, list):
        raise TournamentError("scientific execution lacks action receipts")
    actions: dict[str, Mapping[str, Any]] = {}
    for record in raw_actions:
        if not isinstance(record, Mapping):
            raise TournamentError("scientific action receipt must be an object")
        action = str(record.get("action", ""))
        if action in actions:
            raise TournamentError(f"scientific execution repeats action {action}")
        actions[action] = record
    if "prepare" not in actions or "fit" not in actions:
        raise TournamentError(
            "model selection requires frozen prepare and fit action artifacts"
        )
    prepare_manifest = _sha256_identifier(
        actions["prepare"].get("output_manifest_sha256"),
        "prepare output manifest SHA-256",
    )
    fit_manifest = _sha256_identifier(
        actions["fit"].get("output_manifest_sha256"),
        "fit output manifest SHA-256",
    )
    # One fit action, one read-only output manifest.  The three fit-state roles
    # bind that same bundle by policy; this makes the claim explicit rather
    # than leaving three identical-looking fields to be read as independent.
    return {
        "fit_state_artifact_policy": FIT_STATE_ARTIFACT_POLICY,
        "fit_state_bundle_sha256": fit_manifest,
        "role_hashes": {
            "checkpoint_sha256": fit_manifest,
            "preprocessing_sha256": prepare_manifest,
            "task_head_sha256": fit_manifest,
            "calibration_sha256": fit_manifest,
        },
    }


def _model_candidates_from_scientific_runs(
    scientific_runs: Sequence[Mapping[str, Any]],
    dispositions: Mapping[str, Mapping[str, Any]],
    receipt_index: Mapping[str, tuple[Path, Mapping[str, Any]]],
    expected_seed_count: int,
) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str, str], list[Mapping[str, Any]]] = {}
    for record in scientific_runs:
        key = (
            str(record["task_id"]),
            str(record["model_id"]),
            str(record["adaptation_regime"]),
            _sha256_identifier(
                record.get("candidate_configuration_sha256"),
                "candidate_configuration_sha256",
            ),
        )
        groups.setdefault(key, []).append(record)
    candidates: list[dict[str, Any]] = []
    for (task_id, model_id, regime, configuration_sha256), records in sorted(
        groups.items()
    ):
        disposition = dispositions.get(model_id)
        if disposition is None:
            raise TournamentError(f"scientific model lacks frozen disposition: {model_id}")
        family_id = _nonempty_identifier(disposition.get("family_id"), "family_id")
        metric_names = set(records[0]["metrics"])
        if any(set(record["metrics"]) != metric_names for record in records):
            raise TournamentError(f"scientific metrics differ across runs for {model_id}")
        if any(set(record["standard_errors"]) != metric_names for record in records):
            raise TournamentError(f"scientific standard errors are incomplete for {model_id}")
        candidate_configuration = records[0].get("candidate_configuration")
        if (
            not isinstance(candidate_configuration, Mapping)
            or canonical_sha256(candidate_configuration) != configuration_sha256
            or any(
                record.get("candidate_configuration") != candidate_configuration
                for record in records
            )
        ):
            raise TournamentError(
                f"scientific candidate configuration does not rederive for {model_id}"
            )
        comparison_hashes = {
            _sha256_identifier(
                record.get("comparison_identity_sha256"),
                f"{model_id}.comparison_identity_sha256",
            )
            for record in records
        }
        if len(comparison_hashes) != 1:
            raise TournamentError(
                f"scientific runs for {model_id} mix folds, evaluators, inputs, "
                "or row/unit sets"
            )
        comparison_identity = records[0].get("comparison_identity")
        if (
            not isinstance(comparison_identity, Mapping)
            or canonical_sha256(comparison_identity) != next(iter(comparison_hashes))
            or any(record.get("comparison_identity") != comparison_identity for record in records)
        ):
            raise TournamentError(
                f"scientific comparison identity does not rederive for {model_id}"
            )
        # Enforce the composite fit-state policy literally, per record, before
        # any digest is folded into an ensemble hash.
        declared_policies = {
            record.get("fit_state_artifact_policy") for record in records
        }
        if declared_policies != {FIT_STATE_ARTIFACT_POLICY} or not (
            declared_policies <= SUPPORTED_FIT_STATE_ARTIFACT_POLICIES
        ):
            raise TournamentError(
                "scientific run declares an unsupported fit-state artifact "
                f"policy for {model_id}"
            )
        for record in records:
            artifacts = record.get("selection_artifact_hashes")
            bundle = record.get("fit_state_bundle_sha256")
            if not isinstance(artifacts, Mapping):
                raise TournamentError(
                    f"scientific run lacks selection artifacts for {model_id}"
                )
            if not (
                artifacts.get("checkpoint_sha256")
                == artifacts.get("task_head_sha256")
                == artifacts.get("calibration_sha256")
                == bundle
            ):
                raise TournamentError(
                    f"{model_id} fit-state roles do not bind one composite "
                    "immutable fit bundle"
                )
            if artifacts.get("preprocessing_sha256") == bundle:
                raise TournamentError(
                    f"{model_id} preprocessing must come from the separate "
                    "prepare manifest"
                )
        artifact_fields = (
            "checkpoint_sha256",
            "preprocessing_sha256",
            "task_head_sha256",
            "calibration_sha256",
        )
        selection_artifact_ensemble_sha256s: dict[str, str] = {}
        for field_name in artifact_fields:
            per_run: dict[str, str] = {}
            for record in records:
                artifacts = record.get("selection_artifact_hashes")
                if not isinstance(artifacts, Mapping):
                    raise TournamentError(
                        f"scientific run lacks selection artifacts for {model_id}"
                    )
                per_run[str(record["run_id"])] = _sha256_identifier(
                    artifacts.get(field_name),
                    f"{model_id}.{field_name}",
                )
            selection_artifact_ensemble_sha256s[field_name] = canonical_sha256(
                dict(sorted(per_run.items()))
            )
        # Per-seed values are retained only as a stability diagnostic; ranking
        # uses the endpoint of the five-seed ensemble, which is the output file
        # held-back inference actually deploys.
        per_seed_metrics: dict[str, float] = {}
        per_seed_standard_errors: dict[str, float] = {}
        for metric in sorted(metric_names):
            values = [
                _finite_number(record["metrics"][metric], f"{model_id}.{metric}")
                for record in records
            ]
            errors = [
                _finite_number(
                    record["standard_errors"][metric],
                    f"{model_id}.{metric}.standard_error",
                )
                for record in records
            ]
            if any(error < 0.0 for error in errors):
                raise TournamentError("scientific standard errors cannot be negative")
            per_seed_metrics[metric] = fsum(values) / len(values)
            per_seed_standard_errors[metric] = _aggregate_standard_error(
                values, errors
            )
        metrics = dict(per_seed_metrics)
        standard_errors = dict(per_seed_standard_errors)
        ensemble_record: dict[str, Any] | None = None
        if len(records) == expected_seed_count == _ENSEMBLE_SEED_COUNT:
            per_seed_sources: list[tuple[int, Path, Mapping[str, Any]]] = []
            roles: set[str] = set()
            for record in records:
                bound = receipt_index.get(str(record["run_id"]))
                if bound is None:
                    per_seed_sources = []
                    break
                per_seed_sources.append(
                    (int(record["seed"]), bound[0], bound[1])
                )
                roles.add(bound[2])
            if per_seed_sources and len(roles) != 1:
                raise TournamentError(
                    f"{model_id} mixes candidate and baseline roles across seeds"
                )
            if per_seed_sources:
                role = next(iter(roles))
                receipt = per_seed_sources[0][2]
                endpoint_parameters = dict(receipt["endpoint_parameters"])
                ensemble_rows, alignment = _development_ensemble_rows(
                    task_id=task_id,
                    per_seed=per_seed_sources,
                    endpoint_parameters=endpoint_parameters,
                )
                staging = Path(
                    tempfile.mkdtemp(prefix=".candidate-ensemble.")
                )
                try:
                    table_path = staging / "ensemble_endpoint.tsv"
                    write_text_exclusive(
                        table_path,
                        _tsv_text(
                            _JOINED_ENDPOINT_FIELDS[task_id], ensemble_rows
                        ),
                        mode=0o440,
                    )
                    try:
                        ensemble = recompute_endpoint_bootstrap(
                            evaluator_id=str(receipt["endpoint_evaluator_id"]),
                            table_path=table_path,
                            n_resamples=_SELECTION_BOOTSTRAP_RESAMPLES,
                            seed=_SELECTION_BOOTSTRAP_SEED,
                            parameters=endpoint_parameters,
                        )
                    except EndpointPowerError as exc:
                        raise TournamentError(
                            "candidate development ensemble failed to "
                            f"recompute for {model_id}: {exc}"
                        ) from exc
                finally:
                    shutil.rmtree(staging, ignore_errors=True)
                if role == "candidate":
                    ensemble_metric = ensemble.candidate_primary_metric
                    ensemble_bootstrap = ensemble.candidate_primary_bootstrap
                else:
                    ensemble_metric = ensemble.baseline_primary_metric
                    ensemble_bootstrap = ensemble.baseline_primary_bootstrap
                metrics["absolute_primary_metric"] = ensemble_metric
                standard_errors["absolute_primary_metric"] = (
                    _sample_standard_error(ensemble_bootstrap)
                )
                ensemble_record = {
                    "development_ensemble_policy_id": (
                        DEVELOPMENT_ENSEMBLE_POLICY_ID
                    ),
                    "ensemble_role": role,
                    "ensemble_alignment_sha256": alignment[
                        "ensemble_alignment_sha256"
                    ],
                    "ensemble_rows_sha256": alignment["ensemble_rows_sha256"],
                    "ensemble_evaluator_id": ensemble.evaluator_id,
                    "ensemble_bootstrap_seed": ensemble.seed,
                    "ensemble_n_resamples": ensemble.n_resamples,
                    "per_seed_absolute_primary_metric": per_seed_metrics.get(
                        "absolute_primary_metric"
                    ),
                    "per_seed_absolute_primary_standard_error": (
                        per_seed_standard_errors.get("absolute_primary_metric")
                    ),
                    "used_for_ranking": False,
                }
        identity = {
            "task_id": task_id,
            "model_id": model_id,
            "adaptation_regime": regime,
            "candidate_configuration_sha256": configuration_sha256,
        }
        candidates.append(
            {
                "candidate_id": canonical_sha256(identity),
                "task_id": task_id,
                "model_id": model_id,
                "adaptation_regime": regime,
                "candidate_configuration": dict(
                    candidate_configuration
                ),
                "candidate_configuration_sha256": configuration_sha256,
                "family": family_id,
                "eligible": (
                    disposition.get("disposition") != "blocked"
                    and model_id not in NON_SELECTABLE_NEGATIVE_CONTROL_IDS
                ),
                "metrics": metrics,
                "standard_errors": standard_errors,
                "standard_error_policy": CANDIDATE_STANDARD_ERROR_POLICY,
                "fit_state_artifact_policy": FIT_STATE_ARTIFACT_POLICY,
                "fit_state_bundle_ensemble_sha256": (
                    selection_artifact_ensemble_sha256s["checkpoint_sha256"]
                ),
                "metric_basis": (
                    DEVELOPMENT_ENSEMBLE_POLICY_ID
                    if ensemble_record is not None
                    else "per_seed_mean_not_champion_eligible_v1"
                ),
                "ensemble_alignment_sha256": (
                    ensemble_record["ensemble_alignment_sha256"]
                    if ensemble_record is not None
                    else None
                ),
                "ensemble_rows_sha256": (
                    ensemble_record["ensemble_rows_sha256"]
                    if ensemble_record is not None
                    else None
                ),
                "seed_stability": ensemble_record,
                "comparison_identity": dict(comparison_identity),
                "comparison_identity_sha256": next(iter(comparison_hashes)),
                "selection_artifact_ensemble_sha256s": (
                    selection_artifact_ensemble_sha256s
                ),
                "complexity_status": "unavailable_no_comparable_measured_basis",
                "run_ids": sorted(str(record["run_id"]) for record in records),
                "seed_count": len(records),
            }
        )
    return candidates


def _selection_ledger_identity(
    *,
    ledger_authority: str,
    campaign_universe_binding: Mapping[str, Any],
    reviewed_universe_id: str,
    campaign_universe_kind: str,
    campaign_universe_sha256_value: str,
    campaign_universe_wave_ids: Sequence[str],
    expected_seed_count: int,
    selection_universe_sha256: str,
    core_registry_contract_sha256: str,
    source_campaign_bindings: Sequence[Mapping[str, Any]],
    scientific_receipt_bindings: Sequence[Mapping[str, Any]],
    variant_secondary_receipt_bindings: Sequence[Mapping[str, Any]],
    terminal_run_bindings: Sequence[Mapping[str, Any]],
    selection_universe: Sequence[Mapping[str, Any]],
    task_seed_sets: Sequence[Mapping[str, Any]],
    dataset_locks: Mapping[str, Any],
    model_dispositions: Sequence[Mapping[str, Any]],
    task_dispositions: Sequence[Mapping[str, Any]],
    runs: Sequence[Mapping[str, Any]],
    scientific_runs: Sequence[Mapping[str, Any]],
    variant_secondary_runs: Sequence[Mapping[str, Any]],
    model_candidates: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "schema_version": SELECTION_CANDIDATE_LEDGER_SCHEMA_VERSION,
        "ledger_authority": ledger_authority,
        "lock_construction_allowed": (
            ledger_authority == LEDGER_AUTHORITY_FINALIST
        ),
        "champion_selection_allowed": (
            ledger_authority == LEDGER_AUTHORITY_FINALIST
        ),
        "campaign_universe_binding": dict(campaign_universe_binding),
        "reviewed_universe_id": reviewed_universe_id,
        "campaign_universe_kind": campaign_universe_kind,
        "campaign_universe_sha256": campaign_universe_sha256_value,
        "campaign_universe_wave_ids": list(campaign_universe_wave_ids),
        "expected_seed_count": expected_seed_count,
        "selection_universe_sha256": selection_universe_sha256,
        "registry_binding_kind": "core_registry_contract_sha256",
        "core_registry_contract_sha256": core_registry_contract_sha256,
        # Retained under the plan-shaped name that selection.py consumes; it is
        # an exact duplicate of core_registry_contract_sha256, never the
        # planner's full config-tree snapshot.
        "registry_snapshot_sha256": core_registry_contract_sha256,
        "source_campaign_bindings": [dict(item) for item in source_campaign_bindings],
        "scientific_receipt_bindings": [
            dict(item) for item in scientific_receipt_bindings
        ],
        "variant_secondary_receipt_bindings": [
            dict(item) for item in variant_secondary_receipt_bindings
        ],
        "terminal_run_bindings": [dict(item) for item in terminal_run_bindings],
        "selection_universe": [dict(item) for item in selection_universe],
        "task_seed_sets": [dict(item) for item in task_seed_sets],
        "dataset_locks": dict(dataset_locks),
        "model_dispositions": [dict(item) for item in model_dispositions],
        "task_dispositions": [dict(item) for item in task_dispositions],
        "runs": [dict(item) for item in runs],
        "scientific_runs": [dict(item) for item in scientific_runs],
        "variant_secondary_runs": [
            dict(item) for item in variant_secondary_runs
        ],
        "model_candidates": [dict(item) for item in model_candidates],
        "scientific_status": "primary_succeeded_auxiliary_non_scoring",
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
    }


def campaign_universe_sha256(
    candidate_campaign_dirs: Iterable[str | Path],
) -> str:
    """Hash the exact reviewed finalist campaign-plan universe."""

    paths = tuple(
        _safe_resolve(path, "selection source campaign")
        for path in candidate_campaign_dirs
    )
    if not paths or len(paths) != len(set(paths)):
        raise TournamentError(
            "campaign universe requires distinct frozen candidate directories"
        )
    records: list[dict[str, str]] = []
    observed_plans: set[str] = set()
    for path in paths:
        _, _, binding = _candidate_plan_source(path)
        plan_sha256 = str(binding["plan_sha256"])
        if plan_sha256 in observed_plans:
            raise TournamentError("campaign universe repeats a plan")
        observed_plans.add(plan_sha256)
        records.append(
            {
                "plan_sha256": plan_sha256,
                "manifest_sha256": str(binding["manifest_sha256"]),
            }
        )
    records.sort(key=lambda item: item["plan_sha256"])
    return canonical_sha256(records)


def _derive_campaign_universe(
    candidate_campaign_dirs: Iterable[str | Path],
    *,
    universe_kind: str,
) -> dict[str, Any]:
    """Derive one pre-scoring campaign universe of the requested kind.

    The screening kind allows only the prospectively frozen three-seed
    ``frozen_screen`` wave; the finalist kind allows only the full specialist,
    adaptation, and authorized conditional waves at exactly five seeds.
    """

    policy = _CAMPAIGN_UNIVERSE_POLICIES.get(universe_kind)
    if policy is None:
        raise TournamentError(
            f"unsupported campaign universe kind {universe_kind!r}"
        )
    allowed_waves, required_seeds, schema_version = policy
    paths = tuple(
        _safe_resolve(path, f"{universe_kind} source campaign")
        for path in candidate_campaign_dirs
    )
    if not paths or len(paths) != len(set(paths)):
        raise TournamentError(
            f"{universe_kind} universe requires distinct frozen campaign directories"
        )
    campaigns: dict[str, dict[str, Any]] = {}
    plans: dict[str, dict[str, Any]] = {}
    for path in paths:
        _, plan, binding = _candidate_plan_source(path)
        plan_sha256 = str(binding["plan_sha256"])
        if plan_sha256 in campaigns:
            raise TournamentError(
                f"{universe_kind} universe repeats a campaign plan"
            )
        wave = str(plan.get("campaign", {}).get("wave", ""))
        if wave not in allowed_waves:
            raise TournamentError(
                f"{universe_kind} universe rejects campaign wave {wave!r}; "
                f"allowed waves are {sorted(allowed_waves)}"
            )
        campaigns[plan_sha256] = binding
        plans[plan_sha256] = plan
    core_registry_hashes = {
        _sha256_identifier(
            plan.get("core_registry_contract_sha256"),
            "core_registry_contract_sha256",
        )
        for plan in plans.values()
    }
    if len(core_registry_hashes) != 1:
        raise TournamentError(
            f"{universe_kind} campaigns do not share one core scientific "
            "registry contract"
        )

    run_records: dict[str, dict[str, Any]] = {}
    configuration_records: dict[
        tuple[str, str, str, str], list[dict[str, Any]]
    ] = {}
    for plan in plans.values():
        raw_runs = plan.get("runs")
        if not isinstance(raw_runs, list):
            raise TournamentError(
                f"{universe_kind} source campaign lacks a run array"
            )
        for raw_run in raw_runs:
            if not isinstance(raw_run, Mapping):
                raise TournamentError(
                    f"{universe_kind} source run must be an object"
                )
            metadata = raw_run.get("metadata")
            if (
                isinstance(metadata, Mapping)
                and metadata.get("primary_endpoint_scoring_allowed") is False
            ):
                continue
            if raw_run.get("stage") in {
                "prediction_first_stress",
                "sealed_inference",
            }:
                raise TournamentError(
                    "a campaign universe cannot pre-register stress or sealed runs"
                )
            run_id = _sha256_identifier(
                raw_run.get("run_id"), "campaign universe run_id"
            )
            if run_id in run_records:
                raise TournamentError(
                    f"{universe_kind} universe repeats a primary run"
                )
            configuration = _candidate_configuration_identity(raw_run)
            configuration_sha256 = canonical_sha256(configuration)
            record = {
                "run_id": run_id,
                "task_id": raw_run["task_id"],
                "model_id": raw_run["model_id"],
                "adaptation_regime": raw_run["adaptation_regime"],
                "candidate_configuration_sha256": configuration_sha256,
                "seed": raw_run["seed"],
                "fold": raw_run["fold"],
                "stage": raw_run["stage"],
            }
            run_records[run_id] = record
            key = (
                str(record["task_id"]),
                str(record["model_id"]),
                str(record["adaptation_regime"]),
                configuration_sha256,
            )
            configuration_records.setdefault(key, []).append(record)
    if not run_records:
        raise TournamentError(
            f"{universe_kind} universe contains no primary-scoring runs"
        )

    task_seed_authorities: dict[str, tuple[int, ...]] = {}
    for (task_id, model_id, regime, configuration_sha256), records in sorted(
        configuration_records.items()
    ):
        seeds = tuple(sorted(int(record["seed"]) for record in records))
        if len(seeds) != required_seeds or len(set(seeds)) != required_seeds:
            raise TournamentError(
                f"every {universe_kind} configuration must schedule exactly "
                f"{required_seeds} distinct seeds: task={task_id}, "
                f"model={model_id}, regime={regime}, "
                f"configuration={configuration_sha256}"
            )
        prior = task_seed_authorities.get(task_id)
        if prior is not None and prior != seeds:
            raise TournamentError(
                f"all {universe_kind} configurations for {task_id} must share "
                f"one exact {required_seeds}-seed authority"
            )
        task_seed_authorities[task_id] = seeds
    source_bindings = [campaigns[key] for key in sorted(campaigns)]
    identity = {
        "schema_version": schema_version,
        "universe_kind": universe_kind,
        "expected_seed_count": required_seeds,
        "campaign_universe_sha256": canonical_sha256(
            [
                {
                    "plan_sha256": item["plan_sha256"],
                    "manifest_sha256": item["manifest_sha256"],
                }
                for item in source_bindings
            ]
        ),
        "core_registry_contract_sha256": next(iter(core_registry_hashes)),
        "source_campaign_bindings": source_bindings,
        "selection_universe": [
            run_records[run_id] for run_id in sorted(run_records)
        ],
        "task_seed_sets": [
            {
                "task_id": task_id,
                "seeds": list(task_seed_authorities[task_id]),
            }
            for task_id in sorted(task_seed_authorities)
        ],
        "created_before_development_scoring": True,
        "development_outcomes_used": False,
        "sealed_results_used": False,
    }
    return identity


def freeze_campaign_universe(
    *,
    candidate_campaign_dirs: Iterable[str | Path],
    universe_kind: str,
    output_root: str | Path,
) -> Path:
    """Freeze one complete pre-scoring campaign universe."""

    identity = _derive_campaign_universe(
        candidate_campaign_dirs, universe_kind=universe_kind
    )
    payload = {"universe_id": canonical_sha256(identity), **identity}
    root = _safe_resolve(output_root, "campaign universe output root")
    root.mkdir(parents=True, exist_ok=True)
    target = root / f"campaign-universe--{payload['universe_id']}"
    if target.exists():
        raise TournamentError(f"campaign universe exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".campaign-universe.", dir=root))
    try:
        write_json_exclusive(
            staging / _CAMPAIGN_UNIVERSE_DOCUMENT, payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": _CAMPAIGN_UNIVERSE_ARTIFACT_CLASS,
                "universe_id": payload["universe_id"],
            },
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_campaign_universe(target, expected_universe_kind=universe_kind)
    return target


def verify_campaign_universe(
    path: str | Path, *, expected_universe_kind: str | None = None
) -> dict[str, Any]:
    """Recursively verify one pre-scoring campaign universe authority."""

    root = _safe_resolve(path, "campaign universe")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / _CAMPAIGN_UNIVERSE_DOCUMENT).read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(
            f"invalid frozen campaign universe: {exc}"
        ) from exc
    if not isinstance(payload, Mapping):
        raise TournamentError("campaign universe must be a JSON object")
    universe_kind = payload.get("universe_kind")
    policy = (
        _CAMPAIGN_UNIVERSE_POLICIES.get(universe_kind)
        if isinstance(universe_kind, str)
        else None
    )
    if policy is None:
        raise TournamentError("campaign universe declares no supported kind")
    if (
        expected_universe_kind is not None
        and universe_kind != expected_universe_kind
    ):
        raise TournamentError(
            f"campaign universe is {universe_kind!r}, not the required "
            f"{expected_universe_kind!r}"
        )
    _, required_seeds, schema_version = policy
    if payload.get("schema_version") != schema_version:
        raise TournamentError("unsupported campaign universe schema")
    if payload.get("expected_seed_count") != required_seeds:
        raise TournamentError("campaign universe seed-count authority differs")
    claimed = _sha256_identifier(payload.get("universe_id"), "universe_id")
    identity = dict(payload)
    identity.pop("universe_id", None)
    if canonical_sha256(identity) != claimed:
        raise TournamentError("campaign universe identity hash mismatch")
    bindings = payload.get("source_campaign_bindings")
    if not isinstance(bindings, list) or not bindings:
        raise TournamentError("campaign universe lacks source campaigns")
    expected = _derive_campaign_universe(
        [
            binding.get("path")
            for binding in bindings
            if isinstance(binding, Mapping)
        ],
        universe_kind=universe_kind,
    )
    if identity != expected:
        raise TournamentError("campaign universe does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": _CAMPAIGN_UNIVERSE_ARTIFACT_CLASS,
        "universe_id": claimed,
    }:
        raise TournamentError("campaign universe manifest metadata mismatch")
    return dict(payload)


def freeze_finalist_campaign_universe(
    *,
    candidate_campaign_dirs: Iterable[str | Path],
    output_root: str | Path,
) -> Path:
    """Freeze the complete five-seed finalist campaign universe before scoring."""

    return freeze_campaign_universe(
        candidate_campaign_dirs=candidate_campaign_dirs,
        universe_kind=UNIVERSE_KIND_FINALIST,
        output_root=output_root,
    )


def verify_finalist_campaign_universe(path: str | Path) -> dict[str, Any]:
    """Verify one pre-scoring finalist campaign authority."""

    return verify_campaign_universe(
        path, expected_universe_kind=UNIVERSE_KIND_FINALIST
    )


def freeze_screening_campaign_universe(
    *,
    candidate_campaign_dirs: Iterable[str | Path],
    output_root: str | Path,
) -> Path:
    """Freeze the three-seed frozen-screen universe before screening scores."""

    return freeze_campaign_universe(
        candidate_campaign_dirs=candidate_campaign_dirs,
        universe_kind=UNIVERSE_KIND_SCREENING,
        output_root=output_root,
    )


def verify_screening_campaign_universe(path: str | Path) -> dict[str, Any]:
    """Verify one pre-scoring screening campaign authority."""

    return verify_campaign_universe(
        path, expected_universe_kind=UNIVERSE_KIND_SCREENING
    )


def _derive_selection_ledger(
    *,
    ledger_authority: str,
    campaign_universe_dir: str | Path,
    reviewed_universe_id: str,
    receipt_paths: Iterable[str | Path],
    variant_secondary_receipt_paths: Iterable[str | Path] = (),
    terminal_execution_attempt_paths: Iterable[str | Path] = (),
) -> dict[str, Any]:
    policy = _LEDGER_AUTHORITY_POLICIES.get(ledger_authority)
    if policy is None:
        raise TournamentError(
            f"unsupported selection ledger authority {ledger_authority!r}"
        )
    expected_kind, expected_seed_count = policy
    universe_root = _safe_resolve(campaign_universe_dir, "campaign universe")
    universe = verify_campaign_universe(
        universe_root, expected_universe_kind=expected_kind
    )
    universe_binding, _ = _frozen_document_binding(
        universe_root,
        filename=_CAMPAIGN_UNIVERSE_DOCUMENT,
        artifact_class=_CAMPAIGN_UNIVERSE_ARTIFACT_CLASS,
    )
    universe_binding = {
        **universe_binding,
        "universe_id": universe["universe_id"],
    }
    if (
        _sha256_identifier(reviewed_universe_id, "reviewed_universe_id")
        != universe["universe_id"]
    ):
        raise TournamentError(
            "campaign universe was not independently reviewed under this ID"
        )
    if (
        universe.get("created_before_development_scoring") is not True
        or universe.get("development_outcomes_used") is not False
        or universe.get("sealed_results_used") is not False
    ):
        raise TournamentError(
            "campaign universe does not declare a pre-scoring provenance"
        )
    explicit_campaign_paths = tuple(
        _safe_resolve(binding["path"], "selection source campaign")
        for binding in universe["source_campaign_bindings"]
    )
    if not explicit_campaign_paths or len(explicit_campaign_paths) != len(
        set(explicit_campaign_paths)
    ):
        raise TournamentError(
            "selection ledger requires a distinct explicit campaign universe"
        )
    campaigns: dict[str, dict[str, Any]] = {}
    plans: dict[str, dict[str, Any]] = {}
    for path in explicit_campaign_paths:
        _, plan, binding = _candidate_plan_source(path)
        plan_sha256 = str(binding["plan_sha256"])
        if plan_sha256 in campaigns:
            raise TournamentError("selection campaign universe repeats a plan")
        campaigns[plan_sha256] = binding
        plans[plan_sha256] = plan
    if [campaigns[key] for key in sorted(campaigns)] != [
        dict(item) for item in universe["source_campaign_bindings"]
    ]:
        raise TournamentError("campaign universe source bindings drifted")

    paths = tuple(
        _safe_resolve(path, "ScientificRunReceipt") for path in receipt_paths
    )
    if not paths or len(paths) != len(set(paths)):
        raise TournamentError(
            "selection candidate ledger requires distinct ScientificRunReceipt directories"
        )
    receipt_pairs = tuple(_receipt_binding(path) for path in paths)
    receipt_pairs = tuple(sorted(receipt_pairs, key=lambda item: item[0]["scientific_receipt_id"]))
    if len({item[0]["scientific_receipt_id"] for item in receipt_pairs}) != len(receipt_pairs):
        raise TournamentError("selection candidate ledger repeats a scientific receipt")
    secondary_paths = tuple(
        _safe_resolve(path, "variant secondary run receipt")
        for path in variant_secondary_receipt_paths
    )
    if len(secondary_paths) != len(set(secondary_paths)):
        raise TournamentError(
            "selection candidate ledger requires distinct secondary receipt directories"
        )
    secondary_pairs = tuple(
        sorted(
            (_variant_secondary_receipt_binding(path) for path in secondary_paths),
            key=lambda item: item[0]["secondary_run_receipt_id"],
        )
    )
    if len(
        {item[0]["secondary_run_receipt_id"] for item in secondary_pairs}
    ) != len(secondary_pairs):
        raise TournamentError("selection candidate ledger repeats a secondary receipt")

    for _, receipt in receipt_pairs:
        for field in ("candidate_campaign_binding", "baseline_campaign_binding"):
            plan, binding = _campaign_from_receipt_binding(receipt[field])
            plan_sha256 = binding["plan_sha256"]
            if campaigns.get(plan_sha256) != binding or plans.get(plan_sha256) != plan:
                raise TournamentError(
                    "scientific receipt campaign is outside the explicit selection universe"
                )
    for _, receipt in secondary_pairs:
        plan, binding = _campaign_from_receipt_binding(
            receipt["candidate_campaign_binding"]
        )
        plan_sha256 = binding["plan_sha256"]
        if campaigns.get(plan_sha256) != binding or plans.get(plan_sha256) != plan:
            raise TournamentError(
                "secondary receipt campaign is outside the explicit selection universe"
            )
    terminal_paths = tuple(
        _safe_resolve(path, "terminal run execution attempt")
        for path in terminal_execution_attempt_paths
    )
    if len(terminal_paths) != len(set(terminal_paths)):
        raise TournamentError("selection ledger repeats a terminal run attempt path")
    terminal_pairs = tuple(
        sorted(
            (_terminal_execution_attempt_binding(path) for path in terminal_paths),
            key=lambda item: item[0]["run_id"],
        )
    )
    terminal_run_ids: set[str] = set()
    for binding, payload in terminal_pairs:
        run_id = str(binding["run_id"])
        if run_id in terminal_run_ids:
            raise TournamentError("selection ledger repeats a terminal run disposition")
        terminal_run_ids.add(run_id)
        plan_sha256 = str(binding["plan_sha256"])
        plan = plans.get(plan_sha256)
        campaign = campaigns.get(plan_sha256)
        if plan is None or campaign is None:
            raise TournamentError(
                "terminal run campaign is outside the explicit selection universe"
            )
        if payload.get("candidate_path") != campaign["path"]:
            raise TournamentError("terminal run changed its candidate campaign path")
        retry_policy = plan.get("retry_policy")
        max_attempts = (
            retry_policy.get("max_attempts")
            if isinstance(retry_policy, Mapping)
            else None
        )
        if (
            isinstance(max_attempts, bool)
            or not isinstance(max_attempts, int)
            or int(binding["attempt"]) < max_attempts
        ):
            raise TournamentError(
                "terminal run disposition does not prove retry exhaustion"
            )
    registry_hashes = {
        _sha256_identifier(
            plan.get("core_registry_contract_sha256"),
            "core_registry_contract_sha256",
        )
        for plan in plans.values()
    }
    if len(registry_hashes) != 1:
        raise TournamentError(
            "all campaigns in one selection ledger must share one core scientific "
            "registry contract"
        )

    dataset_locks: dict[str, dict[str, Any]] = {}
    dispositions: dict[str, dict[str, Any]] = {}
    task_dispositions: dict[str, dict[str, Any]] = {}
    for plan in plans.values():
        raw_locks = plan.get("dataset_locks")
        if not isinstance(raw_locks, Mapping):
            raise TournamentError("candidate plan dataset_locks must be an object")
        for dataset_id, raw in raw_locks.items():
            if not isinstance(raw, Mapping):
                raise TournamentError("candidate plan dataset lock must be an object")
            _merge_exact_record(
                dataset_locks,
                identifier=_nonempty_identifier(dataset_id, "dataset_id"),
                value=raw,
                label="dataset lock",
            )
        for raw in plan.get("model_dispositions", ()):
            if not isinstance(raw, Mapping):
                raise TournamentError("model disposition must be an object")
            _merge_exact_record(
                dispositions,
                identifier=_nonempty_identifier(raw.get("model_id"), "model_id"),
                value=raw,
                label="model disposition",
            )
        for raw in plan.get("task_dispositions", ()):
            if not isinstance(raw, Mapping):
                raise TournamentError("task disposition must be an object")
            _merge_exact_record(
                task_dispositions,
                identifier=_nonempty_identifier(raw.get("task_id"), "task_id"),
                value=raw,
                label="task disposition",
            )

    runs: dict[str, dict[str, Any]] = {}
    scientific_runs: dict[str, dict[str, Any]] = {}
    variant_secondary_runs: dict[str, dict[str, Any]] = {}
    run_keys: dict[tuple[str, str, str, str, int, str], str] = {}
    for binding, receipt in receipt_pairs:
        candidate_plan = plans[receipt["candidate_campaign_binding"]["plan_sha256"]]
        baseline_plan = plans[receipt["baseline_campaign_binding"]["plan_sha256"]]
        candidate_run = _run_from_plan(candidate_plan, str(receipt["run_id"]))
        baseline_run = _run_from_plan(baseline_plan, str(receipt["baseline_run_id"]))
        for label, run in (
            ("candidate", candidate_run),
            ("baseline", baseline_run),
        ):
            metadata = run.get("metadata")
            routing = metadata.get("dataset_routing") if isinstance(metadata, Mapping) else None
            if run.get("stage") == "prediction_first_stress" or (
                isinstance(routing, Mapping)
                and any(
                    isinstance(route, Mapping)
                    and route.get("task_partition") == "prediction_first_stress"
                    for route in routing.values()
                )
            ):
                raise TournamentError(
                    f"{label} run uses prediction-first stress evidence, which is "
                    "forbidden in model selection"
                )
        endpoint = receipt["endpoint_summary"]
        for field_name in ("split_id", "fold", "dataset_ids", "immutable_inputs"):
            if candidate_run.get(field_name) != baseline_run.get(field_name):
                raise TournamentError(
                    "candidate and baseline scientific runs differ in frozen "
                    f"comparison field {field_name}"
                )
        comparison_identity = _scientific_comparison_identity(
            run=candidate_run,
            baseline_run=baseline_run,
            receipt=receipt,
        )
        comparison_identity_sha256 = canonical_sha256(comparison_identity)
        candidate_configuration = _candidate_configuration_identity(candidate_run)
        baseline_configuration = _candidate_configuration_identity(baseline_run)
        candidate_selection_artifacts = _selection_artifact_hashes(
            receipt["candidate_execution_receipt_binding"]
        )
        baseline_selection_artifacts = _selection_artifact_hashes(
            receipt["baseline_execution_receipt_binding"]
        )
        candidate_science = {
            "run_id": candidate_run["run_id"],
            "task_id": candidate_run["task_id"],
            "model_id": candidate_run["model_id"],
            "adaptation_regime": candidate_run["adaptation_regime"],
            "seed": candidate_run["seed"],
            "fold": candidate_run["fold"],
            "family_id": receipt["family_id"],
            "endpoint_evaluator_id": receipt["endpoint_evaluator_id"],
            "endpoint_evaluator_sha256": receipt["endpoint_evaluator_sha256"],
            "comparison_identity": comparison_identity,
            "comparison_identity_sha256": comparison_identity_sha256,
            "selection_artifact_hashes": candidate_selection_artifacts["role_hashes"],
            "fit_state_artifact_policy": candidate_selection_artifacts[
                "fit_state_artifact_policy"
            ],
            "fit_state_bundle_sha256": candidate_selection_artifacts[
                "fit_state_bundle_sha256"
            ],
            "candidate_configuration": candidate_configuration,
            "candidate_configuration_sha256": canonical_sha256(
                candidate_configuration
            ),
            "metrics": {
                "absolute_primary_metric": receipt["metrics"][
                    "absolute_primary_metric"
                ]
            },
            "standard_errors": {
                "absolute_primary_metric": receipt["standard_errors"][
                    "absolute_primary_metric"
                ]
            },
            "compute_cost_status": receipt["compute_cost_status"],
            "source_scientific_receipt_ids": [receipt["scientific_receipt_id"]],
        }
        baseline_disposition = _model_disposition_for_plan(
            baseline_plan, str(baseline_run["model_id"])
        )
        baseline_science = {
            "run_id": baseline_run["run_id"],
            "task_id": baseline_run["task_id"],
            "model_id": baseline_run["model_id"],
            "adaptation_regime": baseline_run["adaptation_regime"],
            "seed": baseline_run["seed"],
            "fold": baseline_run["fold"],
            "family_id": baseline_disposition["family_id"],
            "endpoint_evaluator_id": receipt["endpoint_evaluator_id"],
            "endpoint_evaluator_sha256": receipt["endpoint_evaluator_sha256"],
            "comparison_identity": comparison_identity,
            "comparison_identity_sha256": comparison_identity_sha256,
            "selection_artifact_hashes": baseline_selection_artifacts["role_hashes"],
            "fit_state_artifact_policy": baseline_selection_artifacts[
                "fit_state_artifact_policy"
            ],
            "fit_state_bundle_sha256": baseline_selection_artifacts[
                "fit_state_bundle_sha256"
            ],
            "candidate_configuration": baseline_configuration,
            "candidate_configuration_sha256": canonical_sha256(
                baseline_configuration
            ),
            "metrics": {
                "absolute_primary_metric": endpoint["baseline_primary_metric"],
            },
            "standard_errors": {
                "absolute_primary_metric": endpoint[
                    "baseline_primary_standard_error"
                ],
            },
            "compute_cost_status": "unavailable_no_comparable_measured_basis",
            "source_scientific_receipt_ids": [receipt["scientific_receipt_id"]],
        }
        for run, science in (
            (candidate_run, candidate_science),
            (baseline_run, baseline_science),
        ):
            run_id = str(run["run_id"])
            _merge_exact_record(
                runs, identifier=run_id, value=run, label="scientifically succeeded run"
            )
            key = (
                str(run["task_id"]),
                str(run["model_id"]),
                str(run["adaptation_regime"]),
                canonical_sha256(_candidate_configuration_identity(run)),
                int(run["seed"]),
                str(run["fold"]),
            )
            prior_id = run_keys.get(key)
            if prior_id is not None and prior_id != run_id:
                raise TournamentError(
                    "selection ledger has multiple runs for one task/model/regime/seed/fold"
                )
            run_keys[key] = run_id
            prior_science = scientific_runs.get(run_id)
            if prior_science is None:
                scientific_runs[run_id] = science
            else:
                prior_sources = set(prior_science["source_scientific_receipt_ids"])
                new_sources = set(science["source_scientific_receipt_ids"])
                comparable_prior = {**prior_science, "source_scientific_receipt_ids": []}
                comparable_new = {**science, "source_scientific_receipt_ids": []}
                if comparable_prior != comparable_new:
                    raise TournamentError(
                        f"scientific result changed across receipts for run {run_id}"
                    )
                prior_science["source_scientific_receipt_ids"] = sorted(
                    prior_sources | new_sources
                )

    for binding, receipt in secondary_pairs:
        plan = plans[receipt["candidate_campaign_binding"]["plan_sha256"]]
        run = _run_from_plan(plan, str(receipt["run_id"]))
        run_id = str(run["run_id"])
        if run_id in scientific_runs:
            raise TournamentError(
                "variant secondary run cannot also be a primary scientific run"
            )
        if run_id in variant_secondary_runs:
            raise TournamentError(
                f"variant secondary run has multiple auxiliary receipts: {run_id}"
            )
        _merge_exact_record(
            runs,
            identifier=run_id,
            value=run,
            label="verified auxiliary run",
        )
        key = (
            str(run["task_id"]),
            str(run["model_id"]),
            str(run["adaptation_regime"]),
            canonical_sha256(_candidate_configuration_identity(run)),
            int(run["seed"]),
            str(run["fold"]),
        )
        prior_id = run_keys.get(key)
        if prior_id is not None and prior_id != run_id:
            raise TournamentError(
                "selection ledger has multiple runs for one task/model/regime/seed/fold"
            )
        run_keys[key] = run_id
        variant_secondary_runs[run_id] = {
            "run_id": run_id,
            "task_id": receipt["task_id"],
            "model_id": receipt["model_id"],
            "adaptation_regime": receipt["adaptation_regime"],
            "seed": receipt["seed"],
            "fold": receipt["fold"],
            "family_id": receipt["family_id"],
            "dataset_ids": list(receipt["dataset_ids"]),
            "dataset_registry_sha256s": dict(
                receipt["dataset_registry_sha256s"]
            ),
            "immutable_inputs": dict(receipt["immutable_inputs"]),
            "variant_capability": dict(receipt["variant_capability"]),
            "variant_capability_sha256": receipt[
                "variant_capability_sha256"
            ],
            "variant_capability_registry_sha256": receipt[
                "variant_capability_registry_sha256"
            ],
            "endpoint_id": receipt["endpoint_id"],
            "evaluator_id": receipt["evaluator_id"],
            "score_transform_id": receipt["score_transform_id"],
            "source_secondary_run_receipt_id": binding[
                "secondary_run_receipt_id"
            ],
            "outcomes_read": False,
            "metrics_computed": False,
            "auxiliary_status": "succeeded_non_scoring",
        }

    expected_primary_runs: dict[str, Mapping[str, Any]] = {}
    for plan in plans.values():
        plan_wave = _nonempty_identifier(
            plan.get("campaign", {}).get("wave"), "campaign.wave"
        )
        raw_plan_runs = plan.get("runs")
        if not isinstance(raw_plan_runs, list):
            raise TournamentError("selection source campaign lacks a run array")
        for run in raw_plan_runs:
            if not isinstance(run, Mapping):
                raise TournamentError("selection source run must be an object")
            metadata = run.get("metadata")
            if (
                isinstance(metadata, Mapping)
                and metadata.get("primary_endpoint_scoring_allowed") is False
            ):
                continue
            if run.get("stage") in {"prediction_first_stress", "sealed_inference"}:
                raise TournamentError(
                    "selection campaign universe cannot include stress or sealed runs"
                )
            # planner sets RunSpec.stage = campaign wave, so the wave is a
            # rederived fact here rather than a declaration.
            if str(run.get("stage")) != plan_wave:
                raise TournamentError(
                    "selection universe run stage does not match its campaign wave"
                )
            run_id = _sha256_identifier(run.get("run_id"), "selection universe run_id")
            if run_id in expected_primary_runs:
                raise TournamentError("selection campaign universe repeats a primary run")
            expected_primary_runs[run_id] = run

    successful_primary_ids = set(scientific_runs)
    if successful_primary_ids.intersection(terminal_run_ids):
        raise TournamentError("a selection run cannot both succeed and fail")
    unexpected_terminal = terminal_run_ids.difference(expected_primary_runs)
    if unexpected_terminal:
        raise TournamentError(
            "terminal dispositions include runs outside the primary selection universe"
        )
    missing_primary = set(expected_primary_runs).difference(
        successful_primary_ids | terminal_run_ids
    )
    if missing_primary:
        raise TournamentError(
            "selection ledger omits scheduled primary runs: "
            + ",".join(sorted(missing_primary))
        )
    terminal_by_run = {str(binding["run_id"]): binding for binding, _ in terminal_pairs}
    selection_universe: list[dict[str, Any]] = []
    for run_id, run in sorted(expected_primary_runs.items()):
        configuration = _candidate_configuration_identity(run)
        terminal = terminal_by_run.get(run_id)
        selection_universe.append(
            {
                "run_id": run_id,
                "task_id": run["task_id"],
                "model_id": run["model_id"],
                "adaptation_regime": run["adaptation_regime"],
                "candidate_configuration_sha256": canonical_sha256(configuration),
                "seed": run["seed"],
                "fold": run["fold"],
                "terminal_status": (
                    "succeeded" if terminal is None else terminal["status"]
                ),
                "terminal_receipt_document_sha256": (
                    None if terminal is None else terminal["document_sha256"]
                ),
            }
        )

    universe_by_task_configuration: dict[
        tuple[str, str, str, str], list[dict[str, Any]]
    ] = {}
    for record in selection_universe:
        key = (
            str(record["task_id"]),
            str(record["model_id"]),
            str(record["adaptation_regime"]),
            str(record["candidate_configuration_sha256"]),
        )
        universe_by_task_configuration.setdefault(key, []).append(record)
    task_seed_authorities: dict[str, tuple[int, ...]] = {}
    for (task_id, model_id, regime, configuration_sha256), records in sorted(
        universe_by_task_configuration.items()
    ):
        seeds = tuple(sorted(int(record["seed"]) for record in records))
        if (
            len(seeds) != expected_seed_count
            or len(set(seeds)) != expected_seed_count
        ):
            raise TournamentError(
                f"every {ledger_authority} configuration must schedule exactly "
                f"{expected_seed_count} distinct seeds: task={task_id}, "
                f"model={model_id}, regime={regime}, "
                f"configuration={configuration_sha256}"
            )
        prior = task_seed_authorities.get(task_id)
        if prior is not None and prior != seeds:
            raise TournamentError(
                f"all {ledger_authority} configurations for {task_id} must "
                f"share one exact {expected_seed_count}-seed authority"
            )
        task_seed_authorities[task_id] = seeds
    task_seed_sets = [
        {"task_id": task_id, "seeds": list(task_seed_authorities[task_id])}
        for task_id in sorted(task_seed_authorities)
    ]

    # The anti-cherry-pick control: the complete scheduled run set and the exact
    # per-task seed sets were fixed in an independently reviewed output file before
    # any development score existed.  A caller can therefore neither omit a
    # stronger campaign nor choose a convenient seed set after the fact.
    observed_projection = [
        {key: record[key] for key in _UNIVERSE_IDENTITY_KEYS}
        for record in selection_universe
    ]
    expected_projection = sorted(
        (
            {key: record[key] for key in _UNIVERSE_IDENTITY_KEYS}
            for record in universe["selection_universe"]
        ),
        key=lambda item: str(item["run_id"]),
    )
    if observed_projection != expected_projection:
        raise TournamentError(
            "scheduled selection universe differs from the pre-scoring "
            "campaign universe authority"
        )
    if task_seed_sets != [dict(item) for item in universe["task_seed_sets"]]:
        raise TournamentError(
            "task seed sets differ from the pre-scoring seed authority"
        )
    selection_universe_sha256 = canonical_sha256(observed_projection)

    scientific_records = [scientific_runs[key] for key in sorted(scientific_runs)]
    # A run appears in the ledger both as a candidate and, for the fixed
    # baseline, as its own candidate row.  Index BOTH roles, otherwise baseline
    # rows never receive an ensemble endpoint and the strongest-baseline
    # comparator becomes unrankable.  One baseline can be paired with several
    # candidates, so pick deterministically by receipt id; the baseline column
    # and the row alignment are identical across those receipts.
    receipt_index: dict[str, tuple[Path, Mapping[str, Any], str]] = {}
    for binding, receipt in sorted(
        receipt_pairs, key=lambda item: str(item[0]["scientific_receipt_id"])
    ):
        root = _safe_resolve(binding["path"], "ScientificRunReceipt")
        for role, key in (
            ("candidate", "run_id"),
            ("baseline", "baseline_run_id"),
        ):
            run_id = str(receipt[key])
            receipt_index.setdefault(run_id, (root, receipt, role))
    model_candidates = _model_candidates_from_scientific_runs(
        scientific_records, dispositions, receipt_index, expected_seed_count
    )
    universe_groups: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}
    for record in selection_universe:
        key = (
            str(record["task_id"]),
            str(record["model_id"]),
            str(record["adaptation_regime"]),
            str(record["candidate_configuration_sha256"]),
        )
        universe_groups.setdefault(key, []).append(record)
    for candidate in model_candidates:
        key = (
            str(candidate["task_id"]),
            str(candidate["model_id"]),
            str(candidate["adaptation_regime"]),
            str(candidate["candidate_configuration_sha256"]),
        )
        scheduled = sorted(
            universe_groups.get(key, []), key=lambda item: str(item["run_id"])
        )
        scheduled_ids = [str(item["run_id"]) for item in scheduled]
        succeeded_ids = list(candidate["run_ids"])
        if not set(succeeded_ids).issubset(scheduled_ids):
            raise TournamentError(
                "scientific candidate runs fall outside the campaign universe"
            )
        scheduled_seed_universe_complete = (
            bool(scheduled)
            and len({int(item["seed"]) for item in scheduled}) == len(scheduled)
            and all(item["terminal_status"] == "succeeded" for item in scheduled)
            and succeeded_ids == scheduled_ids
        )
        five_seed_complete = (
            scheduled_seed_universe_complete
            and len(scheduled) == expected_seed_count
        )
        candidate["scheduled_run_ids"] = scheduled_ids
        candidate["terminal_dispositions"] = [
            {
                "run_id": item["run_id"],
                "status": item["terminal_status"],
                "receipt_document_sha256": item[
                    "terminal_receipt_document_sha256"
                ],
            }
            for item in scheduled
            if item["terminal_status"] != "succeeded"
        ]
        # Best-model standing is inexpressible on a screening-authority ledger.
        candidate["five_seed_universe_complete"] = bool(
            five_seed_complete and ledger_authority == LEDGER_AUTHORITY_FINALIST
        )
        candidate["scheduled_seed_universe_complete"] = (
            scheduled_seed_universe_complete
        )
        candidate["eligible"] = bool(
            candidate["eligible"] and scheduled_seed_universe_complete
        )
    return _selection_ledger_identity(
        ledger_authority=ledger_authority,
        campaign_universe_binding=universe_binding,
        reviewed_universe_id=universe["universe_id"],
        campaign_universe_kind=universe["universe_kind"],
        campaign_universe_sha256_value=universe["campaign_universe_sha256"],
        campaign_universe_wave_ids=sorted(
            {str(plan.get("campaign", {}).get("wave", "")) for plan in plans.values()}
        ),
        expected_seed_count=expected_seed_count,
        selection_universe_sha256=selection_universe_sha256,
        core_registry_contract_sha256=next(iter(registry_hashes)),
        source_campaign_bindings=[campaigns[key] for key in sorted(campaigns)],
        scientific_receipt_bindings=[item[0] for item in receipt_pairs],
        variant_secondary_receipt_bindings=[item[0] for item in secondary_pairs],
        terminal_run_bindings=[item[0] for item in terminal_pairs],
        selection_universe=selection_universe,
        task_seed_sets=task_seed_sets,
        dataset_locks={key: dataset_locks[key] for key in sorted(dataset_locks)},
        model_dispositions=[dispositions[key] for key in sorted(dispositions)],
        task_dispositions=[task_dispositions[key] for key in sorted(task_dispositions)],
        runs=[runs[key] for key in sorted(runs)],
        scientific_runs=scientific_records,
        variant_secondary_runs=[
            variant_secondary_runs[key] for key in sorted(variant_secondary_runs)
        ],
        model_candidates=model_candidates,
    )


def _freeze_selection_candidate_ledger(
    *,
    ledger_authority: str,
    campaign_universe_dir: str | Path,
    reviewed_universe_id: str,
    scientific_run_receipt_dirs: Iterable[str | Path],
    variant_secondary_run_receipt_dirs: Iterable[str | Path] = (),
    terminal_execution_attempt_dirs: Iterable[str | Path] = (),
    output_root: str | Path,
) -> Path:
    identity = _derive_selection_ledger(
        ledger_authority=ledger_authority,
        campaign_universe_dir=campaign_universe_dir,
        reviewed_universe_id=reviewed_universe_id,
        receipt_paths=scientific_run_receipt_dirs,
        variant_secondary_receipt_paths=variant_secondary_run_receipt_dirs,
        terminal_execution_attempt_paths=terminal_execution_attempt_dirs,
    )
    payload = {"ledger_id": canonical_sha256(identity), **identity}
    output = _safe_resolve(output_root, "selection candidate ledger output root")
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"selection-candidates--{payload['ledger_id']}"
    if target.exists():
        raise TournamentError(f"selection candidate ledger already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".selection-candidates.", dir=output))
    try:
        write_json_exclusive(staging / "selection_candidate_ledger.json", payload, mode=0o440)
        freeze_tree(
            staging,
            {"artifact_class": "selection_candidate_ledger", "ledger_id": payload["ledger_id"]},
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_selection_candidate_ledger(target)
    return target


def freeze_finalist_candidate_ledger(
    *,
    finalist_campaign_universe_dir: str | Path,
    reviewed_universe_id: str,
    scientific_run_receipt_dirs: Iterable[str | Path],
    variant_secondary_run_receipt_dirs: Iterable[str | Path] = (),
    terminal_execution_attempt_dirs: Iterable[str | Path] = (),
    output_root: str | Path,
) -> Path:
    """Freeze the only ledger class that may support a best-model SelectionLock.

    ``reviewed_universe_id`` is the human review handshake.  Transcribe it from
    the independently reviewed universe output file; it is deliberately not
    defaulted from the directory being passed.
    """

    return _freeze_selection_candidate_ledger(
        ledger_authority=LEDGER_AUTHORITY_FINALIST,
        campaign_universe_dir=finalist_campaign_universe_dir,
        reviewed_universe_id=reviewed_universe_id,
        scientific_run_receipt_dirs=scientific_run_receipt_dirs,
        variant_secondary_run_receipt_dirs=variant_secondary_run_receipt_dirs,
        terminal_execution_attempt_dirs=terminal_execution_attempt_dirs,
        output_root=output_root,
    )


def freeze_screening_candidate_ledger(
    *,
    screening_campaign_universe_dir: str | Path,
    reviewed_universe_id: str,
    scientific_run_receipt_dirs: Iterable[str | Path],
    variant_secondary_run_receipt_dirs: Iterable[str | Path] = (),
    terminal_execution_attempt_dirs: Iterable[str | Path] = (),
    output_root: str | Path,
) -> Path:
    """Freeze a frozen-screen ledger that may shortlist but never settle the selection."""

    return _freeze_selection_candidate_ledger(
        ledger_authority=LEDGER_AUTHORITY_SCREENING,
        campaign_universe_dir=screening_campaign_universe_dir,
        reviewed_universe_id=reviewed_universe_id,
        scientific_run_receipt_dirs=scientific_run_receipt_dirs,
        variant_secondary_run_receipt_dirs=variant_secondary_run_receipt_dirs,
        terminal_execution_attempt_dirs=terminal_execution_attempt_dirs,
        output_root=output_root,
    )


_LEDGER_VERIFICATION_CACHE: OrderedDict[tuple[str, str], dict[str, Any]] = (
    OrderedDict()
)
_LEDGER_VERIFICATION_CACHE_MAX = 8


def _clear_ledger_verification_cache_for_testing() -> None:
    _LEDGER_VERIFICATION_CACHE.clear()


def verify_selection_candidate_ledger(path: str | Path) -> dict[str, Any]:
    """Recursively rederive a multi-campaign scientific candidate ledger.

    Rederivation now recomputes one endpoint bootstrap per candidate, and the
    downstream chain verifies the same read-only ledger many times over (record
    construction, record validation, record verification, and each finalist metric
    bundle).  A frozen tree cannot change, so the result is memoized on the
    resolved path plus its ARTIFACTS.json digest: a different or tampered tree
    yields a different key and is fully rederived.  Nothing about WHAT is
    verified changes; only repeated identical work is skipped.
    """

    root = _safe_resolve(path, "selection candidate ledger")
    # ARTIFACTS.json is a manifest LISTING member digests, so editing
    # selection_candidate_ledger.json in place leaves this file's own bytes
    # unchanged.  Keying on it alone would let a tampered tree be served from
    # cache without ever running verify_frozen_tree.  Verify the tree first --
    # it is a cheap walk, and the expensive part being skipped is the
    # per-candidate endpoint bootstrap, not this check -- and fold the document
    # digest into the key as well.
    try:
        verify_frozen_tree(root)
        cache_key = (
            root.as_posix(),
            sha256_file(root / "ARTIFACTS.json"),
            sha256_file(root / "selection_candidate_ledger.json"),
        )
    except (ArtifactError, OSError) as exc:
        if isinstance(exc, ArtifactError):
            raise TournamentError(
                f"invalid frozen selection candidate ledger: {exc}"
            ) from exc
        cache_key = None
    if cache_key is not None:
        cached = _LEDGER_VERIFICATION_CACHE.get(cache_key)
        if cached is not None:
            _LEDGER_VERIFICATION_CACHE.move_to_end(cache_key)
            return json.loads(json.dumps(cached))
    payload = _verify_selection_candidate_ledger_uncached(root)
    if cache_key is not None:
        _LEDGER_VERIFICATION_CACHE[cache_key] = json.loads(json.dumps(payload))
        _LEDGER_VERIFICATION_CACHE.move_to_end(cache_key)
        while len(_LEDGER_VERIFICATION_CACHE) > _LEDGER_VERIFICATION_CACHE_MAX:
            _LEDGER_VERIFICATION_CACHE.popitem(last=False)
    return payload


def _verify_selection_candidate_ledger_uncached(
    path: str | Path,
) -> dict[str, Any]:
    root = _safe_resolve(path, "selection candidate ledger")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "selection_candidate_ledger.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid frozen selection candidate ledger: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise TournamentError("selection candidate ledger must be a JSON object")
    if payload.get("schema_version") != SELECTION_CANDIDATE_LEDGER_SCHEMA_VERSION:
        raise TournamentError("unsupported selection candidate ledger schema")
    campaign_bindings = payload.get("source_campaign_bindings")
    if not isinstance(campaign_bindings, list) or not campaign_bindings:
        raise TournamentError(
            "selection candidate ledger lacks its explicit campaign universe"
        )
    # The campaign set is no longer read from the ledger's own list; it is
    # recovered from the recursively verified pre-scoring universe output file.
    universe_binding = payload.get("campaign_universe_binding")
    if not isinstance(universe_binding, Mapping):
        raise TournamentError(
            "selection candidate ledger lacks its campaign universe binding"
        )
    universe_source = _verify_bound_frozen_document(
        {
            key: universe_binding.get(key)
            for key in ("path", "manifest_sha256", "document_sha256")
        },
        filename=_CAMPAIGN_UNIVERSE_DOCUMENT,
        artifact_class=_CAMPAIGN_UNIVERSE_ARTIFACT_CLASS,
    )
    universe = verify_campaign_universe(
        universe_source,
        expected_universe_kind=payload.get("campaign_universe_kind"),
    )
    if (
        universe_binding.get("universe_id") != universe["universe_id"]
        or payload.get("reviewed_universe_id") != universe["universe_id"]
        or payload.get("campaign_universe_sha256")
        != universe["campaign_universe_sha256"]
    ):
        raise TournamentError(
            "selection ledger campaign universe binding changed"
        )
    bindings = payload.get("scientific_receipt_bindings")
    if not isinstance(bindings, list) or not bindings:
        raise TournamentError("selection candidate ledger lacks scientific receipts")
    receipt_paths = []
    for index, binding in enumerate(bindings):
        if not isinstance(binding, Mapping):
            raise TournamentError(f"scientific receipt binding {index} must be an object")
        source = _verify_bound_frozen_document(
            {
                key: binding.get(key)
                for key in ("path", "manifest_sha256", "document_sha256")
            },
            filename="scientific_run_receipt.json",
            artifact_class="scientific_run_receipt",
        )
        receipt = verify_scientific_run_receipt(source)
        if (
            binding.get("scientific_receipt_id") != receipt["scientific_receipt_id"]
            or binding.get("run_id") != receipt["run_id"]
            or binding.get("baseline_run_id") != receipt["baseline_run_id"]
        ):
            raise TournamentError("scientific receipt ledger binding changed")
        receipt_paths.append(source)
    secondary_bindings = payload.get("variant_secondary_receipt_bindings")
    if not isinstance(secondary_bindings, list):
        raise TournamentError(
            "selection candidate ledger secondary receipt bindings must be an array"
        )
    secondary_paths: list[Path] = []
    for index, binding in enumerate(secondary_bindings):
        if not isinstance(binding, Mapping):
            raise TournamentError(
                f"variant secondary receipt binding {index} must be an object"
            )
        source = _verify_bound_frozen_document(
            {
                key: binding.get(key)
                for key in ("path", "manifest_sha256", "document_sha256")
            },
            filename="variant_secondary_run_receipt.json",
            artifact_class="variant_secondary_run_receipt",
        )
        receipt = verify_variant_secondary_run_receipt(source)
        if (
            binding.get("secondary_run_receipt_id")
            != receipt["secondary_run_receipt_id"]
            or binding.get("run_id") != receipt["run_id"]
            or binding.get("model_id") != receipt["model_id"]
        ):
            raise TournamentError(
                "variant secondary receipt ledger binding changed"
            )
        secondary_paths.append(source)
    terminal_bindings = payload.get("terminal_run_bindings")
    if not isinstance(terminal_bindings, list):
        raise TournamentError(
            "selection candidate ledger terminal bindings must be an array"
        )
    terminal_paths: list[Path] = []
    for index, binding in enumerate(terminal_bindings):
        if not isinstance(binding, Mapping) or not isinstance(
            binding.get("path"), str
        ):
            raise TournamentError(f"terminal run binding {index} is invalid")
        terminal_paths.append(
            _safe_resolve(binding["path"], "terminal run execution attempt")
        )
    identity = _derive_selection_ledger(
        ledger_authority=payload.get("ledger_authority"),
        campaign_universe_dir=universe_source,
        reviewed_universe_id=payload.get("reviewed_universe_id"),
        receipt_paths=receipt_paths,
        variant_secondary_receipt_paths=secondary_paths,
        terminal_execution_attempt_paths=terminal_paths,
    )
    claimed = _sha256_identifier(payload.get("ledger_id"), "ledger_id")
    if dict(payload) != {"ledger_id": claimed, **identity}:
        raise TournamentError("selection candidate ledger does not rederive")
    if canonical_sha256(identity) != claimed:
        raise TournamentError("selection candidate ledger identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "selection_candidate_ledger",
        "ledger_id": claimed,
    }:
        raise TournamentError("selection candidate ledger manifest metadata mismatch")
    return dict(payload)


def _selection_lock_binding(
    path: str | Path,
) -> tuple[SelectionLock, dict[str, Any], Path]:
    root = _safe_resolve(path, "SelectionLock")
    try:
        lock = verify_selection_lock(root)
    except (SelectionError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid frozen SelectionLock: {exc}") from exc
    binding, _ = _frozen_document_binding(
        root, filename="selection_lock.json", artifact_class="selection_lock"
    )
    return lock, {**binding, "selection_lock_id": lock.lock_id}, root


def _ledger_binding(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, Any], Path]:
    root = _safe_resolve(path, "selection candidate ledger")
    ledger = verify_selection_candidate_ledger(root)
    binding, _ = _frozen_document_binding(
        root,
        filename="selection_candidate_ledger.json",
        artifact_class="selection_candidate_ledger",
    )
    return ledger, {**binding, "ledger_id": ledger["ledger_id"]}, root


def _verify_lock_ledger_binding(
    lock: SelectionLock,
    ledger: Mapping[str, Any],
    ledger_root: Path,
) -> None:
    if (
        ledger.get("ledger_authority") != LEDGER_AUTHORITY_FINALIST
        or ledger.get("lock_construction_allowed") is not True
        or ledger.get("champion_selection_allowed") is not True
        or ledger.get("campaign_universe_kind") != UNIVERSE_KIND_FINALIST
        or ledger.get("expected_seed_count") != _FINALIST_SEED_COUNT
    ):
        raise TournamentError(
            "a screening-authority ledger cannot support a SelectionLock"
        )
    candidate = lock.metadata.get("candidate_binding")
    expected = {
        "path": ledger_root.as_posix(),
        "manifest_sha256": sha256_file(ledger_root / "ARTIFACTS.json"),
        "artifact_class": "selection_candidate_ledger",
        "source_identity_sha256": ledger["ledger_id"],
    }
    if candidate != expected:
        raise TournamentError("SelectionLock does not recursively bind this candidate ledger")
    if (
        lock.campaign_id != ledger["ledger_id"]
        or lock.plan_sha256 != ledger["ledger_id"]
        or lock.metrics_sha256 != ledger["ledger_id"]
        or lock.candidate_manifest_sha256 != expected["manifest_sha256"]
        or lock.registry_sha256 != ledger["core_registry_contract_sha256"]
    ):
        raise TournamentError("SelectionLock ledger identities do not match")
    expected_universe = {
        "universe_id": ledger["reviewed_universe_id"],
        "universe_kind": ledger["campaign_universe_kind"],
        "campaign_universe_sha256": ledger["campaign_universe_sha256"],
        "manifest_sha256": ledger["campaign_universe_binding"]["manifest_sha256"],
        "document_sha256": ledger["campaign_universe_binding"]["document_sha256"],
        "selection_universe_sha256": ledger["selection_universe_sha256"],
    }
    if lock.metadata.get("finalist_campaign_universe_binding") != expected_universe:
        raise TournamentError(
            "SelectionLock does not bind this ledger's finalist campaign universe"
        )


def _recomputed_scientific_endpoint(
    receipt_root: Path,
    receipt: Mapping[str, Any],
    *,
    n_resamples: int,
    seed: int,
) -> EndpointBootstrapResult:
    task_id = str(receipt["task_id"])
    joined_ref = ArtifactRef.from_dict(receipt["joined_endpoint_table"])
    joined_path = joined_ref.validate(receipt_root, require_relative=True)
    try:
        return recompute_endpoint_bootstrap(
            evaluator_id=str(receipt["endpoint_evaluator_id"]),
            table_path=joined_path,
            n_resamples=n_resamples,
            seed=seed,
            parameters=receipt["endpoint_parameters"],
        )
    except EndpointPowerError as exc:
        raise TournamentError(
            f"finalist scientific endpoint does not rederive for {task_id}: {exc}"
        ) from exc


def _finalist_receipts(
    *,
    lock: SelectionLock,
    ledger: Mapping[str, Any],
    task_id: str,
) -> list[tuple[dict[str, Any], dict[str, Any], Path, EndpointBootstrapResult]]:
    decision = _task_decision(lock, task_id)
    selected_by_seed = {
        int(record["seed"]): str(record["run_id"])
        for record in decision["selected_runs"]
    }
    baseline_by_seed = {
        int(record["seed"]): str(record["run_id"])
        for record in decision["baseline_runs"]
    }
    if set(selected_by_seed) != set(baseline_by_seed) or len(selected_by_seed) != 5:
        raise TournamentError("finalist metric bundle requires five aligned locked seeds")
    receipt_bindings = ledger["scientific_receipt_bindings"]
    matches: list[tuple[dict[str, Any], dict[str, Any], Path, EndpointBootstrapResult]] = []
    used_receipts: set[str] = set()
    for seed in sorted(selected_by_seed):
        candidates = [
            binding
            for binding in receipt_bindings
            if binding.get("run_id") == selected_by_seed[seed]
            and binding.get("baseline_run_id") == baseline_by_seed[seed]
        ]
        if len(candidates) != 1:
            raise TournamentError(
                f"locked finalist seed {seed} requires exactly one aligned ScientificRunReceipt"
            )
        binding = dict(candidates[0])
        root = _safe_resolve(binding["path"], "ScientificRunReceipt")
        receipt = verify_scientific_run_receipt(root)
        if (
            receipt["task_id"] != task_id
            or receipt["seed"] != seed
            or receipt["run_id"] != selected_by_seed[seed]
            or receipt["baseline_run_id"] != baseline_by_seed[seed]
        ):
            raise TournamentError("finalist ScientificRunReceipt changed its locked pairing")
        receipt_id = str(receipt["scientific_receipt_id"])
        if receipt_id in used_receipts:
            raise TournamentError("finalist metric bundle repeats a ScientificRunReceipt")
        used_receipts.add(receipt_id)
        matches.append(
            (
                binding,
                receipt,
                root,
                _recomputed_scientific_endpoint(
                    root,
                    receipt,
                    n_resamples=_POWER_BOOTSTRAP_RESAMPLES,
                    seed=_POWER_BOOTSTRAP_SEED,
                ),
            )
        )
    return matches


def _seed_stability_record(
    endpoints: Sequence[EndpointBootstrapResult], seeds: Sequence[int]
) -> dict[str, Any]:
    """Per-seed diagnostics.  These never order candidates."""

    per_seed = {
        str(int(seed)): float(endpoint.observed_effect)
        for seed, endpoint in zip(seeds, endpoints, strict=True)
    }
    per_seed_metric = {
        str(int(seed)): float(endpoint.candidate_primary_metric)
        for seed, endpoint in zip(seeds, endpoints, strict=True)
    }
    values = list(per_seed_metric.values())
    mean = fsum(values) / len(values)
    variance = (
        fsum((value - mean) ** 2 for value in values) / (len(values) - 1)
        if len(values) > 1
        else 0.0
    )
    positive = sum(1 for value in per_seed.values() if value > 0.0)
    return {
        "per_seed_observed_effect": dict(sorted(per_seed.items())),
        "per_seed_absolute_primary_metric": dict(sorted(per_seed_metric.items())),
        "per_seed_mean_absolute_primary_metric": mean,
        "seed_stability_sd": sqrt(variance),
        "positive_seed_count": positive,
        "evaluated_seed_count": len(values),
        "used_for_ranking": False,
        "seeds_are_biological_replicates": False,
    }


def _aggregate_finalist_endpoint(
    receipts: Sequence[
        tuple[Mapping[str, Any], Mapping[str, Any], Path, EndpointBootstrapResult]
    ],
    *,
    ensemble: EndpointBootstrapResult,
    alignment: Mapping[str, Any],
    stability: Mapping[str, Any],
) -> dict[str, Any]:
    """Report the endpoint of the five-seed ensemble, not a mean of five.

    The per-seed endpoints are still verified so their evaluator, unit set, row
    set, strata, bootstrap seed, and replicate count are provably identical, but
    the reported statistic is the single task-native endpoint evaluated once on
    the ensemble prediction.
    """

    endpoints = [item[3] for item in receipts]
    if len(endpoints) != _ENSEMBLE_SEED_COUNT:
        raise TournamentError("finalist endpoint aggregation requires exactly five seeds")
    invariant_fields = (
        "evaluator_id",
        "n_units",
        "independent_unit_counts",
        "unit_set_sha256",
        "n_rows",
        "row_set_sha256",
        "strata",
        "seed",
        "n_resamples",
    )
    first = endpoints[0]
    for endpoint in endpoints[1:]:
        if any(getattr(endpoint, field) != getattr(first, field) for field in invariant_fields):
            raise TournamentError(
                "finalist seed endpoints differ in evaluator or row/unit identity"
            )
    if any(
        getattr(ensemble, field) != getattr(first, field)
        for field in (
            "evaluator_id",
            "n_units",
            "unit_set_sha256",
            "n_rows",
            "row_set_sha256",
            "strata",
            "seed",
            "n_resamples",
        )
    ):
        raise TournamentError(
            "development ensemble endpoint differs from its seed evaluations in "
            "evaluator or row/unit identity"
        )
    return {
        "schema_version": _FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION,
        "aggregation_method": _DEVELOPMENT_ENSEMBLE_AGGREGATION,
        "development_ensemble_policy_id": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "evaluator_id": ensemble.evaluator_id,
        "observed_effect": ensemble.observed_effect,
        "candidate_primary_metric": ensemble.candidate_primary_metric,
        "baseline_primary_metric": ensemble.baseline_primary_metric,
        "candidate_primary_bootstrap": list(ensemble.candidate_primary_bootstrap),
        "baseline_primary_bootstrap": list(ensemble.baseline_primary_bootstrap),
        "paired_difference_bootstrap": list(ensemble.paired_difference_bootstrap),
        "candidate_primary_standard_error": _sample_standard_error(
            ensemble.candidate_primary_bootstrap
        ),
        "baseline_primary_standard_error": _sample_standard_error(
            ensemble.baseline_primary_bootstrap
        ),
        "paired_effect_standard_error": _sample_standard_error(
            ensemble.paired_difference_bootstrap
        ),
        "ensemble_alignment_sha256": alignment["ensemble_alignment_sha256"],
        "ensemble_rows_sha256": alignment["ensemble_rows_sha256"],
        "seed_stability": dict(stability),
        "n_units": ensemble.n_units,
        "independent_unit_counts": dict(ensemble.independent_unit_counts),
        "unit_set_sha256": ensemble.unit_set_sha256,
        "n_rows": ensemble.n_rows,
        "row_set_sha256": ensemble.row_set_sha256,
        "strata": list(ensemble.strata),
        "bootstrap_seed": ensemble.seed,
        "n_resamples": ensemble.n_resamples,
        "selected_seed_count": len(endpoints),
    }


def _finalist_metric_sources(
    *,
    selection_lock_dir: str | Path,
    selection_candidate_ledger_dir: str | Path,
    task_id: str,
) -> dict[str, Any]:
    lock, lock_binding, _ = _selection_lock_binding(selection_lock_dir)
    ledger, ledger_binding, ledger_root = _ledger_binding(
        selection_candidate_ledger_dir
    )
    _verify_lock_ledger_binding(lock, ledger, ledger_root)
    receipts = _finalist_receipts(lock=lock, ledger=ledger, task_id=task_id)
    evaluator_ids = {str(item[1]["endpoint_evaluator_id"]) for item in receipts}
    evaluator_hashes = {str(item[1]["endpoint_evaluator_sha256"]) for item in receipts}
    parameter_hashes = {
        canonical_sha256(item[1]["endpoint_parameters"]) for item in receipts
    }
    task_contract_hashes = {
        str(item[1]["task_contract_sha256"]) for item in receipts
    }
    outcome_ids = {
        str(item[1]["development_outcome_binding"]["outcome_bundle_id"])
        for item in receipts
    }
    outcome_manifests = {
        str(item[1]["development_outcome_binding"]["manifest_sha256"])
        for item in receipts
    }
    if not (
        len(evaluator_ids)
        == len(evaluator_hashes)
        == len(parameter_hashes)
        == len(task_contract_hashes)
        == len(outcome_ids)
        == len(outcome_manifests)
        == 1
    ):
        raise TournamentError(
            "finalist ScientificRunReceipts differ in evaluator, parameters, or outcomes"
        )
    decision = _task_decision(lock, task_id)
    endpoint_parameters = dict(receipts[0][1]["endpoint_parameters"])
    ensemble_rows, alignment = _development_ensemble_rows(
        task_id=task_id,
        per_seed=[
            (int(receipt["seed"]), root, receipt)
            for _binding, receipt, root, _endpoint in receipts
        ],
        endpoint_parameters=endpoint_parameters,
    )
    stability = _seed_stability_record(
        [item[3] for item in sorted(receipts, key=lambda x: int(x[1]["seed"]))],
        [int(item[1]["seed"]) for item in sorted(receipts, key=lambda x: int(x[1]["seed"]))],
    )
    if stability["positive_seed_count"] < 4:
        raise TournamentError(
            "finalist development candidate fails the four-of-five positive "
            "seed-direction stability requirement"
        )
    ensemble_staging = Path(
        tempfile.mkdtemp(prefix=".development-ensemble.")
    )
    try:
        ensemble_path = ensemble_staging / "ensemble_endpoint.tsv"
        write_text_exclusive(
            ensemble_path,
            _tsv_text(_JOINED_ENDPOINT_FIELDS[task_id], ensemble_rows),
            mode=0o440,
        )
        try:
            ensemble = recompute_endpoint_bootstrap(
                evaluator_id=next(iter(evaluator_ids)),
                table_path=ensemble_path,
                n_resamples=_POWER_BOOTSTRAP_RESAMPLES,
                seed=_POWER_BOOTSTRAP_SEED,
                parameters=endpoint_parameters,
            )
        except EndpointPowerError as exc:
            raise TournamentError(
                f"development ensemble endpoint failed to recompute: {exc}"
            ) from exc
        ensemble_text = ensemble_path.read_text(encoding="utf-8")
    finally:
        shutil.rmtree(ensemble_staging, ignore_errors=True)
    return {
        "lock": lock,
        "selection_lock_binding": lock_binding,
        "ledger_binding": ledger_binding,
        "receipts": receipts,
        "task_id": task_id,
        "selected_run_ids": sorted(_locked_run_ids(decision, "selected")),
        "baseline_run_ids": sorted(_locked_run_ids(decision, "baseline")),
        "endpoint_evaluator_id": next(iter(evaluator_ids)),
        "endpoint_evaluator_sha256": next(iter(evaluator_hashes)),
        "endpoint_parameters": endpoint_parameters,
        "task_contract_sha256": next(iter(task_contract_hashes)),
        "development_outcome_bundle_id": next(iter(outcome_ids)),
        "development_outcome_manifest_sha256": next(iter(outcome_manifests)),
        "ensemble_rows": ensemble_rows,
        "ensemble_table_text": ensemble_text,
        "ensemble_alignment": alignment,
        "endpoint_result": _aggregate_finalist_endpoint(
            receipts,
            ensemble=ensemble,
            alignment=alignment,
            stability=stability,
        ),
    }


def _finalist_metric_identity(
    *,
    sources: Mapping[str, Any],
    endpoint_result_ref: ArtifactRef,
    ensemble_table_ref: ArtifactRef,
) -> dict[str, Any]:
    receipt_bindings = [
        {
            **dict(binding),
            "seed": receipt["seed"],
            "endpoint_summary_sha256": canonical_sha256(receipt["endpoint_summary"]),
        }
        for binding, receipt, _, _ in sources["receipts"]
    ]
    receipt_bindings.sort(key=lambda item: int(item["seed"]))
    result = sources["endpoint_result"]
    return {
        "schema_version": FINALIST_DEVELOPMENT_METRIC_BUNDLE_SCHEMA_VERSION,
        "selection_lock_binding": dict(sources["selection_lock_binding"]),
        "selection_candidate_ledger_binding": dict(sources["ledger_binding"]),
        "task_id": sources["task_id"],
        "selected_run_ids": list(sources["selected_run_ids"]),
        "baseline_run_ids": list(sources["baseline_run_ids"]),
        "scientific_receipt_bindings": receipt_bindings,
        "endpoint_evaluator_id": sources["endpoint_evaluator_id"],
        "endpoint_evaluator_sha256": sources["endpoint_evaluator_sha256"],
        "endpoint_parameters": dict(sources["endpoint_parameters"]),
        "endpoint_parameters_sha256": canonical_sha256(sources["endpoint_parameters"]),
        "task_contract_sha256": sources["task_contract_sha256"],
        "development_outcome_bundle_id": sources["development_outcome_bundle_id"],
        "development_outcome_manifest_sha256": sources[
            "development_outcome_manifest_sha256"
        ],
        "endpoint_result": endpoint_result_ref.to_dict(),
        "endpoint_result_sha256": canonical_sha256(result),
        "development_ensemble_policy_id": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "ensemble_table": ensemble_table_ref.to_dict(),
        "ensemble_alignment_sha256": result["ensemble_alignment_sha256"],
        "ensemble_rows_sha256": result["ensemble_rows_sha256"],
        "seed_stability": dict(result["seed_stability"]),
        "observed_effect": result["observed_effect"],
        "candidate_primary_metric": result["candidate_primary_metric"],
        "baseline_primary_metric": result["baseline_primary_metric"],
        "n_units": result["n_units"],
        "unit_set_sha256": result["unit_set_sha256"],
        "n_rows": result["n_rows"],
        "row_set_sha256": result["row_set_sha256"],
        "power_method_id": "empirical_paired_endpoint_bootstrap_v1",
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
    }


def freeze_finalist_development_metric_bundle(
    *,
    selection_lock_dir: str | Path,
    selection_candidate_ledger_dir: str | Path,
    task_id: str,
    output_root: str | Path,
) -> Path:
    """Freeze the exact five-seed finalist development endpoint distribution."""

    sources = _finalist_metric_sources(
        selection_lock_dir=selection_lock_dir,
        selection_candidate_ledger_dir=selection_candidate_ledger_dir,
        task_id=task_id,
    )
    output = _safe_resolve(output_root, "finalist metric output root")
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".finalist-development-metrics.", dir=output))
    try:
        result_path = staging / "finalist_endpoint_distribution.json"
        write_json_exclusive(result_path, sources["endpoint_result"], mode=0o440)
        result_ref = ArtifactRef.from_path(
            result_path,
            relative_to=staging,
            media_type="application/json",
            role=f"finalist_development_endpoint_distribution:{task_id}",
        )
        ensemble_path = staging / "ensemble_endpoint.tsv"
        write_text_exclusive(
            ensemble_path, sources["ensemble_table_text"], mode=0o440
        )
        ensemble_ref = ArtifactRef.from_path(
            ensemble_path,
            relative_to=staging,
            media_type="text/tab-separated-values",
            role=f"development_ensemble_endpoint_table:{task_id}",
        )
        identity = _finalist_metric_identity(
            sources=sources,
            endpoint_result_ref=result_ref,
            ensemble_table_ref=ensemble_ref,
        )
        payload = {"metric_bundle_id": canonical_sha256(identity), **identity}
        write_json_exclusive(
            staging / "finalist_development_metric_bundle.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "finalist_development_metric_bundle",
                "metric_bundle_id": payload["metric_bundle_id"],
            },
        )
        target = output / f"finalist-development-metrics--{payload['metric_bundle_id']}"
        if target.exists():
            raise TournamentError(f"finalist development metric bundle exists: {target}")
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_finalist_development_metric_bundle(target, reverify_sources=True)
    return target


def verify_finalist_development_metric_bundle(
    path: str | Path, *, reverify_sources: bool = True
) -> dict[str, Any]:
    """Verify and optionally recursively rederive a five-seed finalist metric bundle."""

    root = _safe_resolve(path, "finalist development metric bundle")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "finalist_development_metric_bundle.json").read_text(
                encoding="utf-8"
            )
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid finalist development metric bundle: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != FINALIST_DEVELOPMENT_METRIC_BUNDLE_SCHEMA_VERSION:
        raise TournamentError("unsupported finalist development metric bundle schema")
    claimed = _sha256_identifier(payload.get("metric_bundle_id"), "metric_bundle_id")
    try:
        result_ref = ArtifactRef.from_dict(payload.get("endpoint_result"))
    except (ContractError, TypeError) as exc:
        raise TournamentError(f"invalid finalist endpoint ArtifactRef: {exc}") from exc
    if (
        result_ref.media_type != "application/json"
        or result_ref.role
        != f"finalist_development_endpoint_distribution:{payload.get('task_id')}"
    ):
        raise TournamentError("finalist endpoint distribution ArtifactRef is invalid")
    result_path = result_ref.validate(root, require_relative=True)
    try:
        result = json.loads(result_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"cannot read finalist endpoint distribution: {exc}") from exc
    if not isinstance(result, Mapping) or canonical_sha256(result) != payload.get(
        "endpoint_result_sha256"
    ):
        raise TournamentError("finalist endpoint distribution identity changed")
    result_fields = {
        "schema_version",
        "aggregation_method",
        "evaluator_id",
        "observed_effect",
        "candidate_primary_metric",
        "baseline_primary_metric",
        "candidate_primary_bootstrap",
        "baseline_primary_bootstrap",
        "paired_difference_bootstrap",
        "candidate_primary_standard_error",
        "baseline_primary_standard_error",
        "paired_effect_standard_error",
        "n_units",
        "independent_unit_counts",
        "unit_set_sha256",
        "n_rows",
        "row_set_sha256",
        "strata",
        "bootstrap_seed",
        "n_resamples",
        "selected_seed_count",
        "development_ensemble_policy_id",
        "ensemble_alignment_sha256",
        "ensemble_rows_sha256",
        "seed_stability",
    }
    if set(result) != result_fields:
        raise TournamentError("finalist endpoint distribution has the wrong exact schema")
    _sha256_identifier(
        result.get("ensemble_alignment_sha256"), "ensemble_alignment_sha256"
    )
    _sha256_identifier(
        result.get("ensemble_rows_sha256"), "ensemble_rows_sha256"
    )
    stability = result.get("seed_stability")
    if (
        not isinstance(stability, Mapping)
        or stability.get("used_for_ranking") is not False
        or stability.get("seeds_are_biological_replicates") is not False
        or stability.get("evaluated_seed_count") != _ENSEMBLE_SEED_COUNT
        or not isinstance(stability.get("positive_seed_count"), int)
        or isinstance(stability.get("positive_seed_count"), bool)
        or int(stability["positive_seed_count"]) < 4
    ):
        raise TournamentError(
            "finalist seed stability record is missing, ranking-bearing, or "
            "fails the four-of-five positive seed-direction requirement"
        )
    if (
        result.get("schema_version")
        != _FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION
        or result.get("aggregation_method")
        != _DEVELOPMENT_ENSEMBLE_AGGREGATION
        or result.get("development_ensemble_policy_id")
        != DEVELOPMENT_ENSEMBLE_POLICY_ID
        or result.get("bootstrap_seed") != _POWER_BOOTSTRAP_SEED
        or result.get("n_resamples") != _POWER_BOOTSTRAP_RESAMPLES
        or result.get("selected_seed_count") != 5
        or result.get("evaluator_id") != payload.get("endpoint_evaluator_id")
    ):
        raise TournamentError("finalist endpoint distribution contract changed")
    bootstrap_fields = (
        "candidate_primary_bootstrap",
        "baseline_primary_bootstrap",
        "paired_difference_bootstrap",
    )
    distributions: list[list[float]] = []
    for field in bootstrap_fields:
        raw_values = result.get(field)
        if not isinstance(raw_values, list) or len(raw_values) != _POWER_BOOTSTRAP_RESAMPLES:
            raise TournamentError(f"finalist {field} must contain exactly 10,000 draws")
        distributions.append(
            [_finite_number(value, f"{field}[]") for value in raw_values]
        )
    candidate_values, baseline_values, paired_values = distributions
    if any(
        abs((candidate - baseline) - paired) > 1e-12
        for candidate, baseline, paired in zip(
            candidate_values, baseline_values, paired_values, strict=True
        )
    ):
        raise TournamentError("finalist paired draws differ from candidate-baseline draws")
    for field, values in zip(bootstrap_fields, distributions, strict=True):
        error_field = {
            "candidate_primary_bootstrap": "candidate_primary_standard_error",
            "baseline_primary_bootstrap": "baseline_primary_standard_error",
            "paired_difference_bootstrap": "paired_effect_standard_error",
        }[field]
        if abs(_sample_standard_error(values) - _finite_number(result.get(error_field), error_field)) > 1e-12:
            raise TournamentError(f"finalist {error_field} does not rederive")
    candidate_metric = _finite_number(
        result.get("candidate_primary_metric"), "candidate_primary_metric"
    )
    baseline_metric = _finite_number(
        result.get("baseline_primary_metric"), "baseline_primary_metric"
    )
    if abs(
        (candidate_metric - baseline_metric)
        - _finite_number(result.get("observed_effect"), "observed_effect")
    ) > 1e-12:
        raise TournamentError("finalist observed effect does not rederive")
    for field in (
        "observed_effect",
        "candidate_primary_metric",
        "baseline_primary_metric",
        "n_units",
        "unit_set_sha256",
        "n_rows",
        "row_set_sha256",
    ):
        if payload.get(field) != result.get(field):
            raise TournamentError(f"finalist bundle/result disagree on {field}")
    identity = dict(payload)
    identity.pop("metric_bundle_id", None)
    if canonical_sha256(identity) != claimed:
        raise TournamentError("finalist development metric identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "finalist_development_metric_bundle",
        "metric_bundle_id": claimed,
    }:
        raise TournamentError("finalist development metric manifest metadata mismatch")
    if reverify_sources:
        lock_binding = payload.get("selection_lock_binding")
        ledger_binding = payload.get("selection_candidate_ledger_binding")
        if not isinstance(lock_binding, Mapping) or not isinstance(ledger_binding, Mapping):
            raise TournamentError("finalist metric source bindings are missing")
        lock_root = _verify_bound_frozen_document(
            {
                key: lock_binding.get(key)
                for key in ("path", "manifest_sha256", "document_sha256")
            },
            filename="selection_lock.json",
            artifact_class="selection_lock",
        )
        ledger_root = _verify_bound_frozen_document(
            {
                key: ledger_binding.get(key)
                for key in ("path", "manifest_sha256", "document_sha256")
            },
            filename="selection_candidate_ledger.json",
            artifact_class="selection_candidate_ledger",
        )
        sources = _finalist_metric_sources(
            selection_lock_dir=lock_root,
            selection_candidate_ledger_dir=ledger_root,
            task_id=str(payload.get("task_id")),
        )
        try:
            ensemble_ref = ArtifactRef.from_dict(payload.get("ensemble_table"))
        except (ContractError, TypeError) as exc:
            raise TournamentError(
                f"invalid development ensemble ArtifactRef: {exc}"
            ) from exc
        if (
            ensemble_ref.media_type != "text/tab-separated-values"
            or ensemble_ref.role
            != f"development_ensemble_endpoint_table:{payload.get('task_id')}"
        ):
            raise TournamentError("development ensemble ArtifactRef is invalid")
        ensemble_path = ensemble_ref.validate(root, require_relative=True)
        try:
            published_text = ensemble_path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError) as exc:
            raise TournamentError(
                f"cannot read the published development ensemble table: {exc}"
            ) from exc
        if published_text != sources["ensemble_table_text"]:
            raise TournamentError(
                "published development ensemble table does not rederive"
            )
        published_rows = _read_exact_tsv(
            ensemble_path,
            _JOINED_ENDPOINT_FIELDS[str(payload["task_id"])],
            "published development ensemble table",
        )
        if canonical_sha256(
            {
                "task_id": str(payload["task_id"]),
                "fields": list(_JOINED_ENDPOINT_FIELDS[str(payload["task_id"])]),
                "rows": [
                    [row[field] for field in _JOINED_ENDPOINT_FIELDS[
                        str(payload["task_id"])
                    ]]
                    for row in published_rows
                ],
            }
        ) != payload.get("ensemble_rows_sha256"):
            raise TournamentError(
                "published development ensemble rows do not rederive their hash"
            )
        try:
            republished = recompute_endpoint_bootstrap(
                evaluator_id=str(payload["endpoint_evaluator_id"]),
                table_path=ensemble_path,
                n_resamples=_POWER_BOOTSTRAP_RESAMPLES,
                seed=_POWER_BOOTSTRAP_SEED,
                parameters=payload["endpoint_parameters"],
            )
        except EndpointPowerError as exc:
            raise TournamentError(
                f"published development ensemble does not re-evaluate: {exc}"
            ) from exc
        if (
            abs(
                republished.candidate_primary_metric
                - _finite_number(
                    result["candidate_primary_metric"], "candidate_primary_metric"
                )
            )
            > 1e-12
            or abs(
                republished.observed_effect
                - _finite_number(result["observed_effect"], "observed_effect")
            )
            > 1e-12
        ):
            raise TournamentError(
                "published development ensemble endpoint differs from the frozen "
                "distribution"
            )
        expected = {
            "metric_bundle_id": claimed,
            **_finalist_metric_identity(
                sources=sources,
                endpoint_result_ref=result_ref,
                ensemble_table_ref=ensemble_ref,
            ),
        }
        if dict(payload) != expected or dict(result) != sources["endpoint_result"]:
            raise TournamentError(
                "finalist development metric bundle does not rederive from sources"
            )
    return dict(payload)


def _coerce_objectives(
    objectives: Mapping[str, Any] | Sequence[Objective],
) -> tuple[Objective, ...]:
    result: list[Objective] = []
    if isinstance(objectives, Mapping):
        for name, specification in objectives.items():
            maximize = True
            tolerance = 0.0
            if isinstance(specification, bool):
                maximize = specification
            elif isinstance(specification, str):
                direction = specification.strip().lower()
                if direction not in {"max", "maximize", "min", "minimize"}:
                    raise TournamentError(f"invalid direction for objective {name!r}")
                maximize = direction in {"max", "maximize"}
            elif isinstance(specification, Mapping):
                direction = specification.get(
                    "direction", specification.get("maximize", True)
                )
                if isinstance(direction, str):
                    direction = direction.strip().lower()
                    if direction not in {"max", "maximize", "min", "minimize"}:
                        raise TournamentError(f"invalid direction for objective {name!r}")
                    maximize = direction in {"max", "maximize"}
                elif isinstance(direction, bool):
                    maximize = direction
                else:
                    raise TournamentError(f"invalid direction for objective {name!r}")
                tolerance = _finite_number(
                    specification.get("tolerance", 0.0), f"{name}.tolerance"
                )
            else:
                raise TournamentError(f"invalid objective specification for {name!r}")
            result.append(Objective(str(name), maximize, tolerance))
    else:
        if isinstance(objectives, (str, bytes)):
            raise TournamentError("objectives must be a mapping or Objective sequence")
        for objective in objectives:
            if not isinstance(objective, Objective):
                raise TournamentError("objective sequences must contain Objective values")
            result.append(objective)
    if not result:
        raise TournamentError("at least one Pareto objective is required")
    names: set[str] = set()
    checked: list[Objective] = []
    for objective in result:
        name = _nonempty_identifier(objective.name, "objective.name")
        if name in names:
            raise TournamentError(f"duplicate Pareto objective: {name!r}")
        names.add(name)
        tolerance = _finite_number(objective.tolerance, f"{name}.tolerance")
        if tolerance < 0.0:
            raise TournamentError("objective tolerance cannot be negative")
        checked.append(Objective(name, bool(objective.maximize), tolerance))
    return tuple(checked)


def _dominates(candidate: Any, other: Any, objectives: Sequence[Objective]) -> bool:
    at_least_as_good = True
    strictly_better = False
    for objective in objectives:
        candidate_value = _candidate_metric(candidate, objective.name)
        other_value = _candidate_metric(other, objective.name)
        if objective.maximize:
            at_least_as_good &= candidate_value >= other_value - objective.tolerance
            strictly_better |= candidate_value > other_value + objective.tolerance
        else:
            at_least_as_good &= candidate_value <= other_value + objective.tolerance
            strictly_better |= candidate_value < other_value - objective.tolerance
    return at_least_as_good and strictly_better


def pareto_frontier(
    candidates: Iterable[Any],
    objectives: Mapping[str, Any] | Sequence[Objective],
    *,
    eligible_only: bool = True,
) -> tuple[str, ...]:
    records = _validated_candidates(candidates)
    parsed_objectives = _coerce_objectives(objectives)
    if eligible_only:
        records = [record for record in records if _candidate_eligible(record)]
    records.sort(key=_candidate_id)
    return tuple(
        _candidate_id(candidate)
        for candidate in records
        if not any(
            _dominates(other, candidate, parsed_objectives)
            for other in records
            if other is not candidate
        )
    )


def one_standard_error_candidates(
    candidates: Iterable[Any],
    *,
    primary_metric: str,
    higher_is_better: bool = True,
) -> tuple[str, ...]:
    records = [
        record for record in _validated_candidates(candidates) if _candidate_eligible(record)
    ]
    if not records:
        return ()
    records.sort(key=_candidate_id)
    if higher_is_better:
        best = min(
            records,
            key=lambda record: (
                -_candidate_metric(record, primary_metric),
                _candidate_id(record),
            ),
        )
        threshold = _candidate_metric(best, primary_metric) - _candidate_standard_error(
            best, primary_metric
        )
        selected = [
            record
            for record in records
            if _candidate_metric(record, primary_metric) >= threshold
        ]
        selected.sort(
            key=lambda record: (
                -_candidate_metric(record, primary_metric),
                _candidate_id(record),
            )
        )
    else:
        best = min(
            records,
            key=lambda record: (
                _candidate_metric(record, primary_metric),
                _candidate_id(record),
            ),
        )
        threshold = _candidate_metric(best, primary_metric) + _candidate_standard_error(
            best, primary_metric
        )
        selected = [
            record
            for record in records
            if _candidate_metric(record, primary_metric) <= threshold
        ]
        selected.sort(
            key=lambda record: (
                _candidate_metric(record, primary_metric),
                _candidate_id(record),
            )
        )
    return tuple(_candidate_id(record) for record in selected)


def _locked_string_array(value: Any) -> list[str] | None:
    """Read a fixed string array that may be a list or a requirements tuple.

    ``contracts._as_metadata`` freezes every nested array inside a
    ``SelectionLock`` task decision into a tuple, so ``isinstance(x, list)``
    and ``x != [...]`` are both False for perfectly valid fixed records.
    Returning a list lets the callers keep comparing values exactly while
    ignoring container type.
    """

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        return None
    items = list(value)
    if any(not isinstance(item, str) or not item for item in items):
        return None
    return items


def _stable_variant_primary_capability(
    task_id: str, decision: Mapping[str, Any]
) -> None:
    binding = decision.get("variant_primary_capability")
    if task_id != VARIANT_TASK:
        if binding is not None:
            raise SelectionLockError(
                f"{task_id}.variant_primary_capability must be null"
            )
        return
    root_fields = {
        "schema_version",
        "task_id",
        "primary_endpoint_id",
        "primary_evaluator_id",
        "capability_registry_sha256",
        "selected",
        "baseline",
    }
    if not isinstance(binding, Mapping) or set(binding) != root_fields:
        raise SelectionLockError(
            f"{task_id}.variant_primary_capability has an invalid schema"
        )
    expected = {
        "schema_version": "masld-bench-locked-variant-primary-capability-v1",
        "task_id": VARIANT_TASK,
        "primary_endpoint_id": "signed_cell_type_eqtl_effect",
        "primary_evaluator_id": "variant_ld_block_fisher_z_spearman_gain_v1",
    }
    if any(binding.get(field) != value for field, value in expected.items()):
        raise SelectionLockError(
            f"{task_id}.variant_primary_capability identity is not canonical"
        )
    _sha256_identifier(
        binding["capability_registry_sha256"],
        f"{task_id}.variant_primary_capability.capability_registry_sha256",
    )
    role_fields = {"model_id", "capability", "capability_sha256"}
    capability_fields = {
        "model_id",
        "role",
        "native_outputs",
        "allowed_endpoints",
        "primary_eligible",
        "requires_fitted_head",
        "requires_observed_target_context",
        "is_mandatory_baseline",
    }
    for role, model_field in (
        ("selected", "selected_model_id"),
        ("baseline", "baseline_model_id"),
    ):
        role_binding = binding.get(role)
        if not isinstance(role_binding, Mapping) or set(role_binding) != role_fields:
            raise SelectionLockError(
                f"{task_id} locked {role} capability binding is invalid"
            )
        model_id = _nonempty_identifier(
            role_binding["model_id"],
            f"{task_id}.variant_primary_capability.{role}.model_id",
        )
        if model_id != decision[model_field]:
            raise SelectionLockError(
                f"{task_id} locked {role} capability names the wrong model"
            )
        capability = role_binding.get("capability")
        if not isinstance(capability, Mapping) or set(capability) != capability_fields:
            raise SelectionLockError(
                f"{task_id} locked {role} capability record is invalid"
            )
        native_outputs = _locked_string_array(capability.get("native_outputs"))
        allowed_endpoints = _locked_string_array(
            capability.get("allowed_endpoints")
        )
        if (
            capability.get("role")
            not in {"baseline", "conditional_only", "primary_candidate"}
            or not native_outputs
            or native_outputs != sorted(set(native_outputs))
            or not allowed_endpoints
            or allowed_endpoints != sorted(set(allowed_endpoints))
            or any(
                not isinstance(capability.get(field), bool)
                for field in (
                    "primary_eligible",
                    "requires_fitted_head",
                    "requires_observed_target_context",
                    "is_mandatory_baseline",
                )
            )
        ):
            raise SelectionLockError(
                f"{task_id} locked {role} capability record is not canonical"
            )
        if (
            capability.get("model_id") != model_id
            or capability.get("primary_eligible") is not True
            or capability.get("requires_observed_target_context") is not False
            or "signed_cell_type_eqtl_effect"
            not in allowed_endpoints
            or not set(native_outputs).intersection(
                {"gene_expression_delta", "rna_coverage_delta"}
            )
        ):
            raise SelectionLockError(
                f"{task_id} locked {role} model cannot enter signed-effect scoring"
            )
        capability_sha256 = _sha256_identifier(
            role_binding["capability_sha256"],
            f"{task_id}.variant_primary_capability.{role}.capability_sha256",
        )
        if capability_sha256 != canonical_sha256(capability):
            raise SelectionLockError(
                f"{task_id} locked {role} capability SHA-256 changed"
            )


def _stable_variant_secondary_evaluation(
    task_id: str, decision: Mapping[str, Any]
) -> tuple[str, ...]:
    binding = decision.get("variant_secondary_evaluation")
    if task_id != VARIANT_TASK:
        if binding is not None:
            raise SelectionLockError(
                f"{task_id}.variant_secondary_evaluation must be null"
            )
        return ()
    root_fields = {
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
    if not isinstance(binding, Mapping) or set(binding) != root_fields:
        raise SelectionLockError(
            f"{task_id}.variant_secondary_evaluation has an invalid schema"
        )
    expected = {
        "schema_version": "masld-bench-locked-variant-secondary-evaluation-v1",
        "task_id": VARIANT_TASK,
        "endpoint_id": "eqtl_retrieval",
        "evaluator_id": "variant_ld_block_eqtl_retrieval_auprc_v1",
        "candidate_transform_id": "absolute_signed_effect_v1",
        "primary_baseline_transform_id": "absolute_signed_effect_v1",
    }
    if any(binding.get(field) != value for field, value in expected.items()):
        raise SelectionLockError(
            f"{task_id}.variant_secondary_evaluation identity is not canonical"
        )
    registry_sha256 = _sha256_identifier(
        binding["capability_registry_sha256"],
        f"{task_id}.variant_secondary_evaluation.capability_registry_sha256",
    )
    primary_binding = decision.get("variant_primary_capability")
    if (
        not isinstance(primary_binding, Mapping)
        or registry_sha256 != primary_binding.get("capability_registry_sha256")
    ):
        raise SelectionLockError(
            f"{task_id} primary and secondary capability registries differ"
        )
    mandatory_model_ids = _locked_string_array(binding.get("mandatory_model_ids"))
    expected_models = ["abc", "nearest_gene", "re2g"]
    expected_roles = {
        "abc": "link_only",
        "nearest_gene": "baseline",
        "re2g": "link_only",
    }
    if mandatory_model_ids != expected_models:
        raise SelectionLockError(
            f"{task_id} secondary evaluation lacks the exact mandatory roster"
        )
    raw_comparators = binding.get("comparators")
    comparators = (
        list(raw_comparators)
        if isinstance(raw_comparators, Sequence)
        and not isinstance(raw_comparators, (str, bytes))
        else None
    )
    if comparators is None or len(comparators) != len(expected_models):
        raise SelectionLockError(
            f"{task_id} secondary evaluation must bind each mandatory comparator once"
        )
    comparator_fields = {
        "model_id",
        "run_id",
        "capability",
        "capability_sha256",
        "score_transform_id",
    }
    capability_fields = {
        "model_id",
        "role",
        "native_outputs",
        "allowed_endpoints",
        "primary_eligible",
        "requires_fitted_head",
        "requires_observed_target_context",
        "is_mandatory_baseline",
    }
    observed_model_ids: list[str] = []
    run_ids: list[str] = []
    for index, comparator in enumerate(comparators):
        if not isinstance(comparator, Mapping) or set(comparator) != comparator_fields:
            raise SelectionLockError(
                f"{task_id} secondary comparator {index} has an invalid schema"
            )
        model_id = _nonempty_identifier(
            comparator["model_id"],
            f"{task_id}.variant_secondary_evaluation.comparators[{index}].model_id",
        )
        observed_model_ids.append(model_id)
        run_ids.append(
            _sha256_identifier(
                comparator["run_id"],
                f"{task_id}.variant_secondary_evaluation.comparators[{index}].run_id",
            )
        )
        if comparator.get("score_transform_id") != "identity_link_score_v1":
            raise SelectionLockError(
                f"{task_id} secondary comparator {model_id} score transform changed"
            )
        capability = comparator.get("capability")
        if not isinstance(capability, Mapping) or set(capability) != capability_fields:
            raise SelectionLockError(
                f"{task_id} secondary comparator {model_id} capability is invalid"
            )
        native_outputs = _locked_string_array(capability.get("native_outputs"))
        allowed_endpoints = _locked_string_array(
            capability.get("allowed_endpoints")
        )
        if (
            capability.get("model_id") != model_id
            or capability.get("role") != expected_roles.get(model_id)
            or native_outputs != ["enhancer_gene_link_score"]
            or allowed_endpoints != ["enhancer_gene_link", "eqtl_retrieval"]
            or any(
                not isinstance(capability.get(field), bool)
                for field in (
                    "primary_eligible",
                    "requires_fitted_head",
                    "requires_observed_target_context",
                    "is_mandatory_baseline",
                )
            )
            or capability.get("primary_eligible") is not False
            or capability.get("requires_fitted_head") is not False
            or capability.get("requires_observed_target_context") is not False
            or capability.get("is_mandatory_baseline") is not True
        ):
            raise SelectionLockError(
                f"{task_id} secondary comparator {model_id} is not an eligible link baseline"
            )
        capability_sha256 = _sha256_identifier(
            comparator["capability_sha256"],
            f"{task_id}.variant_secondary_evaluation.comparators[{index}].capability_sha256",
        )
        if capability_sha256 != canonical_sha256(capability):
            raise SelectionLockError(
                f"{task_id} secondary comparator {model_id} capability SHA-256 changed"
            )
    if observed_model_ids != expected_models or len(set(run_ids)) != len(run_ids):
        raise SelectionLockError(
            f"{task_id} secondary comparators are duplicated or not canonical"
        )
    primary_model_ids = {
        str(decision["selected_model_id"]), str(decision["baseline_model_id"])
    }
    if primary_model_ids.intersection(observed_model_ids):
        raise SelectionLockError(
            f"{task_id} signed models and secondary comparators must be distinct"
        )
    return tuple(run_ids)


def _stable_task_decision(
    lock: SelectionLock, decision: Mapping[str, Any]
) -> Mapping[str, Any]:
    task_id = str(decision.get("task_id", ""))
    if task_id not in EXPECTED_TASK_IDS:
        raise SelectionLockError(f"SelectionLock contains unsupported task_id {task_id!r}")
    if set(decision) != _LOCKED_TASK_DECISION_FIELDS:
        missing = sorted(_LOCKED_TASK_DECISION_FIELDS - set(decision))
        extra = sorted(set(decision) - _LOCKED_TASK_DECISION_FIELDS)
        raise SelectionLockError(
            f"task decision {task_id!r} does not match the stable schema: "
            f"missing={missing!r}, extra={extra!r}"
        )
    for field in (
        "selected_model_id",
        "baseline_model_id",
        "selected_adaptation_regime",
        "baseline_adaptation_regime",
        "promotion_gate",
    ):
        _nonempty_identifier(decision[field], f"{task_id}.{field}")
    if decision["selected_model_id"] in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
        raise SelectionLockError(
            f"{task_id} selects a negative control as its scientific model"
        )
    if decision["baseline_model_id"] in NON_SELECTABLE_NEGATIVE_CONTROL_IDS:
        raise SelectionLockError(
            f"{task_id} uses a negative control as its champion baseline"
        )
    if decision["selected_model_id"] == decision["baseline_model_id"]:
        raise SelectionLockError(
            f"{task_id} selected and baseline model IDs must be distinct"
        )
    for field in _TASK_HASH_FIELDS:
        _sha256_identifier(decision[field], f"{task_id}.{field}")
    if not isinstance(decision["open_champion"], bool):
        raise SelectionLockError(f"{task_id}.open_champion must be boolean")
    if (
        decision["open_champion"]
        and decision["selected_model_id"]
        in OPEN_CHAMPION_INELIGIBLE_CONTROL_IDS
    ):
        raise SelectionLockError(
            f"{task_id} gives an open champion claim to a control or ablation"
        )
    if not isinstance(decision["thresholds"], Mapping):
        raise SelectionLockError(f"{task_id}.thresholds must be a mapping")

    seeds = decision["seeds"]
    if not isinstance(seeds, Sequence) or isinstance(seeds, (str, bytes)):
        raise SelectionLockError(f"{task_id}.seeds must be an array")
    frozen_seeds = tuple(seeds)
    if (
        len(frozen_seeds) != 5
        or any(not isinstance(seed, int) or isinstance(seed, bool) for seed in frozen_seeds)
        or len(set(frozen_seeds)) != 5
        or frozen_seeds != tuple(sorted(frozen_seeds))
    ):
        raise SelectionLockError(f"{task_id}.seeds must be five unique sorted integers")

    run_ids_by_role: dict[str, tuple[str, ...]] = {}
    for role in ("selected", "baseline"):
        records = decision[f"{role}_runs"]
        if (
            not isinstance(records, Sequence)
            or isinstance(records, (str, bytes))
            or len(records) != 5
        ):
            raise SelectionLockError(
                f"{task_id}.{role}_runs must bind all five frozen seeds"
            )
        run_ids: list[str] = []
        for index, record in enumerate(records):
            if not isinstance(record, Mapping) or set(record) != {"seed", "run_id"}:
                raise SelectionLockError(
                    f"{task_id}.{role}_runs[{index}] must contain only seed and run_id"
                )
            if record["seed"] != frozen_seeds[index]:
                raise SelectionLockError(
                    f"{task_id}.{role}_runs are not aligned to frozen seeds"
                )
            run_ids.append(
                _sha256_identifier(
                    record["run_id"], f"{task_id}.{role}_runs[{index}].run_id"
                )
            )
        if len(set(run_ids)) != 5:
            raise SelectionLockError(f"{task_id}.{role}_runs repeat a run_id")
        run_ids_by_role[role] = tuple(run_ids)
    if set(run_ids_by_role["selected"]) & set(run_ids_by_role["baseline"]):
        raise SelectionLockError(
            f"{task_id} selected and strongest-baseline run sets must be disjoint"
        )
    if not set(run_ids_by_role["selected"]).issubset(lock.selected_run_ids):
        raise SelectionLockError(f"{task_id} selected runs are absent from selected_run_ids")
    if not set(run_ids_by_role["selected"] + run_ids_by_role["baseline"]).issubset(
        lock.candidate_run_ids
    ):
        raise SelectionLockError(f"{task_id} selected or baseline runs are not candidates")

    candidate_datasets = _unique_strings(
        decision["candidate_dataset_ids"],
        f"{task_id}.candidate_dataset_ids",
        allow_empty=True,
        sorted_required=True,
    )
    sealed_datasets = _unique_strings(
        decision["sealed_dataset_ids"],
        f"{task_id}.sealed_dataset_ids",
        allow_empty=True,
        sorted_required=True,
    )
    if not set(sealed_datasets).issubset(candidate_datasets):
        raise SelectionLockError(f"{task_id}.sealed_dataset_ids are not candidates")
    bindings = decision["dataset_registry_sha256s"]
    if not isinstance(bindings, Mapping) or set(bindings) != set(candidate_datasets):
        raise SelectionLockError(
            f"{task_id}.dataset_registry_sha256s must exactly cover candidate datasets"
        )
    for dataset_id, digest in bindings.items():
        _sha256_identifier(digest, f"{task_id}.dataset_registry_sha256s[{dataset_id}]")
    _stable_variant_primary_capability(task_id, decision)
    secondary_run_ids = _stable_variant_secondary_evaluation(task_id, decision)
    if not set(secondary_run_ids).issubset(lock.candidate_run_ids):
        raise SelectionLockError(
            f"{task_id} secondary comparator runs are not candidates"
        )
    if set(secondary_run_ids).intersection(
        run_ids_by_role["selected"] + run_ids_by_role["baseline"]
    ):
        raise SelectionLockError(
            f"{task_id} secondary comparator runs overlap signed-model runs"
        )
    return decision


def _strict_selection_lock(selection_lock: Any) -> SelectionLock:
    if isinstance(selection_lock, SelectionLock):
        lock = selection_lock
    elif isinstance(selection_lock, Mapping):
        try:
            lock = SelectionLock.from_dict(selection_lock)
        except (ContractError, TypeError, ValueError) as exc:
            raise SelectionLockError(
                f"selection lock mapping does not satisfy SelectionLock: {exc}"
            ) from exc
    else:
        raise SelectionLockError("selection_lock must be SelectionLock or its exact mapping")
    try:
        lock.require_locked(for_external_scoring=True)
    except (ContractError, TypeError, ValueError) as exc:
        raise SelectionLockError(str(exc)) from exc
    if lock.schema_version != SELECTION_LOCK_SCHEMA_VERSION:
        raise SelectionLockError("selection lock has the wrong schema_version")
    if lock.release_state is not ReleaseState.CANDIDATE:
        raise SelectionLockError("selection lock must have release_state=candidate")
    if lock.multiplicity_plan != LOCKED_MULTIPLICITY_PLAN:
        raise SelectionLockError("selection lock has the wrong multiplicity plan")
    if lock.terminal_policy != LOCKED_TERMINAL_POLICY:
        raise SelectionLockError("selection lock has the wrong terminal policy")
    identity = lock.to_dict()
    claimed_lock_id = identity.pop("lock_id")
    if canonical_sha256(identity) != claimed_lock_id:
        raise SelectionLockError("selection lock identity hash mismatch")
    if lock.candidate_run_ids != tuple(sorted(lock.candidate_run_ids)):
        raise SelectionLockError("candidate_run_ids must be canonically sorted")
    if lock.selected_run_ids != tuple(sorted(lock.selected_run_ids)):
        raise SelectionLockError("selected_run_ids must be canonically sorted")
    task_ids = tuple(str(item.get("task_id", "")) for item in lock.task_decisions)
    if task_ids != tuple(sorted(task_ids)):
        raise SelectionLockError("task_decisions must be canonically sorted by task_id")
    for decision in lock.task_decisions:
        _stable_task_decision(lock, decision)
    expected_selected = sorted(
        str(record["run_id"])
        for decision in lock.task_decisions
        for record in decision["selected_runs"]
    )
    if tuple(expected_selected) != lock.selected_run_ids:
        raise SelectionLockError(
            "aggregate selected_run_ids do not equal task-decision selected runs"
        )
    return lock


def require_selection_lock(
    selection_lock: Any, *, outcomes_must_be_locked: bool = True
) -> str:
    if outcomes_must_be_locked is not True:
        raise SelectionLockError(
            "promotion and development authorization require a pre-unblind SelectionLock"
        )
    return _strict_selection_lock(selection_lock).lock_id


def family_model_slots(
    candidates: Iterable[Any],
    *,
    primary_metric: str,
    higher_is_better: bool = True,
    max_per_family: int = 2,
) -> dict[str, tuple[str, ...]]:
    """Per family, the best candidate of each of the top ``max_per_family`` models.

    This is the auditable projection of the unique-model slot rule: two recipes
    of one model can never consume both slots, because the slots are chosen over
    models and only then resolved to that model's best recipe.  The Pareto and
    one-standard-error unions added by :func:`select_family_top_two` are
    deliberately NOT model-deduplicated, so the full shortlist may still contain
    two recipes of one model; those arrive on merit, not on a family slot.
    """

    if (
        not isinstance(max_per_family, int)
        or isinstance(max_per_family, bool)
        or max_per_family < 1
    ):
        raise TournamentError("max_per_family must be a positive integer")
    records = _validated_candidates(candidates)
    grouped: dict[str, list[Any]] = {}
    for record in records:
        if _candidate_eligible(record):
            grouped.setdefault(_candidate_family(record), []).append(record)

    def sort_key(record: Any) -> tuple[float, str]:
        value = _candidate_metric(record, primary_metric)
        return (-value if higher_is_better else value, _candidate_id(record))

    slots: dict[str, tuple[str, ...]] = {}
    for family in sorted(grouped):
        best_by_model: dict[str, Any] = {}
        for record in grouped[family]:
            model_id = _candidate_model_id(record)
            prior = best_by_model.get(model_id)
            if prior is None or sort_key(record) < sort_key(prior):
                best_by_model[model_id] = record
        ranked = sorted(
            best_by_model.values(),
            key=lambda record: (
                sort_key(record)[0],
                _candidate_model_id(record),
                _candidate_id(record),
            ),
        )
        slots[family] = tuple(
            _candidate_id(record) for record in ranked[:max_per_family]
        )
    return slots


def select_family_top_two(
    candidates: Iterable[Any],
    *,
    primary_metric: str,
    objectives: Mapping[str, Any] | Sequence[Objective] | None = None,
    higher_is_better: bool = True,
    max_per_family: int = 2,
) -> dict[str, tuple[str, ...]]:
    """Return the frozen screening union at model level, before seed expansion.

    The shortlist is the union of the best two models in every family, the
    global Pareto frontier, and every model within one standard error of the
    global primary-metric winner.  Five seeded runs are created only later by
    SelectionLock.
    """

    if (
        not isinstance(max_per_family, int)
        or isinstance(max_per_family, bool)
        or max_per_family < 1
    ):
        raise TournamentError("max_per_family must be a positive integer")
    records = _validated_candidates(candidates)
    if objectives is None:
        objectives = {primary_metric: "max" if higher_is_better else "min"}
    grouped: dict[str, list[Any]] = {}
    for record in records:
        if _candidate_eligible(record):
            grouped.setdefault(_candidate_family(record), []).append(record)
    global_frontier = set(pareto_frontier(records, objectives))
    global_one_se = set(
        one_standard_error_candidates(
            records,
            primary_metric=primary_metric,
            higher_is_better=higher_is_better,
        )
    )
    family_top = {
        candidate_id
        for family_ids in family_model_slots(
            records,
            primary_metric=primary_metric,
            higher_is_better=higher_is_better,
            max_per_family=max_per_family,
        ).values()
        for candidate_id in family_ids
    }
    selected_ids = family_top | global_frontier | global_one_se
    selections: dict[str, tuple[str, ...]] = {}
    for family in sorted(grouped):
        family_records = [
            record
            for record in grouped[family]
            if _candidate_id(record) in selected_ids
        ]
        family_records.sort(
            key=lambda record: (
                -_candidate_metric(record, primary_metric)
                if higher_is_better
                else _candidate_metric(record, primary_metric),
                _candidate_model_id(record),
                _candidate_id(record),
            )
        )
        selections[family] = tuple(_candidate_id(record) for record in family_records)
    return selections


def _shortlist_candidate(candidate: Any) -> dict[str, Any]:
    metrics = _read(candidate, "metrics", None)
    errors = _read(candidate, "standard_errors", None)
    if not isinstance(metrics, Mapping) or not isinstance(errors, Mapping):
        raise TournamentError(
            "development-shortlist candidates require metrics and standard_errors mappings"
        )
    normalized_metrics = {
        _nonempty_identifier(str(key), "metric name"): _finite_number(
            value, f"metric {key}"
        )
        for key, value in metrics.items()
    }
    normalized_errors = {
        _nonempty_identifier(str(key), "standard-error metric"): _finite_number(
            value, f"standard error {key}"
        )
        for key, value in errors.items()
    }
    if any(value < 0.0 for value in normalized_errors.values()):
        raise TournamentError("standard errors cannot be negative")
    task_id = _nonempty_identifier(_read(candidate, "task_id"), "task_id")
    model_id = _nonempty_identifier(_read(candidate, "model_id"), "model_id")
    regime = _nonempty_identifier(
        _read(candidate, "adaptation_regime"), "adaptation_regime"
    )
    configuration_sha256 = _sha256_identifier(
        _read(candidate, "candidate_configuration_sha256"),
        "candidate_configuration_sha256",
    )
    complexity_status = _nonempty_identifier(
        _read(candidate, "complexity_status"), "complexity_status"
    )
    return {
        "candidate_id": _candidate_id(candidate),
        "task_id": task_id,
        "model_id": model_id,
        "adaptation_regime": regime,
        "candidate_configuration_sha256": configuration_sha256,
        "family": _candidate_family(candidate),
        "metrics": dict(sorted(normalized_metrics.items())),
        "standard_errors": dict(sorted(normalized_errors.items())),
        "complexity_status": complexity_status,
        "eligible": _candidate_eligible(candidate),
    }


def freeze_development_shortlist(
    *,
    selection_candidate_ledger_dir: str | Path,
    task_id: str,
    source_wave: str,
    primary_metric: str,
    objectives: Mapping[str, Any] | Sequence[Objective],
    output_root: str | Path,
    higher_is_better: bool = True,
    max_per_family: int = 2,
) -> Path:
    """Freeze a pre-selection development leaderboard and its model-level union.

    ``source_wave`` must name a wave that is allowed to shortlist and must be
    exactly the wave the bound ledger was built from, at exactly that wave's
    prospectively frozen seed set.  A one-seed smoke ledger names no such wave
    and can therefore never produce a specialist finalist shortlist.
    """

    if task_id not in EXPECTED_TASK_IDS:
        raise TournamentError(f"unsupported development task {task_id!r}")
    contract = SHORTLIST_SOURCE_WAVE_CONTRACTS.get(
        _nonempty_identifier(source_wave, "source_wave")
    )
    if contract is None:
        raise TournamentError(
            f"wave {source_wave!r} may not create a development shortlist"
        )
    ledger_root = _safe_resolve(
        selection_candidate_ledger_dir, "selection candidate ledger"
    )
    ledger = verify_selection_candidate_ledger(ledger_root)
    ledger_binding, _ = _frozen_document_binding(
        ledger_root,
        filename="selection_candidate_ledger.json",
        artifact_class="selection_candidate_ledger",
    )
    ledger_binding = {**ledger_binding, "ledger_id": ledger["ledger_id"]}
    if list(ledger.get("campaign_universe_wave_ids") or []) != [source_wave]:
        raise TournamentError(
            "development shortlist ledger is not exactly the declared source wave"
        )
    expected_seeds = list(contract["expected_seeds"])
    task_seeds = [
        list(item.get("seeds") or [])
        for item in ledger.get("task_seed_sets") or []
        if isinstance(item, Mapping) and item.get("task_id") == task_id
    ]
    if task_seeds != [expected_seeds]:
        raise TournamentError(
            f"development shortlist requires the frozen {source_wave} seed set "
            f"{expected_seeds} for {task_id}"
        )
    primary_metric = _nonempty_identifier(primary_metric, "primary_metric")
    parsed_objectives = _coerce_objectives(objectives)
    records = tuple(
        sorted(
            (
                _shortlist_candidate(item)
                for item in ledger["model_candidates"]
                if item.get("task_id") == task_id
            ),
            key=lambda item: item["candidate_id"],
        )
    )
    if not records:
        raise TournamentError("development leaderboard cannot be empty")
    selected = select_family_top_two(
        records,
        primary_metric=primary_metric,
        objectives=parsed_objectives,
        higher_is_better=higher_is_better,
        max_per_family=max_per_family,
    )
    selected_ids = sorted(
        identifier for family_ids in selected.values() for identifier in family_ids
    )
    # Reject BEFORE anything read-only is written.  Publishing an invalid target
    # first would strand that identity permanently, because a frozen directory
    # can never be replaced or deleted.
    if not selected_ids:
        raise TournamentError(
            "development shortlist selected no eligible candidates"
        )
    identity = {
        "schema_version": DEVELOPMENT_SHORTLIST_SCHEMA_VERSION,
        "task_id": task_id,
        "source_wave": source_wave,
        "authorized_downstream_wave": contract["authorizes_wave"],
        "source_task_seeds": expected_seeds,
        "family_model_slots": {
            family: list(ids)
            for family, ids in sorted(
                family_model_slots(
                    records,
                    primary_metric=primary_metric,
                    higher_is_better=higher_is_better,
                    max_per_family=max_per_family,
                ).items()
            )
        },
        "source_family_ids": sorted({record["family"] for record in records}),
        "selection_candidate_ledger_binding": ledger_binding,
        "created_before_outcome_unblind": True,
        "selection_rule": (
            "wave_bound_top_two_unique_models_per_family_union_pareto_union_"
            "one_standard_error_v3"
        ),
        "primary_metric": primary_metric,
        "higher_is_better": higher_is_better,
        "max_per_family": max_per_family,
        "objectives": [asdict(item) for item in parsed_objectives],
        "candidates": list(records),
        "selected_candidate_ids": selected_ids,
    }
    payload = {"shortlist_id": canonical_sha256(identity), **identity}
    root = _safe_resolve(output_root, "development shortlist output root")
    root.mkdir(parents=True, exist_ok=True)
    target = root / f"development-shortlist--{payload['shortlist_id']}"
    if target.exists():
        raise TournamentError(f"development shortlist already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".development-shortlist.", dir=root))
    try:
        write_json_exclusive(
            staging / "development_shortlist.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "development_shortlist",
                "shortlist_id": payload["shortlist_id"],
            },
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_development_shortlist(target)
    return target


def verify_development_shortlist(path: str | Path) -> dict[str, Any]:
    root = _safe_resolve(path, "development shortlist")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "development_shortlist.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid frozen development shortlist: {exc}") from exc
    fields = {
        "shortlist_id",
        "schema_version",
        "task_id",
        "source_wave",
        "authorized_downstream_wave",
        "source_task_seeds",
        "family_model_slots",
        "source_family_ids",
        "selection_candidate_ledger_binding",
        "created_before_outcome_unblind",
        "selection_rule",
        "primary_metric",
        "higher_is_better",
        "max_per_family",
        "objectives",
        "candidates",
        "selected_candidate_ids",
    }
    if not isinstance(payload, Mapping) or set(payload) != fields:
        raise TournamentError("development shortlist has the wrong schema")
    if payload["schema_version"] != DEVELOPMENT_SHORTLIST_SCHEMA_VERSION:
        raise TournamentError("unsupported development shortlist schema_version")
    identity = dict(payload)
    claimed = _sha256_identifier(identity.pop("shortlist_id"), "shortlist_id")
    if canonical_sha256(identity) != claimed:
        raise TournamentError("development shortlist identity hash mismatch")
    if payload["created_before_outcome_unblind"] is not True:
        raise TournamentError("development shortlist is post-unblind")
    if payload["task_id"] not in EXPECTED_TASK_IDS:
        raise TournamentError("development shortlist has an unsupported task_id")
    source_families = _unique_strings(
        payload["source_family_ids"], "source_family_ids", sorted_required=True
    )
    binding = payload["selection_candidate_ledger_binding"]
    if not isinstance(binding, Mapping) or set(binding) != {
        "path",
        "manifest_sha256",
        "document_sha256",
        "ledger_id",
    }:
        raise TournamentError("development shortlist ledger binding has the wrong schema")
    ledger_root = _verify_bound_frozen_document(
        {key: binding[key] for key in ("path", "manifest_sha256", "document_sha256")},
        filename="selection_candidate_ledger.json",
        artifact_class="selection_candidate_ledger",
    )
    ledger = verify_selection_candidate_ledger(ledger_root)
    if binding["ledger_id"] != ledger["ledger_id"]:
        raise TournamentError("development shortlist ledger binding changed")
    _nonempty_identifier(payload["primary_metric"], "primary_metric")
    if not isinstance(payload["higher_is_better"], bool):
        raise TournamentError("higher_is_better must be boolean")
    if (
        not isinstance(payload["max_per_family"], int)
        or isinstance(payload["max_per_family"], bool)
        or payload["max_per_family"] < 1
    ):
        raise TournamentError("max_per_family must be a positive integer")
    if (
        payload["selection_rule"]
        != (
            "wave_bound_top_two_unique_models_per_family_union_pareto_union_"
            "one_standard_error_v3"
        )
    ):
        raise TournamentError("development shortlist has the wrong selection rule")
    candidates = payload["candidates"]
    if not isinstance(candidates, list) or not candidates:
        raise TournamentError("development shortlist candidates must be non-empty")
    if any(not isinstance(item, Mapping) for item in candidates):
        raise TournamentError("development shortlist candidates must be objects")
    if [_shortlist_candidate(item) for item in candidates] != candidates:
        raise TournamentError("development shortlist candidates are not normalized")
    ledger_candidates = sorted(
        (
            _shortlist_candidate(item)
            for item in ledger["model_candidates"]
            if item.get("task_id") == payload["task_id"]
        ),
        key=lambda item: item["candidate_id"],
    )
    if candidates != ledger_candidates:
        raise TournamentError("development shortlist candidates changed from its ledger")
    if list(source_families) != sorted({item["family"] for item in candidates}):
        raise TournamentError("development shortlist source families do not rederive")
    if [item.get("candidate_id") for item in candidates] != sorted(
        item.get("candidate_id") for item in candidates
    ):
        raise TournamentError("development candidates must be canonically sorted")
    objectives_raw = payload["objectives"]
    if not isinstance(objectives_raw, list):
        raise TournamentError("development objectives must be an array")
    parsed_objectives = tuple(
        Objective(
            name=item["name"],
            maximize=item["maximize"],
            tolerance=item["tolerance"],
        )
        for item in objectives_raw
        if isinstance(item, Mapping) and set(item) == {"name", "maximize", "tolerance"}
    )
    if len(parsed_objectives) != len(objectives_raw):
        raise TournamentError("development objectives have the wrong schema")
    selected = select_family_top_two(
        candidates,
        primary_metric=payload["primary_metric"],
        objectives=parsed_objectives,
        higher_is_better=payload["higher_is_better"],
        max_per_family=payload["max_per_family"],
    )
    expected = sorted(item for values in selected.values() for item in values)
    if not expected:
        raise TournamentError("development shortlist selected no eligible candidates")
    if payload["selected_candidate_ids"] != expected:
        raise TournamentError("development shortlist selection does not rederive")
    source_wave = _nonempty_identifier(payload["source_wave"], "source_wave")
    contract = SHORTLIST_SOURCE_WAVE_CONTRACTS.get(source_wave)
    if contract is None:
        raise TournamentError(
            f"wave {source_wave!r} may not create a development shortlist"
        )
    if (
        payload["authorized_downstream_wave"] != contract["authorizes_wave"]
        or list(payload["source_task_seeds"]) != list(contract["expected_seeds"])
    ):
        raise TournamentError("development shortlist wave contract changed")
    if list(ledger.get("campaign_universe_wave_ids") or []) != [source_wave]:
        raise TournamentError(
            "development shortlist ledger is not exactly the declared source wave"
        )
    if [
        list(item.get("seeds") or [])
        for item in ledger.get("task_seed_sets") or []
        if isinstance(item, Mapping) and item.get("task_id") == payload["task_id"]
    ] != [list(contract["expected_seeds"])]:
        raise TournamentError(
            "development shortlist seed authority does not rederive from its ledger"
        )
    expected_slots = {
        family: list(ids)
        for family, ids in sorted(
            family_model_slots(
                candidates,
                primary_metric=payload["primary_metric"],
                higher_is_better=payload["higher_is_better"],
                max_per_family=payload["max_per_family"],
            ).items()
        )
    }
    if payload["family_model_slots"] != expected_slots:
        raise TournamentError(
            "development shortlist family model slots do not rederive"
        )
    metadata = manifest.get("metadata")
    if metadata != {"artifact_class": "development_shortlist", "shortlist_id": claimed}:
        raise TournamentError("development shortlist manifest metadata mismatch")
    return dict(payload)


def _join_row_unit_ids(
    rows: Sequence[Mapping[str, Any]],
    *,
    row_id_field: str,
    unit_id_field: str,
    label: str,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[tuple[str, str], ...]]:
    if not rows:
        raise TournamentError(f"{label} contains no rows")
    row_ids: list[str] = []
    unit_ids: list[str] = []
    for index, row in enumerate(rows):
        if row_id_field not in row or unit_id_field not in row:
            raise TournamentError(
                f"{label} row {index} lacks {row_id_field!r} or {unit_id_field!r}"
            )
        row_id = _sha256_identifier(row[row_id_field], f"{label} row ID")
        unit_id = _sha256_identifier(row[unit_id_field], f"{label} unit ID")
        row_ids.append(row_id)
        unit_ids.append(unit_id)
    if len(set(row_ids)) != len(row_ids):
        raise TournamentError(f"{label} repeats a row ID")
    if row_ids != sorted(row_ids):
        raise TournamentError(f"{label} row IDs are not canonically sorted")
    pairs = tuple(sorted(zip(row_ids, unit_ids, strict=True)))
    return tuple(row_ids), tuple(sorted(set(unit_ids))), pairs


def _prediction_join_ids(
    commit: PredictionCommit,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[tuple[str, str], ...]]:
    bundle_path = _safe_resolve(
        commit.prediction_bundle_path, "PredictionCommit PredictionBundle"
    )
    try:
        bundle = PredictionBundle.load_json(bundle_path)
        bundle.validate_artifacts(bundle_path.parent)
    except (ContractError, OSError, ValueError) as exc:
        raise TournamentError(f"invalid PredictionBundle for {commit.run_id}: {exc}") from exc
    if (
        bundle.standardized_table.media_type != "text/tab-separated-values"
        or bundle.row_ids.media_type != "text/tab-separated-values"
    ):
        raise TournamentError("PredictionBundle exact joins require canonical TSV artifacts")
    table_fields = _PREDICTION_TABLE_FIELDS.get(commit.task_id)
    if table_fields is None:
        raise TournamentError(f"PredictionBundle task {commit.task_id} has no row schema")
    table_rows = _read_exact_tsv(
        bundle.standardized_table.validate(bundle_path.parent, require_relative=True),
        table_fields,
        f"PredictionBundle {bundle.bundle_id} standardized table",
    )
    inventory_rows = _read_exact_tsv(
        bundle.row_ids.validate(bundle_path.parent, require_relative=True),
        (commit.row_id_field, commit.unit_id_field),
        f"PredictionBundle {bundle.bundle_id} row inventory",
    )
    table_identity = _join_row_unit_ids(
        table_rows,
        row_id_field=commit.row_id_field,
        unit_id_field=commit.unit_id_field,
        label=f"PredictionBundle {bundle.bundle_id} standardized table",
    )
    inventory_identity = _join_row_unit_ids(
        inventory_rows,
        row_id_field=commit.row_id_field,
        unit_id_field=commit.unit_id_field,
        label=f"PredictionBundle {bundle.bundle_id} row inventory",
    )
    if table_identity != inventory_identity:
        raise TournamentError("PredictionBundle row inventory differs from standardized table")
    if len(table_identity[0]) != commit.n_predictions:
        raise TournamentError(
            f"PredictionBundle {bundle.bundle_id} declares {commit.n_predictions} predictions "
            f"but contains {len(table_identity[0])} rows"
        )
    return table_identity


def _outcome_join_ids(
    *,
    outcome_root: Path,
    outcome_bundle: Any,
    task_id: str,
    row_id_field: str,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[tuple[str, str], ...]]:
    role = f"outcome:{task_id}"
    artifacts = tuple(
        artifact for artifact in outcome_bundle.artifacts if artifact.role == role
    )
    if not artifacts:
        raise TournamentError(f"sealed outcome bundle lacks {role!r} artifacts")
    unit_id_field = _nonempty_identifier(
        outcome_bundle.unit_id_fields.get(task_id), f"outcome unit field for {task_id}"
    )
    if len(artifacts) != 1 or artifacts[0].media_type != "text/tab-separated-values":
        raise TournamentError(
            f"sealed outcome {task_id} requires one canonical TSV outcome artifact"
        )
    artifact_path = artifacts[0].validate(outcome_root, require_relative=True)
    try:
        with artifact_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = tuple(reader.fieldnames or ())
            if (
                len(fields) != len(set(fields))
                or row_id_field not in fields
                or unit_id_field not in fields
            ):
                raise TournamentError(f"sealed outcome {task_id} has an invalid header")
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise TournamentError(f"cannot read sealed outcome {task_id}: {exc}") from exc
    return _join_row_unit_ids(
        rows,
        row_id_field=_nonempty_identifier(row_id_field, f"outcome row field for {task_id}"),
        unit_id_field=unit_id_field,
        label=f"sealed outcome {task_id}",
    )


def _verify_exact_task_units(
    *,
    commits: Sequence[PredictionCommit],
    outcome_root: Path | None = None,
    outcome_bundle: Any | None = None,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[tuple[str, str], ...]]:
    if not commits:
        raise TournamentError("unit-set verification requires PredictionCommits")
    task_ids = {commit.task_id for commit in commits}
    if len(task_ids) != 1:
        raise TournamentError("unit-set verification cannot mix tasks")
    prediction_sets = {
        commit.run_id: _prediction_join_ids(commit) for commit in commits
    }
    first = next(iter(prediction_sets.values()))
    if any(unit_ids != first for unit_ids in prediction_sets.values()):
        raise TournamentError(
            f"PredictionBundles for {next(iter(task_ids))} do not have an identical unit-ID set"
        )
    if any(commit.n_predictions != len(first[0]) for commit in commits):
        raise TournamentError("PredictionCommit counts do not equal the verified row-ID set")
    if (outcome_root is None) != (outcome_bundle is None):
        raise TournamentError("outcome root and bundle must be supplied together")
    if outcome_root is not None and outcome_bundle is not None:
        task_id = next(iter(task_ids))
        outcome_ids = _outcome_join_ids(
            outcome_root=outcome_root,
            outcome_bundle=outcome_bundle,
            task_id=task_id,
            row_id_field=commits[0].row_id_field,
        )
        if outcome_ids != first:
            raise TournamentError(
                f"PredictionBundle and sealed outcome row/unit-ID sets differ for {task_id}"
            )
    return first


def _frozen_document_binding(
    path: str | Path,
    *,
    filename: str,
    artifact_class: str,
) -> tuple[dict[str, str], dict[str, Any]]:
    root = reject_symlink_components(Path(path), label=f"{artifact_class} path")
    if not root.is_absolute() or root.is_symlink() or not root.is_dir():
        raise TournamentError(
            f"{artifact_class} path must be an absolute non-symlink directory"
        )
    root = root.resolve()
    try:
        manifest = verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid frozen {artifact_class}: {exc}") from exc
    metadata = manifest.get("metadata")
    if not isinstance(metadata, Mapping) or metadata.get("artifact_class") != artifact_class:
        raise TournamentError(f"frozen source is not classified as {artifact_class}")
    document = root / filename
    if document.is_symlink() or not document.is_file():
        raise TournamentError(f"frozen {artifact_class} lacks {filename}")
    return (
        {
            "path": root.as_posix(),
            "manifest_sha256": sha256_file(root / "ARTIFACTS.json"),
            "document_sha256": sha256_file(document),
        },
        manifest,
    )


def _verify_bound_frozen_document(
    binding: Any,
    *,
    filename: str,
    artifact_class: str,
) -> Path:
    if not isinstance(binding, Mapping) or set(binding) != {
        "path",
        "manifest_sha256",
        "document_sha256",
    }:
        raise TournamentError(f"{artifact_class} binding has the wrong schema")
    source = Path(_nonempty_identifier(binding["path"], f"{artifact_class}.path"))
    observed, _ = _frozen_document_binding(
        source, filename=filename, artifact_class=artifact_class
    )
    normalized = {
        "path": observed["path"],
        "manifest_sha256": _sha256_identifier(
            binding["manifest_sha256"], f"{artifact_class}.manifest_sha256"
        ),
        "document_sha256": _sha256_identifier(
            binding["document_sha256"], f"{artifact_class}.document_sha256"
        ),
    }
    if observed != normalized:
        raise TournamentError(f"{artifact_class} binding changed")
    return _safe_resolve(source, f"bound {artifact_class}")


def _prediction_commit_binding(
    path: str | Path,
) -> tuple[dict[str, Any], PredictionCommit, Path]:
    root = _safe_resolve(path, "PredictionCommit")
    try:
        commit = verify_prediction_commit(root, reverify_sources=True)
    except (FirewallError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid frozen PredictionCommit: {exc}") from exc
    binding, _ = _frozen_document_binding(
        root,
        filename="prediction_commit.json",
        artifact_class="prediction_commit",
    )
    return (
        {
            **binding,
            "prediction_commit_sha256": commit.prediction_commit_sha256,
            "prediction_bundle_sha256": commit.prediction_bundle_sha256,
            "run_id": commit.run_id,
            "model_role": commit.model_role,
        },
        commit,
        root,
    )


def _verified_task_prediction_commits(
    paths: Iterable[str | Path],
    *,
    lock: SelectionLock,
    task_id: str,
    allow_empty: bool,
) -> tuple[tuple[dict[str, Any], ...], tuple[PredictionCommit, ...], tuple[Path, ...]]:
    results = tuple(_prediction_commit_binding(path) for path in paths)
    if not results and not allow_empty:
        raise TournamentError(f"{task_id} requires frozen PredictionCommits")
    bindings = tuple(sorted((item[0] for item in results), key=lambda item: item["run_id"]))
    commits_by_run = {item[1].run_id: item[1] for item in results}
    roots_by_run = {item[1].run_id: item[2] for item in results}
    if len(commits_by_run) != len(results):
        raise TournamentError(f"{task_id} repeats a PredictionCommit run")
    commits = tuple(commits_by_run[binding["run_id"]] for binding in bindings)
    roots = tuple(roots_by_run[binding["run_id"]] for binding in bindings)
    for commit in commits:
        if commit.selection_lock_id != lock.lock_id or commit.task_id != task_id:
            raise TournamentError(
                f"PredictionCommit {commit.run_id} does not bind {task_id}/{lock.lock_id}"
            )
    return bindings, commits, roots


def _required_commit_runs(
    *, lock: SelectionLock, task_id: str
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    task_decision = _task_decision(lock, task_id)
    return (
        tuple(sorted(_locked_run_ids(task_decision, "selected"))),
        tuple(sorted(_locked_run_ids(task_decision, "baseline"))),
    )


def _secondary_commit_runs(
    *, lock: SelectionLock, task_id: str
) -> tuple[tuple[str, str], ...]:
    decision = _task_decision(lock, task_id)
    binding = decision.get("variant_secondary_evaluation")
    if task_id != VARIANT_TASK:
        return ()
    if not isinstance(binding, Mapping):
        raise TournamentError(
            "variant task lacks a locked secondary-evaluation binding"
        )
    return tuple(
        (str(comparator["run_id"]), str(comparator["model_id"]))
        for comparator in binding["comparators"]
    )


def _assert_prediction_commit_coverage(
    *,
    lock: SelectionLock,
    task_id: str,
    commits: Sequence[PredictionCommit],
    required: str,
) -> None:
    selected, baseline = _required_commit_runs(lock=lock, task_id=task_id)
    secondary = _secondary_commit_runs(lock=lock, task_id=task_id)
    secondary_run_ids = tuple(run_id for run_id, _model_id in secondary)
    observed = tuple(sorted(commit.run_id for commit in commits))
    if required == "selected":
        expected = selected
    elif required == "selected_and_baseline":
        expected = tuple(sorted((*selected, *baseline)))
    elif required == "selected_baseline_secondary":
        expected = tuple(sorted((*selected, *baseline, *secondary_run_ids)))
    else:
        raise TournamentError(f"unsupported prediction coverage mode {required!r}")
    if observed != expected:
        raise TournamentError(
            f"{task_id} PredictionCommits must exactly cover {required.replace('_', ' ')} runs"
        )
    task_decision = _task_decision(lock, task_id)
    expected_models = {
        **{run_id: str(task_decision["selected_model_id"]) for run_id in selected},
        **{run_id: str(task_decision["baseline_model_id"]) for run_id in baseline},
        **{run_id: model_id for run_id, model_id in secondary},
    }
    for commit in commits:
        if commit.run_id in selected:
            expected_role = "selected"
        elif commit.run_id in baseline:
            expected_role = "baseline"
        else:
            expected_role = "secondary_comparator"
        if commit.model_id != expected_models[commit.run_id] or commit.model_role != expected_role:
            raise TournamentError(
                f"PredictionCommit {commit.run_id} has the wrong locked model role"
            )


def _task_decision(lock: SelectionLock, task_id: str) -> Mapping[str, Any]:
    if task_id not in EXPECTED_TASK_IDS:
        raise TournamentError(
            f"unsupported task_id {task_id!r}; expected one of {sorted(EXPECTED_TASK_IDS)!r}"
        )
    matches = [item for item in lock.task_decisions if item.get("task_id") == task_id]
    if len(matches) != 1:
        raise SelectionLockError(
            f"selection lock must contain exactly one task decision for {task_id!r}"
        )
    return _stable_task_decision(lock, matches[0])


def _locked_run_ids(decision: Mapping[str, Any], role: str) -> tuple[str, ...]:
    return tuple(str(record["run_id"]) for record in decision[f"{role}_runs"])


def _exact_keys(value: Any, expected: set[str], name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != expected:
        observed = sorted(value) if isinstance(value, Mapping) else []
        raise SelectionLockError(
            f"{name} must contain exactly {sorted(expected)!r}; observed={observed!r}"
        )
    return value


def _exact_number(mapping: Mapping[str, Any], name: str, expected: Any) -> float:
    value = _finite_number(mapping.get(name), name)
    if value != float(expected):
        raise SelectionLockError(f"{name} must remain exactly {float(expected)}")
    return value


def _power_spec(
    promotion: Mapping[str, Any], policy: Mapping[str, Any], task_id: str
) -> Mapping[str, Any]:
    fields = {
        "metric_id",
        "minimum_effect",
        "alpha",
        "target_power",
        "n_primary_claims",
        "min_units",
        "two_sided",
    }
    power = _exact_keys(promotion.get("power"), fields, f"{task_id}.promotion.power")
    if power["metric_id"] != policy["power_metric"]:
        raise SelectionLockError(f"{task_id} power metric differs from the gate registry")
    _exact_number(power, "minimum_effect", policy["power_minimum_effect"])
    _exact_number(power, "alpha", 0.05)
    _exact_number(power, "target_power", 0.80)
    for field in ("n_primary_claims", "min_units"):
        value = power[field]
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise SelectionLockError(f"{task_id}.promotion.power.{field} must be positive")
    if power["min_units"] < 2:
        raise SelectionLockError(f"{task_id}.promotion.power.min_units must be at least two")
    if power["two_sided"] is not True:
        raise SelectionLockError(f"{task_id} promotion power must be two-sided")
    return power


def _conditional_spec(promotion: Mapping[str, Any], task_id: str) -> None:
    activation = _exact_keys(
        promotion.get("conditional_activation"),
        {"activated", "activation_receipt_sha256"},
        f"{task_id}.promotion.conditional_activation",
    )
    if activation["activated"] is not True:
        raise SelectionLockError(f"{task_id} conditional seal is not activated")
    _sha256_identifier(
        activation["activation_receipt_sha256"],
        f"{task_id}.promotion.activation_receipt_sha256",
    )


def _promotion_contract(
    decision: Mapping[str, Any], task_id: str
) -> Mapping[str, Any]:
    policy = CHAMPION_GATE_THRESHOLDS[task_id]
    if decision["promotion_gate"] != policy["promotion_gate_id"]:
        raise SelectionLockError(
            f"{task_id}.promotion_gate differs from the executable gate registry"
        )
    promotion = decision["thresholds"].get("promotion")
    if not isinstance(promotion, Mapping):
        raise SelectionLockError(f"{task_id}.thresholds.promotion must be a mapping")
    mode = _nonempty_identifier(promotion.get("promotion_mode"), "promotion_mode")
    allowed_modes = _unique_strings(policy["allowed_promotion_modes"], "allowed modes")
    if mode not in allowed_modes:
        raise SelectionLockError(f"promotion_mode {mode!r} is not allowed for {task_id}")
    promotable = mode in PROMOTABLE_MODES
    if decision["open_champion"] is not promotable:
        raise SelectionLockError(
            f"{task_id}.open_champion must exactly match its promotion mode"
        )
    sealed_ids = tuple(decision["sealed_dataset_ids"])
    expected_holdouts = tuple(policy["sealed_holdout_dataset_ids"])
    if promotable:
        if task_id == GRAPH_TASK:
            if not sealed_ids:
                raise SelectionLockError("a graph seal requires a frozen holdout dataset")
        elif sealed_ids != expected_holdouts:
            raise SelectionLockError(
                f"{task_id} sealed datasets must equal {expected_holdouts!r}"
            )
    elif not set(sealed_ids).issubset(expected_holdouts):
        raise SelectionLockError(
            f"non-promotable {task_id} has a seal outside its frozen TaskSpec gate"
        )

    numeric_fields: tuple[str, ...] = ()
    roster_field: str | None = None
    required_count: int | None = None
    extra_fields: set[str] = set()
    if task_id == CELL_TASK:
        numeric_fields = (
            "primary_min",
            "delta_min",
            "delta_ci_low_strict_min",
            "stratum_min",
            "brier_delta_max",
        )
        roster_field = "required_class_ids"
    elif task_id == VARIANT_TASK:
        numeric_fields = (
            "delta_min",
            "delta_ci_low_strict_min",
            "stratum_min",
        )
        roster_field = "required_major_lineage_ids"
        extra_fields.add("eqtl_auprc_noninferiority_margin")
    elif task_id == RNA_ATAC_TASK:
        numeric_fields = (
            "delta_min",
            "delta_ci_low_strict_min",
            "stratum_min",
            "min_improved_strata",
        )
        roster_field = "required_major_lineage_ids"
        required_count = int(policy["required_stratum_count"])
    elif task_id == BULK_TASK:
        numeric_fields = ("delta_min", "delta_ci_low_strict_min")
        roster_field = "required_endpoint_ids"
    elif task_id == GRAPH_TASK and promotable:
        numeric_fields = (
            "delta_min",
            "delta_ci_low_strict_min",
            "brier_delta_max",
        )
        extra_fields.add("sealed_source_available")
    elif task_id == UNIFIED_TASK and promotable:
        numeric_fields = (
            "min_superior_source_families",
            "min_external_superior_source_families",
            "min_eligible_external_source_families",
        )
        extra_fields.update(
            {
                "eligible_specialist_ids",
                "specialist_noninferiority_margins",
                "eligible_source_family_ids",
                "eligible_external_source_family_ids",
                "source_independence_group_by_family",
            }
        )

    fields = {"promotion_mode", *numeric_fields, *extra_fields}
    if roster_field is not None:
        fields.add(roster_field)
    if promotable:
        fields.update(
            {
                "confirmatory_family_id",
                "confirmatory_fwer",
                "sealed_outcome_bundle_id",
                "power",
            }
        )
    if mode == CONDITIONAL_SEALED:
        fields.add("conditional_activation")
    _exact_keys(promotion, fields, f"{task_id}.thresholds.promotion")

    for field in numeric_fields:
        _exact_number(promotion, field, policy[field])
    if roster_field is not None:
        roster = _unique_strings(promotion[roster_field], f"{task_id}.{roster_field}")
        if required_count is not None and len(roster) != required_count:
            raise SelectionLockError(
                f"{task_id} requires exactly {required_count} frozen major lineages"
            )
    if task_id == VARIANT_TASK:
        margin = _finite_number(
            promotion["eqtl_auprc_noninferiority_margin"],
            "eqtl_auprc_noninferiority_margin",
        )
        if margin < 0.0:
            raise SelectionLockError("AUPRC noninferiority margin cannot be negative")
    if promotable:
        _nonempty_identifier(
            promotion["confirmatory_family_id"], "confirmatory_family_id"
        )
        _exact_number(
            promotion,
            "confirmatory_fwer",
            PROMOTION_GATE_POLICY["confirmatory_fwer"],
        )
        _nonempty_identifier(
            promotion["sealed_outcome_bundle_id"], "sealed_outcome_bundle_id"
        )
        _power_spec(promotion, policy, task_id)
    if mode == CONDITIONAL_SEALED:
        _conditional_spec(promotion, task_id)
    if task_id == GRAPH_TASK and promotable:
        if promotion["sealed_source_available"] is not True:
            raise SelectionLockError("typed evidence graph lacks a wholly sealed source")
    if task_id == UNIFIED_TASK and promotable:
        specialists = _unique_strings(
            promotion["eligible_specialist_ids"], "eligible_specialist_ids"
        )
        margins = promotion["specialist_noninferiority_margins"]
        if not isinstance(margins, Mapping) or set(margins) != set(specialists):
            raise SelectionLockError(
                "specialist_noninferiority_margins must exactly cover eligible specialists"
            )
        if any(
            _finite_number(value, f"specialist margin {key}") < 0.0
            for key, value in margins.items()
        ):
            raise SelectionLockError("specialist noninferiority margins cannot be negative")
        source_ids = _unique_strings(
            promotion["eligible_source_family_ids"], "eligible_source_family_ids"
        )
        external_ids = _unique_strings(
            promotion["eligible_external_source_family_ids"],
            "eligible_external_source_family_ids",
        )
        if not set(external_ids).issubset(source_ids):
            raise SelectionLockError("external source families must be eligible sources")
        groups = promotion["source_independence_group_by_family"]
        if not isinstance(groups, Mapping) or set(groups) != set(source_ids):
            raise SelectionLockError(
                "source_independence_group_by_family must exactly cover eligible sources"
            )
        normalized_groups = {
            key: _nonempty_identifier(value, f"source group {key}")
            for key, value in groups.items()
        }
        required_external = int(policy["min_eligible_external_source_families"])
        if (
            len(external_ids) < required_external
            or len({normalized_groups[key] for key in external_ids}) < required_external
        ):
            raise SelectionLockError(
                "conditional unified evaluation requires two eligible independent external families"
            )
    return promotion


def _joint_bundle_check(
    lock: SelectionLock, task_id: str, bundle_id: str
) -> GateCheck:
    if task_id not in {CELL_TASK, VARIANT_TASK}:
        return GateCheck(
            "joint_outcome_bundle",
            True,
            bundle_id,
            {"shared_tasks": []},
            "task has no co-unblinded GSE289173 peer",
        )
    peer = VARIANT_TASK if task_id == CELL_TASK else CELL_TASK
    peer_decision = _task_decision(lock, peer)
    peer_promotion = _promotion_contract(peer_decision, peer)
    peer_bundle = peer_promotion.get("sealed_outcome_bundle_id")
    return GateCheck(
        "joint_outcome_bundle",
        peer_bundle == bundle_id,
        {task_id: bundle_id, peer: peer_bundle},
        {"same_bundle": True},
        "cell and regulatory GSE289173 outcomes must share one co-unblinded bundle",
    )


def _metric_mapping(metrics: Any) -> Mapping[str, Any]:
    if isinstance(metrics, Mapping):
        return metrics
    nested = _read(metrics, "metrics", _MISSING)
    if isinstance(nested, Mapping):
        return nested
    raise TournamentError("promotion metrics must be a mapping")


def _observed_number(metrics: Mapping[str, Any], name: str) -> float | None:
    if name not in metrics or isinstance(metrics[name], bool):
        return None
    try:
        value = float(metrics[name])
    except (TypeError, ValueError):
        return None
    return value if isfinite(value) else None


def _observed_number_mapping(
    metrics: Mapping[str, Any], name: str
) -> dict[str, float] | None:
    value = metrics.get(name)
    if not isinstance(value, Mapping) or not value:
        return None
    result: dict[str, float] = {}
    for key, item in value.items():
        if not isinstance(key, str) or not key.strip() or isinstance(item, bool):
            return None
        try:
            number = float(item)
        except (TypeError, ValueError):
            return None
        if not isfinite(number):
            return None
        result[key] = number
    return result


def _minimum_check(
    metrics: Mapping[str, Any], name: str, threshold: float, *, strict: bool = False
) -> GateCheck:
    observed = _observed_number(metrics, name)
    passed = observed is not None and (
        observed > threshold if strict else observed >= threshold
    )
    operator = ">" if strict else ">="
    return GateCheck(
        name,
        passed,
        observed,
        {operator: threshold},
        f"{name} must be {operator} {threshold}",
    )


def _maximum_check(metrics: Mapping[str, Any], name: str, threshold: float) -> GateCheck:
    observed = _observed_number(metrics, name)
    return GateCheck(
        name,
        observed is not None and observed <= threshold,
        observed,
        {"<=": threshold},
        f"{name} must be <= {threshold}",
    )


def _boolean_check(
    metrics: Mapping[str, Any], name: str, expected: bool = True
) -> GateCheck:
    observed = metrics.get(name)
    return GateCheck(
        name,
        isinstance(observed, bool) and observed is expected,
        observed,
        {"equals": expected},
        f"{name} must be the boolean {expected!r}",
    )


def _exact_mapping_check(
    metrics: Mapping[str, Any],
    name: str,
    required_ids: Sequence[str],
    minimum: float,
) -> GateCheck:
    observed = _observed_number_mapping(metrics, name)
    required = tuple(required_ids)
    exact = observed is not None and set(observed) == set(required)
    return GateCheck(
        name,
        exact and all(value >= minimum for value in observed.values()),
        observed,
        {"exact_ids": list(required), "all_gte": minimum},
        f"{name} must contain exactly the frozen roster and all values >= {minimum}",
    )


def _receipt_check(name: str, value: Any) -> GateCheck:
    observed = value if isinstance(value, str) else None
    return GateCheck(
        name,
        observed is not None and _SHA256.fullmatch(observed) is not None,
        observed,
        {"sha256": True},
        f"{name} must be a 64-character lowercase SHA-256 receipt",
    )


def _identity_check(name: str, observed: Any, expected: Any) -> GateCheck:
    return GateCheck(
        name,
        observed == expected,
        observed,
        {"equals": expected},
        f"{name} must match its frozen binding",
    )


def _holm_check(metrics: Mapping[str, Any]) -> GateCheck:
    observed = _observed_number(metrics, "holm_adjusted_p")
    threshold = float(PROMOTION_GATE_POLICY["confirmatory_fwer"])
    return GateCheck(
        "holm_adjusted_p",
        observed is not None and 0.0 <= observed <= threshold,
        observed,
        {"between_inclusive": [0.0, threshold]},
        "Holm-adjusted p must be a probability no greater than the frozen FWER",
    )


def _track_checks(
    task_id: str,
    metrics: Mapping[str, Any],
    promotion: Mapping[str, Any],
) -> list[GateCheck]:
    policy = CHAMPION_GATE_THRESHOLDS[task_id]
    if task_id == CELL_TASK:
        return [
            _minimum_check(metrics, policy["primary_metric"], float(policy["primary_min"])),
            _minimum_check(metrics, policy["delta_metric"], float(policy["delta_min"])),
            _minimum_check(
                metrics,
                policy["delta_ci_low_metric"],
                float(policy["delta_ci_low_strict_min"]),
                strict=True,
            ),
            _exact_mapping_check(
                metrics,
                policy["stratum_metric"],
                _unique_strings(promotion["required_class_ids"], "required_class_ids"),
                float(policy["stratum_min"]),
            ),
            _maximum_check(metrics, policy["brier_delta_metric"], float(policy["brier_delta_max"])),
            _holm_check(metrics),
        ]
    if task_id == VARIANT_TASK:
        margin = float(promotion["eqtl_auprc_noninferiority_margin"])
        return [
            _boolean_check(metrics, policy["summary_statistics_flag"]),
            _minimum_check(metrics, policy["delta_metric"], float(policy["delta_min"])),
            _minimum_check(
                metrics,
                policy["delta_ci_low_metric"],
                float(policy["delta_ci_low_strict_min"]),
                strict=True,
            ),
            _exact_mapping_check(
                metrics,
                policy["stratum_metric"],
                _unique_strings(
                    promotion["required_major_lineage_ids"],
                    "required_major_lineage_ids",
                ),
                float(policy["stratum_min"]),
            ),
            _minimum_check(metrics, policy["auprc_delta_ci_low_metric"], -margin),
            _holm_check(metrics),
        ]
    if task_id == RNA_ATAC_TASK:
        required = _unique_strings(
            promotion["required_major_lineage_ids"], "required_major_lineage_ids"
        )
        observed = _observed_number_mapping(metrics, policy["stratum_metric"])
        exact = observed is not None and set(observed) == set(required)
        improved = sum(value > 0.0 for value in observed.values()) if exact else None
        return [
            _minimum_check(metrics, policy["delta_metric"], float(policy["delta_min"])),
            _minimum_check(
                metrics,
                policy["delta_ci_low_metric"],
                float(policy["delta_ci_low_strict_min"]),
                strict=True,
            ),
            _exact_mapping_check(
                metrics, policy["stratum_metric"], required, float(policy["stratum_min"])
            ),
            GateCheck(
                "improved_major_lineages",
                improved is not None and improved >= int(policy["min_improved_strata"]),
                improved,
                {">=": int(policy["min_improved_strata"])},
                "at least four of the exact five major lineages must improve",
            ),
            _holm_check(metrics),
        ]
    if task_id == BULK_TASK:
        required = _unique_strings(promotion["required_endpoint_ids"], "required_endpoint_ids")
        reversals = metrics.get(policy["direction_reversal_metric"])
        exact = isinstance(reversals, Mapping) and set(reversals) == set(required)
        no_reversal = exact and all(value is False for value in reversals.values())
        return [
            _minimum_check(metrics, policy["delta_metric"], float(policy["delta_min"])),
            _minimum_check(
                metrics,
                policy["delta_ci_low_metric"],
                float(policy["delta_ci_low_strict_min"]),
                strict=True,
            ),
            GateCheck(
                policy["direction_reversal_metric"],
                no_reversal,
                reversals,
                {"exact_ids": list(required), "all_false": True},
                "every frozen endpoint must be present and retain direction",
            ),
            _holm_check(metrics),
        ]
    if task_id == GRAPH_TASK:
        return [
            _minimum_check(metrics, policy["delta_metric"], float(policy["delta_min"])),
            _minimum_check(
                metrics,
                policy["delta_ci_low_metric"],
                float(policy["delta_ci_low_strict_min"]),
                strict=True,
            ),
            _maximum_check(metrics, policy["brier_delta_metric"], float(policy["brier_delta_max"])),
            _holm_check(metrics),
        ]
    if task_id == UNIFIED_TASK:
        specialists = _unique_strings(
            promotion["eligible_specialist_ids"], "eligible_specialist_ids"
        )
        margins = promotion["specialist_noninferiority_margins"]
        observed_specialists = _observed_number_mapping(
            metrics, "specialist_delta_ci_low_by_id"
        )
        specialist_pass = (
            observed_specialists is not None
            and set(observed_specialists) == set(specialists)
            and all(observed_specialists[key] >= -float(margins[key]) for key in specialists)
        )
        source_ids = _unique_strings(
            promotion["eligible_source_family_ids"], "eligible_source_family_ids"
        )
        external_ids = _unique_strings(
            promotion["eligible_external_source_family_ids"],
            "eligible_external_source_family_ids",
        )
        groups = promotion["source_independence_group_by_family"]
        p_values = _observed_number_mapping(
            metrics, "holm_adjusted_superiority_p_by_source_family"
        )
        valid_p = (
            p_values is not None
            and set(p_values) == set(source_ids)
            and all(0.0 <= value <= 1.0 for value in p_values.values())
        )
        significant = (
            tuple(
                sorted(
                    key
                    for key, value in p_values.items()
                    if value <= float(PROMOTION_GATE_POLICY["confirmatory_fwer"])
                )
            )
            if valid_p
            else ()
        )
        significant_groups = {groups[key] for key in significant}
        external_significant = set(significant) & set(external_ids)
        external_groups = {groups[key] for key in external_ids}
        return [
            GateCheck(
                "specialist_noninferiority",
                specialist_pass,
                observed_specialists,
                {
                    "exact_ids": list(specialists),
                    "ci_low_gte_negative_margin": dict(margins),
                },
                "every eligible specialist must clear its frozen CI margin",
            ),
            GateCheck(
                "holm_superior_source_families",
                valid_p
                and len(significant) >= int(policy["min_superior_source_families"])
                and len(significant_groups) >= int(policy["min_superior_source_families"]),
                {"adjusted_p": p_values, "significant": significant},
                {"distinct_independence_groups_gte": int(policy["min_superior_source_families"])},
                "Holm superiority must span two non-overlapping source families",
            ),
            GateCheck(
                "external_superiority",
                len(external_significant)
                >= int(policy["min_external_superior_source_families"]),
                sorted(external_significant),
                {">=": int(policy["min_external_superior_source_families"])},
                "at least one Holm-superior source must be genuinely external",
            ),
            GateCheck(
                "eligible_external_source_families",
                len(external_ids) >= int(policy["min_eligible_external_source_families"])
                and len(external_groups)
                >= int(policy["min_eligible_external_source_families"]),
                external_ids,
                {"distinct_groups_gte": int(policy["min_eligible_external_source_families"])},
                "unified promotion requires two eligible independent external families",
            ),
        ]
    raise TournamentError(f"task {task_id!r} has no champion gate")


def enforce_terminal_failure_policy(
    history: Iterable[Any],
    *,
    model_version: str,
    task_id: str,
    sealed_outcome_bundle_id: str,
) -> None:
    model_version = _nonempty_identifier(model_version, "model_version")
    if task_id not in EXPECTED_TASK_IDS | {CONDITIONAL_CONTEXT_TASK}:
        raise TournamentError(f"unsupported terminal-policy task {task_id!r}")
    bundle_id = _nonempty_identifier(
        sealed_outcome_bundle_id, "sealed_outcome_bundle_id"
    )
    for record in history:
        prior_consumed = _read(record, "holdout_consumed", False) is True
        prior_bundle = str(
            _read_first(record, ("sealed_outcome_bundle_id", "holdout_id"), "")
        ).strip()
        if prior_consumed and prior_bundle == bundle_id:
            prior_task = str(_read_first(record, ("task_id", "track"), "")).strip()
            joint_peer = {prior_task, task_id} == {CELL_TASK, VARIANT_TASK}
            if not joint_peer:
                raise TerminalFailureError(
                    f"sealed outcome bundle {bundle_id!r} has already been consumed"
                )
        prior_failed = str(_read(record, "disposition", "")).startswith(
            "terminal_failure"
        )
        prior_version = str(_read(record, "model_version", "")).strip()
        if prior_failed and prior_version == model_version:
            raise TerminalFailureError(
                f"terminally failed model version {model_version!r} cannot be relabeled"
            )


def _power_decision_binding(
    path: str | Path,
) -> tuple[PowerDecision, dict[str, Any], Path, Any]:
    root = _safe_resolve(path, "PowerDecision")
    try:
        decision = verify_power_decision(root, reverify_sources=True)
        evidence = verify_development_power_evidence(
            decision.development_power_evidence_dir,
            reverify_sources=True,
        )
    except (FirewallError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid production PowerDecision: {exc}") from exc
    binding, _ = _frozen_document_binding(
        root, filename="power_decision.json", artifact_class="power_decision"
    )
    return (
        decision,
        {
            **binding,
            "power_decision_sha256": decision.power_decision_sha256,
            "power_result_sha256": decision.power_result_sha256,
            "development_power_evidence_dir": decision.development_power_evidence_dir,
            "development_power_evidence_sha256": (
                decision.development_power_evidence_sha256
            ),
            "development_power_evidence_manifest_sha256": (
                decision.development_power_evidence_manifest_sha256
            ),
        },
        root,
        evidence,
    )


def _sealed_outcome_binding(
    path: str | Path,
) -> tuple[Any, dict[str, Any], Path]:
    root = _safe_resolve(path, "sealed outcome bundle")
    try:
        bundle = verify_sealed_outcome_bundle(root)
    except (FirewallError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid sealed outcome bundle: {exc}") from exc
    binding, _ = _frozen_document_binding(
        root,
        filename="sealed_outcome_bundle.json",
        artifact_class="sealed_outcome_bundle",
    )
    return bundle, {
        **binding,
        "outcome_bundle_id": bundle.bundle_id,
        "outcome_bundle_sha256": bundle.bundle_sha256,
    }, root


def _outcome_consumption_binding(
    path: str | Path,
) -> tuple[OutcomeConsumption, dict[str, Any], Path]:
    marker = _safe_resolve(path, "OutcomeConsumption marker")
    try:
        consumption = verify_outcome_consumption(marker)
    except (FirewallError, OSError, ValueError, KeyError, TypeError) as exc:
        raise TournamentError(f"invalid OutcomeConsumption marker: {exc}") from exc
    return consumption, {
        "path": marker.as_posix(),
        "document_sha256": sha256_file(marker),
        "consumption_sha256": consumption.consumption_sha256,
    }, marker


def _commit_coverage_mode(
    *, lock: SelectionLock, task_id: str, commits: Sequence[PredictionCommit]
) -> str:
    selected, baseline = _required_commit_runs(lock=lock, task_id=task_id)
    observed = tuple(sorted(commit.run_id for commit in commits))
    if not observed:
        return "none"
    if observed == selected:
        _assert_prediction_commit_coverage(
            lock=lock, task_id=task_id, commits=commits, required="selected"
        )
        return "selected"
    if observed == tuple(sorted((*selected, *baseline))):
        _assert_prediction_commit_coverage(
            lock=lock,
            task_id=task_id,
            commits=commits,
            required="selected_and_baseline",
        )
        return "selected_and_baseline"
    secondary = tuple(
        run_id for run_id, _model_id in _secondary_commit_runs(lock=lock, task_id=task_id)
    )
    if secondary and observed == tuple(sorted((*selected, *baseline, *secondary))):
        _assert_prediction_commit_coverage(
            lock=lock,
            task_id=task_id,
            commits=commits,
            required="selected_baseline_secondary",
        )
        return "selected_baseline_secondary"
    raise TournamentError(
        f"{task_id} PredictionCommits do not match an allowed locked coverage set"
    )


def _join_marker_identity(
    commit: PredictionCommit,
    join_ids: tuple[
        tuple[str, ...], tuple[str, ...], tuple[tuple[str, str], ...]
    ],
) -> dict[str, Any]:
    row_ids, unit_ids, row_unit_pairs = join_ids
    return {
        "n_join_units": len(unit_ids),
        "n_join_rows": len(row_ids),
        "unit_set_sha256": canonical_sha256(
            {"unit_id_field": commit.unit_id_field, "unit_ids": list(unit_ids)}
        ),
        "row_set_sha256": canonical_sha256(
            {
                "row_id_field": commit.row_id_field,
                "unit_id_field": commit.unit_id_field,
                "row_unit_pairs": [
                    {"row_id": row_id, "unit_id": unit_id}
                    for row_id, unit_id in row_unit_pairs
                ],
            }
        ),
    }


def _validate_opened_sources(
    *,
    lock: SelectionLock,
    lock_binding: Mapping[str, Any],
    task_id: str,
    promotion: Mapping[str, Any],
    commits: Sequence[PredictionCommit],
    power: PowerDecision,
    outcome: Any,
    outcome_binding: Mapping[str, Any],
    outcome_root: Path,
    consumption: OutcomeConsumption,
) -> None:
    expected_datasets = tuple(_task_decision(lock, task_id)["sealed_dataset_ids"])
    task_bindings = [
        item for item in outcome.task_bindings if item.get("task_id") == task_id
    ]
    if len(task_bindings) != 1 or tuple(task_bindings[0]["dataset_ids"]) != expected_datasets:
        raise TournamentError("sealed outcome task binding differs from SelectionLock")
    if (
        consumption.selection_lock_id != lock.lock_id
        or consumption.selection_lock_manifest_sha256
        != lock_binding["manifest_sha256"]
        or consumption.outcome_bundle_id != outcome.bundle_id
        or consumption.outcome_bundle_sha256 != outcome.bundle_sha256
        or consumption.outcome_manifest_sha256
        != outcome_binding["manifest_sha256"]
    ):
        raise TournamentError("OutcomeConsumption does not bind the verified lock/outcome")
    power_bindings = [
        item for item in consumption.power_bindings if item.get("task_id") == task_id
    ]
    if len(power_bindings) != 1:
        raise TournamentError("OutcomeConsumption lacks one task PowerDecision binding")
    join_ids = _verify_exact_task_units(
        commits=commits, outcome_root=outcome_root, outcome_bundle=outcome
    )
    expected_power_binding = {
        "task_id": task_id,
        "power_decision_sha256": power.power_decision_sha256,
        "power_result_sha256": power.power_result_sha256,
        "prediction_set_sha256": power.prediction_set_sha256,
        "prediction_bindings": [dict(item) for item in power.prediction_bindings],
        **_join_marker_identity(commits[0], join_ids),
    }
    if dict(power_bindings[0]) != expected_power_binding:
        raise TournamentError(
            "OutcomeConsumption power and exact row/unit bindings changed"
        )


def _prediction_rows_and_bundle(
    commit: PredictionCommit,
) -> tuple[PredictionBundle, Path, list[dict[str, str]]]:
    bundle_path = _safe_resolve(
        commit.prediction_bundle_path, "PredictionCommit PredictionBundle"
    )
    try:
        bundle = PredictionBundle.load_json(bundle_path)
        bundle.validate_artifacts(bundle_path.parent)
    except (ContractError, OSError, ValueError) as exc:
        raise TournamentError(
            f"invalid PredictionBundle for sealed metrics: {exc}"
        ) from exc
    fields = _PREDICTION_TABLE_FIELDS.get(commit.task_id)
    if fields is None:
        raise TournamentError(f"task {commit.task_id} has no prediction schema")
    rows = _read_exact_tsv(
        bundle.standardized_table.validate(bundle_path.parent, require_relative=True),
        fields,
        f"sealed metric PredictionBundle {bundle.bundle_id}",
    )
    if len(rows) != commit.n_predictions:
        raise TournamentError("sealed metric PredictionBundle count changed")
    return bundle, bundle_path, rows


def _one_bundle_artifact(
    bundle: PredictionBundle,
    *,
    root: Path,
    role: str,
) -> Path:
    matches = [artifact for artifact in bundle.artifacts if artifact.role == role]
    if len(matches) != 1 or matches[0].media_type != "text/tab-separated-values":
        raise TournamentError(f"PredictionBundle requires exactly one TSV {role!r}")
    return matches[0].validate(root, require_relative=True)


def _cell_probability_rows(
    *,
    bundle: PredictionBundle,
    bundle_path: Path,
    prediction_rows: Sequence[Mapping[str, str]],
    class_roster: Sequence[str],
) -> list[dict[str, Any]]:
    role = f"class_probabilities:{CELL_TASK}"
    probability_fields = tuple(f"probability::{class_id}" for class_id in class_roster)
    fields = ("row_hash", "unit_hash", *probability_fields)
    raw_rows = _read_exact_tsv(
        _one_bundle_artifact(bundle, root=bundle_path.parent, role=role),
        fields,
        f"{bundle.bundle_id} class probabilities",
    )
    prediction_identity = [
        (row["row_hash"], row["unit_hash"]) for row in prediction_rows
    ]
    probability_identity = [
        (row["row_hash"], row["unit_hash"]) for row in raw_rows
    ]
    if probability_identity != prediction_identity:
        raise TournamentError("cell probability rows differ from hard predictions")
    return [
        {
            "row_hash": row["row_hash"],
            "unit_hash": row["unit_hash"],
            "probabilities": {
                class_id: row[f"probability::{class_id}"]
                for class_id in class_roster
            },
        }
        for row in raw_rows
    ]


def _outcome_rows_for_role(
    *,
    outcome: Any,
    outcome_root: Path,
    role: str,
    fields: Sequence[str],
) -> list[dict[str, str]]:
    matches = [artifact for artifact in outcome.artifacts if artifact.role == role]
    if len(matches) != 1 or matches[0].media_type != "text/tab-separated-values":
        raise TournamentError(f"sealed outcome requires exactly one TSV {role!r}")
    return _read_exact_tsv(
        matches[0].validate(outcome_root, require_relative=True),
        fields,
        f"sealed outcome {role}",
    )


def _sealed_metric_endpoint_parameters(
    task_id: str, promotion: Mapping[str, Any]
) -> dict[str, Any]:
    if task_id == CELL_TASK:
        return {
            "class_roster": list(
                _unique_strings(
                    promotion["required_class_ids"],
                    "required_class_ids",
                    sorted_required=True,
                )
            )
        }
    if task_id == VARIANT_TASK:
        return {
            "strata": list(
                _unique_strings(
                    promotion["required_major_lineage_ids"],
                    "required_major_lineage_ids",
                    sorted_required=True,
                )
            )
        }
    raise TournamentError(
        f"task {task_id!r} has no positive sealed-metric implementation"
    )


def _sealed_metric_sources(
    *,
    selection_lock_dir: str | Path,
    task_id: str,
    prediction_commit_dirs: Iterable[str | Path],
    power_decision_dir: str | Path,
    sealed_outcome_bundle_dir: str | Path,
    outcome_consumption_path: str | Path,
) -> dict[str, Any]:
    lock, lock_binding, _ = _selection_lock_binding(selection_lock_dir)
    decision = _task_decision(lock, task_id)
    promotion = _promotion_contract(decision, task_id)
    if promotion["promotion_mode"] not in PROMOTABLE_MODES:
        raise TournamentError("sealed metrics require a promotable locked task")
    variant_primary_capability = decision.get("variant_primary_capability")
    variant_secondary_evaluation = decision.get("variant_secondary_evaluation")
    if task_id == VARIANT_TASK:
        if not isinstance(variant_primary_capability, Mapping):
            raise TournamentError(
                "variant sealed metrics require a locked primary-capability binding"
            )
        # dict() is SHALLOW: it converts the outer mappingproxy but leaves every
        # nested tuple and mappingproxy from contracts._as_metadata in place.
        # These values are written to JSON at freeze time (tuple -> list) and
        # compared against this in-memory rederivation at verify time, so a
        # shallow copy can never round-trip.  canonicalize() converts deeply to
        # plain JSON types, which is what the frozen document actually holds.
        variant_primary_capability = _canonical_locked_binding(
            variant_primary_capability, "variant primary capability"
        )
        if not isinstance(variant_secondary_evaluation, Mapping):
            raise TournamentError(
                "variant sealed metrics require a locked secondary-evaluation binding"
            )
        variant_secondary_evaluation = _canonical_locked_binding(
            variant_secondary_evaluation, "variant secondary evaluation"
        )
    elif (
        variant_primary_capability is not None
        or variant_secondary_evaluation is not None
    ):
        raise TournamentError(
            "nonvariant sealed metrics cannot carry variant evaluation bindings"
        )
    endpoint_parameters = _sealed_metric_endpoint_parameters(task_id, promotion)
    bindings, commits, roots = _verified_task_prediction_commits(
        prediction_commit_dirs,
        lock=lock,
        task_id=task_id,
        allow_empty=False,
    )
    required_coverage = (
        "selected_baseline_secondary"
        if task_id == VARIANT_TASK
        else "selected_and_baseline"
    )
    if _commit_coverage_mode(
        lock=lock, task_id=task_id, commits=commits
    ) != required_coverage:
        raise TournamentError(
            "sealed metrics require the complete locked prediction set"
        )
    power, power_binding, power_root, evidence = _power_decision_binding(
        power_decision_dir
    )
    outcome, outcome_binding, outcome_root = _sealed_outcome_binding(
        sealed_outcome_bundle_dir
    )
    consumption, consumption_binding, consumption_marker = (
        _outcome_consumption_binding(outcome_consumption_path)
    )
    if (
        power.selection_lock_id != lock.lock_id
        or power.task_id != task_id
        or power.power_result.get("passed") is not True
    ):
        raise TournamentError("sealed metric PowerDecision does not pass the locked task")
    _validate_opened_sources(
        lock=lock,
        lock_binding=lock_binding,
        task_id=task_id,
        promotion=promotion,
        commits=commits,
        power=power,
        outcome=outcome,
        outcome_binding=outcome_binding,
        outcome_root=outcome_root,
        consumption=consumption,
    )

    seed_by_run: dict[str, tuple[str, int]] = {}
    for role in ("selected", "baseline"):
        for run in decision[f"{role}_runs"]:
            seed_by_run[str(run["run_id"])] = (role, int(run["seed"]))
    selected_rows: dict[int, Sequence[Mapping[str, Any]]] = {}
    baseline_rows: dict[int, Sequence[Mapping[str, Any]]] = {}
    selected_probabilities: dict[int, Sequence[Mapping[str, Any]]] = {}
    baseline_probabilities: dict[int, Sequence[Mapping[str, Any]]] = {}
    secondary_rows_by_model: dict[str, Sequence[Mapping[str, Any]]] = {}
    comparator_by_run = {
        str(comparator["run_id"]): str(comparator["model_id"])
        for comparator in (
            variant_secondary_evaluation["comparators"]
            if task_id == VARIANT_TASK
            else ()
        )
    }
    for commit in commits:
        bundle, bundle_path, rows = _prediction_rows_and_bundle(commit)
        if commit.run_id in comparator_by_run:
            model_id = comparator_by_run[commit.run_id]
            if model_id in secondary_rows_by_model:
                raise TournamentError(
                    f"sealed metrics repeat secondary comparator {model_id}"
                )
            secondary_rows_by_model[model_id] = rows
            continue
        role, seed = seed_by_run[commit.run_id]
        target = selected_rows if role == "selected" else baseline_rows
        if seed in target:
            raise TournamentError(f"sealed metrics repeat {role} seed {seed}")
        target[seed] = rows
        if task_id == CELL_TASK:
            probability_rows = _cell_probability_rows(
                bundle=bundle,
                bundle_path=bundle_path,
                prediction_rows=rows,
                class_roster=endpoint_parameters["class_roster"],
            )
            probability_target = (
                selected_probabilities if role == "selected" else baseline_probabilities
            )
            probability_target[seed] = probability_rows

    outcome_rows = _outcome_rows_for_role(
        outcome=outcome,
        outcome_root=outcome_root,
        role=f"outcome:{task_id}",
        fields=_OUTCOME_TABLE_FIELDS[task_id],
    )
    secondary_rows = None
    if task_id == VARIANT_TASK:
        secondary_rows = _outcome_rows_for_role(
            outcome=outcome,
            outcome_root=outcome_root,
            role=f"secondary_outcome:{VARIANT_TASK}",
            fields=(
                "row_hash",
                "unit_hash",
                "block_hash",
                "stratum",
                "observed_binary",
                "summary_statistics_complete",
            ),
        )
    evaluator_source_bundle = _sealed_evaluator_source_bundle()
    cache_key = canonical_sha256(
        {
            "evaluator_source_bundle": evaluator_source_bundle,
            "task_id": task_id,
            "selected_rows_by_seed": [
                {"seed": seed, "rows": list(selected_rows[seed])}
                for seed in sorted(selected_rows)
            ],
            "baseline_rows_by_seed": [
                {"seed": seed, "rows": list(baseline_rows[seed])}
                for seed in sorted(baseline_rows)
            ],
            "selected_probabilities_by_seed": [
                {"seed": seed, "rows": list(selected_probabilities[seed])}
                for seed in sorted(selected_probabilities)
            ],
            "baseline_probabilities_by_seed": [
                {"seed": seed, "rows": list(baseline_probabilities[seed])}
                for seed in sorted(baseline_probabilities)
            ],
            "outcome_rows": outcome_rows,
            "secondary_outcome_rows": secondary_rows,
            "endpoint_parameters": endpoint_parameters,
            "variant_primary_capability": variant_primary_capability,
            "variant_secondary_evaluation": variant_secondary_evaluation,
            "variant_secondary_rows_by_model": {
                model_id: list(secondary_rows_by_model[model_id])
                for model_id in sorted(secondary_rows_by_model)
            },
            "n_resamples": _POWER_BOOTSTRAP_RESAMPLES,
            "bootstrap_seed": _POWER_BOOTSTRAP_SEED,
        }
    )
    result = _SEALED_METRIC_RESULT_CACHE.get(cache_key)
    if result is None:
        try:
            result = recompute_sealed_metrics(
                task_id=task_id,
                selected_rows_by_seed=selected_rows,
                baseline_rows_by_seed=baseline_rows,
                selected_probabilities_by_seed=(
                    selected_probabilities if task_id == CELL_TASK else None
                ),
                baseline_probabilities_by_seed=(
                    baseline_probabilities if task_id == CELL_TASK else None
                ),
                outcome_rows=outcome_rows,
                secondary_outcome_rows=secondary_rows,
                endpoint_parameters=endpoint_parameters,
                n_resamples=_POWER_BOOTSTRAP_RESAMPLES,
                bootstrap_seed=_POWER_BOOTSTRAP_SEED,
                variant_primary_capability=variant_primary_capability,
                variant_secondary_evaluation=variant_secondary_evaluation,
                variant_secondary_rows_by_model=(
                    secondary_rows_by_model if task_id == VARIANT_TASK else None
                ),
            )
        except (SealedMetricError, ValueError, KeyError, TypeError) as exc:
            raise TournamentError(f"sealed metric evaluation failed: {exc}") from exc
        _SEALED_METRIC_RESULT_CACHE[cache_key] = result
    return {
        "lock": lock,
        "lock_binding": lock_binding,
        "decision": decision,
        "promotion": promotion,
        "prediction_commit_bindings": list(bindings),
        "prediction_commit_roots": roots,
        "power_binding": power_binding,
        "power_root": power_root,
        "power_evidence": evidence,
        "outcome_binding": outcome_binding,
        "outcome_root": outcome_root,
        "consumption_binding": consumption_binding,
        "consumption_marker": consumption_marker,
        "endpoint_parameters": endpoint_parameters,
        "variant_primary_capability": variant_primary_capability,
        "variant_secondary_evaluation": variant_secondary_evaluation,
        "evaluator_source_bundle": evaluator_source_bundle,
        "result": result,
    }


def _canonical_locked_binding(value: Any, label: str) -> dict[str, Any]:
    """Deep-convert a fixed binding to the plain JSON types a frozen doc holds."""

    try:
        converted = canonicalize(value)
    except (HashingError, TypeError, ValueError) as exc:
        raise TournamentError(f"{label} is not canonical JSON: {exc}") from exc
    if not isinstance(converted, dict):
        raise TournamentError(f"{label} must be an object")
    return converted


def _sealed_metric_identity(
    sources: Mapping[str, Any], joined_table: ArtifactRef
) -> dict[str, Any]:
    result = sources["result"]
    return {
        "schema_version": SEALED_METRIC_BUNDLE_SCHEMA_VERSION,
        "task_id": result.task_id,
        "selection_lock_binding": dict(sources["lock_binding"]),
        "selection_lock_id": sources["lock"].lock_id,
        "prediction_commit_bindings": list(sources["prediction_commit_bindings"]),
        "power_decision_binding": dict(sources["power_binding"]),
        "sealed_outcome_binding": dict(sources["outcome_binding"]),
        "outcome_consumption_binding": dict(sources["consumption_binding"]),
        "evaluator_id": result.evaluator_id,
        "evaluator_source_sha256": sources["evaluator_source_bundle"][
            "bundle_sha256"
        ],
        "evaluator_source_bundle": dict(sources["evaluator_source_bundle"]),
        "endpoint_parameters": dict(sources["endpoint_parameters"]),
        "variant_primary_capability": sources["variant_primary_capability"],
        "variant_secondary_evaluation": sources[
            "variant_secondary_evaluation"
        ],
        "n_resamples": result.n_resamples,
        "bootstrap_seed": result.bootstrap_seed,
        "independent_unit": result.independent_unit,
        "independent_unit_count": result.independent_unit_count,
        "row_count": result.row_count,
        "metrics": dict(result.metrics),
        "joined_ensemble_table": joined_table.to_dict(),
        "seeds_define_one_ensemble": True,
        "seeds_are_biological_replicates": False,
    }


def freeze_sealed_metric_bundle(
    *,
    selection_lock_dir: str | Path,
    task_id: str,
    prediction_commit_dirs: Iterable[str | Path],
    power_decision_dir: str | Path,
    sealed_outcome_bundle_dir: str | Path,
    outcome_consumption_path: str | Path,
    output_root: str | Path,
) -> Path:
    """Freeze task-native metrics from the already consumed held-back outcome."""

    sources = _sealed_metric_sources(
        selection_lock_dir=selection_lock_dir,
        task_id=task_id,
        prediction_commit_dirs=prediction_commit_dirs,
        power_decision_dir=power_decision_dir,
        sealed_outcome_bundle_dir=sealed_outcome_bundle_dir,
        outcome_consumption_path=outcome_consumption_path,
    )
    output = _safe_resolve(output_root, "sealed metric bundle output root")
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".sealed-metric.", dir=output))
    try:
        fields = _SEALED_ENSEMBLE_FIELDS[task_id]
        table_path = staging / "joined_ensemble.tsv"
        write_text_exclusive(
            table_path,
            _tsv_text(fields, sources["result"].ensemble_rows),
            mode=0o440,
        )
        joined_ref = ArtifactRef.from_path(
            table_path,
            relative_to=staging,
            media_type="text/tab-separated-values",
            role=f"joined_sealed_ensemble:{task_id}",
        )
        identity = _sealed_metric_identity(sources, joined_ref)
        payload = {"sealed_metric_bundle_id": canonical_sha256(identity), **identity}
        write_json_exclusive(
            staging / "sealed_metric_bundle.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "sealed_metric_bundle",
                "sealed_metric_bundle_id": payload["sealed_metric_bundle_id"],
                "task_id": task_id,
            },
        )
        target = output / (
            f"sealed-metric--{task_id}--{payload['sealed_metric_bundle_id']}"
        )
        if target.exists():
            raise TournamentError(f"sealed metric bundle already exists: {target}")
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_sealed_metric_bundle(target, reverify_sources=True)
    return target


def verify_sealed_metric_bundle(
    path: str | Path, *, reverify_sources: bool = True
) -> dict[str, Any]:
    """Verify and optionally rederive one task-native held-back metric bundle."""

    root = _safe_resolve(path, "sealed metric bundle")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "sealed_metric_bundle.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid sealed metric bundle: {exc}") from exc
    if (
        not isinstance(payload, Mapping)
        or payload.get("schema_version") != SEALED_METRIC_BUNDLE_SCHEMA_VERSION
    ):
        raise TournamentError("unsupported sealed metric bundle schema")
    claimed = _sha256_identifier(
        payload.get("sealed_metric_bundle_id"), "sealed_metric_bundle_id"
    )
    task_id = _nonempty_identifier(payload.get("task_id"), "task_id")
    if task_id not in _SEALED_ENSEMBLE_FIELDS:
        raise TournamentError("sealed metric task has no positive evaluator")
    try:
        joined_ref = ArtifactRef.from_dict(payload.get("joined_ensemble_table"))
    except ContractError as exc:
        raise TournamentError(f"invalid sealed ensemble ArtifactRef: {exc}") from exc
    if (
        joined_ref.role != f"joined_sealed_ensemble:{task_id}"
        or joined_ref.media_type != "text/tab-separated-values"
    ):
        raise TournamentError("sealed ensemble table has the wrong role or media type")
    joined_rows = _read_exact_tsv(
        joined_ref.validate(root, require_relative=True),
        _SEALED_ENSEMBLE_FIELDS[task_id],
        "joined sealed ensemble",
    )
    identity = dict(payload)
    identity.pop("sealed_metric_bundle_id", None)
    if canonical_sha256(identity) != claimed:
        raise TournamentError("sealed metric bundle identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "sealed_metric_bundle",
        "sealed_metric_bundle_id": claimed,
        "task_id": task_id,
    }:
        raise TournamentError("sealed metric manifest metadata mismatch")
    if reverify_sources:
        lock_binding = payload.get("selection_lock_binding")
        prediction_bindings = payload.get("prediction_commit_bindings")
        power_binding = payload.get("power_decision_binding")
        outcome_binding = payload.get("sealed_outcome_binding")
        consumption_binding = payload.get("outcome_consumption_binding")
        if not isinstance(lock_binding, Mapping):
            raise TournamentError("sealed metric lacks a SelectionLock binding")
        if not isinstance(prediction_bindings, list) or not all(
            isinstance(item, Mapping) for item in prediction_bindings
        ):
            raise TournamentError("sealed metric prediction bindings are invalid")
        if not all(
            isinstance(item, Mapping)
            for item in (power_binding, outcome_binding, consumption_binding)
        ):
            raise TournamentError("sealed metric opened-source bindings are invalid")
        sources = _sealed_metric_sources(
            selection_lock_dir=lock_binding["path"],
            task_id=task_id,
            prediction_commit_dirs=[item["path"] for item in prediction_bindings],
            power_decision_dir=power_binding["path"],
            sealed_outcome_bundle_dir=outcome_binding["path"],
            outcome_consumption_path=consumption_binding["path"],
        )
        if joined_rows != [dict(row) for row in sources["result"].ensemble_rows]:
            raise TournamentError("joined sealed ensemble does not rederive")
        expected = _canonical_locked_binding(
            {
                "sealed_metric_bundle_id": claimed,
                **_sealed_metric_identity(sources, joined_ref),
            },
            "rederived sealed metric identity",
        )
        if _canonical_locked_binding(payload, "sealed metric bundle") != expected:
            raise TournamentError("sealed metric bundle does not rederive from sources")
    return dict(payload)


def _sealed_metric_bundle_binding(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, Any], Path]:
    root = _safe_resolve(path, "sealed metric bundle")
    payload = verify_sealed_metric_bundle(root, reverify_sources=True)
    binding, _ = _frozen_document_binding(
        root,
        filename="sealed_metric_bundle.json",
        artifact_class="sealed_metric_bundle",
    )
    return payload, {
        **binding,
        "sealed_metric_bundle_id": payload["sealed_metric_bundle_id"],
        "task_id": payload["task_id"],
    }, root


def _confirmatory_multiplicity_sources(
    *,
    selection_lock_dir: str | Path,
    sealed_metric_bundle_dirs: Iterable[str | Path],
) -> dict[str, Any]:
    lock, lock_binding, _ = _selection_lock_binding(selection_lock_dir)
    expected_tasks = sorted(
        str(decision["task_id"])
        for decision in lock.task_decisions
        if _promotion_contract(decision, str(decision["task_id"]))[
            "promotion_mode"
        ]
        in PROMOTABLE_MODES
    )
    paths = tuple(
        _safe_resolve(path, "sealed metric bundle")
        for path in sealed_metric_bundle_dirs
    )
    if len(paths) != len(set(paths)):
        raise TournamentError("confirmatory multiplicity repeats a metric bundle")
    records = [_sealed_metric_bundle_binding(path) for path in paths]
    records.sort(key=lambda item: str(item[0]["task_id"]))
    observed_tasks = [str(item[0]["task_id"]) for item in records]
    if observed_tasks != expected_tasks:
        raise TournamentError(
            "confirmatory multiplicity must cover every promotable SelectionLock task"
        )
    raw_p: dict[str, float] = {}
    for payload, _, _ in records:
        if payload["selection_lock_binding"] != lock_binding:
            raise TournamentError("sealed metric bundles use different SelectionLocks")
        value = _finite_number(
            payload.get("metrics", {}).get("raw_confirmatory_p"),
            f"{payload['task_id']}.raw_confirmatory_p",
        )
        if not 0.0 <= value <= 1.0:
            raise TournamentError("raw confirmatory p-values must lie in [0, 1]")
        raw_p[str(payload["task_id"])] = value
    adjusted = holm_correction(raw_p)
    if not isinstance(adjusted, dict):
        raise TournamentError("Holm correction returned the wrong result shape")
    return {
        "lock": lock,
        "lock_binding": lock_binding,
        "task_ids": expected_tasks,
        "metric_bindings": [item[1] for item in records],
        "metric_roots": [item[2] for item in records],
        "raw_p_by_task": raw_p,
        "adjusted_p_by_task": {
            task_id: float(adjusted[task_id]) for task_id in expected_tasks
        },
    }


def _confirmatory_multiplicity_identity(
    sources: Mapping[str, Any]
) -> dict[str, Any]:
    return {
        "schema_version": CONFIRMATORY_MULTIPLICITY_BUNDLE_SCHEMA_VERSION,
        "selection_lock_binding": dict(sources["lock_binding"]),
        "selection_lock_id": sources["lock"].lock_id,
        "family_id": "externally_eligible_primary_tracks",
        "method": "Holm",
        "fwer": float(PROMOTION_GATE_POLICY["confirmatory_fwer"]),
        "task_ids": list(sources["task_ids"]),
        "sealed_metric_bindings": list(sources["metric_bindings"]),
        "raw_p_by_task": dict(sources["raw_p_by_task"]),
        "adjusted_p_by_task": dict(sources["adjusted_p_by_task"]),
    }


def freeze_confirmatory_multiplicity_bundle(
    *,
    selection_lock_dir: str | Path,
    sealed_metric_bundle_dirs: Iterable[str | Path],
    output_root: str | Path,
) -> Path:
    """Freeze the exact Holm family across all adoptable fixed tracks."""

    sources = _confirmatory_multiplicity_sources(
        selection_lock_dir=selection_lock_dir,
        sealed_metric_bundle_dirs=sealed_metric_bundle_dirs,
    )
    identity = _confirmatory_multiplicity_identity(sources)
    payload = {"multiplicity_bundle_id": canonical_sha256(identity), **identity}
    output = _safe_resolve(output_root, "confirmatory multiplicity output root")
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"confirmatory-multiplicity--{payload['multiplicity_bundle_id']}"
    if target.exists():
        raise TournamentError(f"confirmatory multiplicity bundle already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".confirmatory-multiplicity.", dir=output))
    try:
        write_json_exclusive(
            staging / "confirmatory_multiplicity_bundle.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "confirmatory_multiplicity_bundle",
                "multiplicity_bundle_id": payload["multiplicity_bundle_id"],
            },
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_confirmatory_multiplicity_bundle(target, reverify_sources=True)
    return target


def verify_confirmatory_multiplicity_bundle(
    path: str | Path, *, reverify_sources: bool = True
) -> dict[str, Any]:
    root = _safe_resolve(path, "confirmatory multiplicity bundle")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "confirmatory_multiplicity_bundle.json").read_text(
                encoding="utf-8"
            )
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid confirmatory multiplicity bundle: {exc}") from exc
    if (
        not isinstance(payload, Mapping)
        or payload.get("schema_version")
        != CONFIRMATORY_MULTIPLICITY_BUNDLE_SCHEMA_VERSION
    ):
        raise TournamentError("unsupported confirmatory multiplicity schema")
    claimed = _sha256_identifier(
        payload.get("multiplicity_bundle_id"), "multiplicity_bundle_id"
    )
    identity = dict(payload)
    identity.pop("multiplicity_bundle_id", None)
    if canonical_sha256(identity) != claimed:
        raise TournamentError("confirmatory multiplicity identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "confirmatory_multiplicity_bundle",
        "multiplicity_bundle_id": claimed,
    }:
        raise TournamentError("confirmatory multiplicity manifest metadata mismatch")
    if reverify_sources:
        lock_binding = payload.get("selection_lock_binding")
        metric_bindings = payload.get("sealed_metric_bindings")
        if not isinstance(lock_binding, Mapping) or not isinstance(
            metric_bindings, list
        ) or not all(isinstance(item, Mapping) for item in metric_bindings):
            raise TournamentError("confirmatory multiplicity bindings are invalid")
        sources = _confirmatory_multiplicity_sources(
            selection_lock_dir=lock_binding["path"],
            sealed_metric_bundle_dirs=[item["path"] for item in metric_bindings],
        )
        expected = {
            "multiplicity_bundle_id": claimed,
            **_confirmatory_multiplicity_identity(sources),
        }
        if dict(payload) != expected:
            raise TournamentError(
                "confirmatory multiplicity does not rederive from sealed metrics"
            )
    return dict(payload)


def _confirmatory_multiplicity_binding(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, Any], Path]:
    root = _safe_resolve(path, "confirmatory multiplicity bundle")
    payload = verify_confirmatory_multiplicity_bundle(root, reverify_sources=True)
    binding, _ = _frozen_document_binding(
        root,
        filename="confirmatory_multiplicity_bundle.json",
        artifact_class="confirmatory_multiplicity_bundle",
    )
    return payload, {
        **binding,
        "multiplicity_bundle_id": payload["multiplicity_bundle_id"],
    }, root
    if outcome.bundle_id != promotion["sealed_outcome_bundle_id"]:
        raise TournamentError("sealed outcome bundle ID differs from SelectionLock gate")


def _champion_gate_sources(
    *,
    selection_lock_dir: str | Path,
    task_id: str,
    prediction_commit_dirs: Iterable[str | Path],
    power_decision_dir: str | Path | None,
    sealed_outcome_bundle_dir: str | Path | None,
    outcome_consumption_path: str | Path | None,
    sealed_metric_bundle_dir: str | Path | None,
    confirmatory_multiplicity_bundle_dir: str | Path | None,
) -> dict[str, Any]:
    lock, lock_binding, _ = _selection_lock_binding(selection_lock_dir)
    decision = _task_decision(lock, task_id)
    promotion = _promotion_contract(decision, task_id)
    promotion_mode = str(promotion["promotion_mode"])
    commit_paths = tuple(
        _safe_resolve(path, "PredictionCommit") for path in prediction_commit_dirs
    )
    commit_bindings, commits, commit_roots = _verified_task_prediction_commits(
        commit_paths,
        lock=lock,
        task_id=task_id,
        allow_empty=True,
    )
    coverage = _commit_coverage_mode(lock=lock, task_id=task_id, commits=commits)
    power = None
    power_binding = None
    power_root = None
    evidence = None
    if power_decision_dir is not None:
        power, power_binding, power_root, evidence = _power_decision_binding(
            power_decision_dir
        )
        if (
            power.selection_lock_id != lock.lock_id
            or power.selection_lock_manifest_sha256
            != lock_binding["manifest_sha256"]
            or power.task_id != task_id
            or power.power_method_id
            != CHAMPION_GATE_THRESHOLDS[task_id].get("power_method_id")
            or evidence.promotion_gate_id != decision["promotion_gate"]
        ):
            raise TournamentError("PowerDecision does not bind the locked task/method/gate")
        # firewall._complete_prediction_set REQUIRES the variant task's power
        # decision to bind every secondary comparator run as well, so the two
        # modules must agree on which coverage sets are eligible here.  The
        # check needs the selected and baseline commits to be complete; the
        # auxiliary comparators are additionally present for the variant task
        # and are non-scoring, so their presence must not be rejected.
        if coverage not in {"selected_and_baseline", "selected_baseline_secondary"}:
            raise TournamentError(
                "a PowerDecision requires explicit complete selected+baseline commits"
            )
        observed_power_commits = {
            commit.prediction_commit_sha256 for commit in commits
        }
        if observed_power_commits != {
            item["prediction_commit_sha256"] for item in power.prediction_bindings
        }:
            raise TournamentError("PowerDecision prediction set differs from gate commits")

    outcome = None
    outcome_binding = None
    outcome_root = None
    consumption = None
    consumption_binding = None
    consumption_marker = None
    if (sealed_outcome_bundle_dir is None) != (outcome_consumption_path is None):
        raise TournamentError(
            "sealed outcome and OutcomeConsumption marker must be supplied together"
        )
    if sealed_outcome_bundle_dir is not None:
        outcome, outcome_binding, outcome_root = _sealed_outcome_binding(
            sealed_outcome_bundle_dir
        )
        consumption, consumption_binding, consumption_marker = (
            _outcome_consumption_binding(outcome_consumption_path)
        )

    metric_bundle = None
    metric_binding = None
    metric_root = None
    multiplicity_bundle = None
    multiplicity_binding = None
    multiplicity_root = None
    if (sealed_metric_bundle_dir is None) != (
        confirmatory_multiplicity_bundle_dir is None
    ):
        raise TournamentError(
            "sealed metric and confirmatory multiplicity bundles must be supplied together"
        )
    if sealed_metric_bundle_dir is not None:
        if outcome is None:
            raise TournamentError("sealed metrics cannot precede outcome consumption")
        metric_bundle, metric_binding, metric_root = _sealed_metric_bundle_binding(
            sealed_metric_bundle_dir
        )
        (
            multiplicity_bundle,
            multiplicity_binding,
            multiplicity_root,
        ) = _confirmatory_multiplicity_binding(
            confirmatory_multiplicity_bundle_dir
        )

    checks: list[GateCheck] = []
    if promotion_mode not in PROMOTABLE_MODES:
        if power is not None or outcome is not None:
            raise TournamentError("non-promotable task cannot consume power or outcome evidence")
        disposition = promotion_mode
        terminal = True
    else:
        if coverage == "none":
            raise TournamentError(
                "promotable task must bind at least all selected PredictionCommits"
            )
        if power is None:
            if outcome is not None:
                raise TournamentError("outcomes cannot open without a PowerDecision")
            disposition = "downgraded_without_unblinding"
            terminal = True
        elif power.power_result.get("passed") is not True:
            if outcome is not None:
                raise TournamentError("underpowered task cannot open sealed outcomes")
            disposition = "downgraded_without_unblinding"
            terminal = True
        elif outcome is None:
            _verify_exact_task_units(commits=commits)
            disposition = "ready_for_sealed_inference"
            terminal = False
        else:
            assert outcome_binding is not None and outcome_root is not None
            assert consumption is not None
            _validate_opened_sources(
                lock=lock,
                lock_binding=lock_binding,
                task_id=task_id,
                promotion=promotion,
                commits=commits,
                power=power,
                outcome=outcome,
                outcome_binding=outcome_binding,
                outcome_root=outcome_root,
                consumption=consumption,
            )
            if metric_bundle is None:
                checks.append(
                    GateCheck(
                        "sealed_metric_bundle",
                        False,
                        None,
                        {
                            "required_schema": SEALED_METRIC_BUNDLE_SCHEMA_VERSION,
                            "n_resamples": _POWER_BOOTSTRAP_RESAMPLES,
                        },
                        "promotion is blocked until task-native sealed metrics are frozen "
                        "and recursively rederived from the consumed outcome chain",
                    )
                )
                disposition = "terminal_failure_missing_sealed_metric_bundle"
            else:
                assert metric_binding is not None
                assert multiplicity_bundle is not None
                assert multiplicity_binding is not None
                expected_metric_sources = {
                    "selection_lock_binding": lock_binding,
                    "prediction_commit_bindings": list(commit_bindings),
                    "power_decision_binding": power_binding,
                    "sealed_outcome_binding": outcome_binding,
                    "outcome_consumption_binding": consumption_binding,
                }
                for field, expected in expected_metric_sources.items():
                    if metric_bundle.get(field) != expected:
                        raise TournamentError(
                            f"sealed metric bundle {field} differs from champion gate"
                        )
                if metric_bundle.get("task_id") != task_id:
                    raise TournamentError("sealed metric bundle binds a different task")
                if multiplicity_bundle.get("selection_lock_binding") != lock_binding:
                    raise TournamentError("multiplicity bundle binds a different SelectionLock")
                task_metric_bindings = [
                    item
                    for item in multiplicity_bundle.get("sealed_metric_bindings", ())
                    if isinstance(item, Mapping) and item.get("task_id") == task_id
                ]
                if task_metric_bindings != [metric_binding]:
                    raise TournamentError(
                        "multiplicity bundle does not bind this exact sealed metric bundle"
                    )
                adjusted = multiplicity_bundle.get("adjusted_p_by_task")
                if not isinstance(adjusted, Mapping) or task_id not in adjusted:
                    raise TournamentError("multiplicity bundle lacks this task")
                metric_values = metric_bundle.get("metrics")
                if not isinstance(metric_values, Mapping):
                    raise TournamentError("sealed metric bundle metrics must be an object")
                gate_metrics = {
                    **dict(metric_values),
                    "holm_adjusted_p": adjusted[task_id],
                }
                checks.extend(_track_checks(task_id, gate_metrics, promotion))
                seed_deltas = metric_values.get("seed_deltas")
                expected_seed_ids = {str(seed) for seed in decision["seeds"]}
                valid_seed_deltas = (
                    isinstance(seed_deltas, Mapping)
                    and set(seed_deltas) == expected_seed_ids
                    and all(
                        not isinstance(value, bool)
                        and isinstance(value, (int, float))
                        and isfinite(float(value))
                        for value in seed_deltas.values()
                    )
                )
                positive_count = (
                    sum(float(value) > 0.0 for value in seed_deltas.values())
                    if valid_seed_deltas
                    else None
                )
                checks.append(
                    GateCheck(
                        "positive_seed_direction",
                        positive_count is not None and positive_count >= 4,
                        {
                            "positive_seed_count": positive_count,
                            "seed_deltas": seed_deltas,
                        },
                        {"exact_seed_ids": sorted(expected_seed_ids), ">=": 4},
                        "the primary gain direction must be positive in at least four of five seeds",
                    )
                )
                disposition = (
                    "promoted"
                    if all(check.passed for check in checks)
                    else "terminal_failure_gate_not_met"
                )
            terminal = True
    return {
        "selection_lock": lock,
        "selection_lock_binding": lock_binding,
        "task_decision": decision,
        "promotion": promotion,
        "prediction_commit_bindings": list(commit_bindings),
        "prediction_commits": commits,
        "prediction_commit_roots": commit_roots,
        "prediction_coverage": coverage,
        "power": power,
        "power_binding": power_binding,
        "power_root": power_root,
        "power_evidence": evidence,
        "outcome": outcome,
        "outcome_binding": outcome_binding,
        "outcome_root": outcome_root,
        "consumption": consumption,
        "consumption_binding": consumption_binding,
        "consumption_marker": consumption_marker,
        "sealed_metric_bundle": metric_bundle,
        "sealed_metric_binding": metric_binding,
        "sealed_metric_root": metric_root,
        "confirmatory_multiplicity_bundle": multiplicity_bundle,
        "confirmatory_multiplicity_binding": multiplicity_binding,
        "confirmatory_multiplicity_root": multiplicity_root,
        "checks": checks,
        "disposition": disposition,
        "terminal": terminal,
    }


def _champion_gate_identity(sources: Mapping[str, Any]) -> dict[str, Any]:
    lock = sources["selection_lock"]
    decision = sources["task_decision"]
    outcome = sources["outcome"]
    consumption = sources["consumption"]
    checks = sources["checks"]
    return {
        "schema_version": GATE_DECISION_SCHEMA_VERSION,
        "selection_lock_binding": dict(sources["selection_lock_binding"]),
        "task_id": decision["task_id"],
        "selected_model_id": decision["selected_model_id"],
        "selected_run_ids": sorted(_locked_run_ids(decision, "selected")),
        "baseline_model_id": decision["baseline_model_id"],
        "baseline_run_ids": sorted(_locked_run_ids(decision, "baseline")),
        "promotion_mode": sources["promotion"]["promotion_mode"],
        "promotion_gate_id": decision["promotion_gate"],
        "prediction_coverage": sources["prediction_coverage"],
        "prediction_commit_bindings": list(sources["prediction_commit_bindings"]),
        "power_decision_binding": (
            dict(sources["power_binding"])
            if sources["power_binding"] is not None
            else None
        ),
        "sealed_outcome_binding": (
            dict(sources["outcome_binding"])
            if sources["outcome_binding"] is not None
            else None
        ),
        "outcome_consumption_binding": (
            dict(sources["consumption_binding"])
            if sources["consumption_binding"] is not None
            else None
        ),
        "sealed_metric_bundle_binding": (
            dict(sources["sealed_metric_binding"])
            if sources["sealed_metric_binding"] is not None
            else None
        ),
        "confirmatory_multiplicity_binding": (
            dict(sources["confirmatory_multiplicity_binding"])
            if sources["confirmatory_multiplicity_binding"] is not None
            else None
        ),
        "checks": [check.to_dict() for check in checks],
        "failures": [check.name for check in checks if not check.passed],
        "disposition": sources["disposition"],
        "passed": sources["disposition"] == "promoted",
        "terminal": sources["terminal"],
        "sealed": outcome is not None,
        "holdout_consumed": consumption is not None,
        "sealed_outcome_bundle_id": outcome.bundle_id if outcome is not None else None,
        "selection_lock_id": lock.lock_id,
        "terminal_policy": LOCKED_TERMINAL_POLICY,
    }


def freeze_champion_gate_decision(
    *,
    selection_lock_dir: str | Path,
    task_id: str,
    output_root: str | Path,
    prediction_commit_dirs: Iterable[str | Path] = (),
    power_decision_dir: str | Path | None = None,
    sealed_outcome_bundle_dir: str | Path | None = None,
    outcome_consumption_path: str | Path | None = None,
    sealed_metric_bundle_dir: str | Path | None = None,
    confirmatory_multiplicity_bundle_dir: str | Path | None = None,
) -> Path:
    """Freeze one task check from recursively verified concrete source output files."""

    sources = _champion_gate_sources(
        selection_lock_dir=selection_lock_dir,
        task_id=task_id,
        prediction_commit_dirs=prediction_commit_dirs,
        power_decision_dir=power_decision_dir,
        sealed_outcome_bundle_dir=sealed_outcome_bundle_dir,
        outcome_consumption_path=outcome_consumption_path,
        sealed_metric_bundle_dir=sealed_metric_bundle_dir,
        confirmatory_multiplicity_bundle_dir=confirmatory_multiplicity_bundle_dir,
    )
    identity = _champion_gate_identity(sources)
    payload = {"decision_sha256": canonical_sha256(identity), **identity}
    output = _safe_resolve(output_root, "champion gate output root")
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"champion-gate--{task_id}--{payload['decision_sha256']}"
    if target.exists():
        raise TournamentError(f"champion gate decision already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".champion-gate.", dir=output))
    try:
        write_json_exclusive(staging / "champion_gate_decision.json", payload, mode=0o440)
        freeze_tree(
            staging,
            {
                "artifact_class": "champion_gate_decision",
                "decision_sha256": payload["decision_sha256"],
                "task_id": task_id,
            },
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_champion_gate_decision(target, reverify_sources=True)
    return target


def verify_champion_gate_decision(
    path: str | Path, *, reverify_sources: bool = True
) -> dict[str, Any]:
    """Verify and optionally recursively rederive one frozen best-model check."""

    root = _safe_resolve(path, "champion gate decision")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "champion_gate_decision.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid champion gate decision: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != GATE_DECISION_SCHEMA_VERSION:
        raise TournamentError("unsupported champion gate decision schema")
    claimed = _sha256_identifier(payload.get("decision_sha256"), "decision_sha256")
    identity = dict(payload)
    identity.pop("decision_sha256", None)
    if canonical_sha256(identity) != claimed:
        raise TournamentError("champion gate decision identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "champion_gate_decision",
        "decision_sha256": claimed,
        "task_id": payload.get("task_id"),
    }:
        raise TournamentError("champion gate manifest metadata mismatch")
    if reverify_sources:
        lock_binding = payload.get("selection_lock_binding")
        if not isinstance(lock_binding, Mapping):
            raise TournamentError("champion gate lacks a SelectionLock binding")
        lock_root = _verify_bound_frozen_document(
            {
                key: lock_binding.get(key)
                for key in ("path", "manifest_sha256", "document_sha256")
            },
            filename="selection_lock.json",
            artifact_class="selection_lock",
        )
        commit_bindings = payload.get("prediction_commit_bindings")
        if not isinstance(commit_bindings, list):
            raise TournamentError("champion gate prediction bindings must be an array")
        commit_paths = [binding.get("path") for binding in commit_bindings if isinstance(binding, Mapping)]
        if len(commit_paths) != len(commit_bindings):
            raise TournamentError("champion gate prediction binding has the wrong schema")
        power_binding = payload.get("power_decision_binding")
        outcome_binding = payload.get("sealed_outcome_binding")
        consumption_binding = payload.get("outcome_consumption_binding")
        metric_binding = payload.get("sealed_metric_bundle_binding")
        multiplicity_binding = payload.get("confirmatory_multiplicity_binding")
        sources = _champion_gate_sources(
            selection_lock_dir=lock_root,
            task_id=str(payload.get("task_id")),
            prediction_commit_dirs=commit_paths,
            power_decision_dir=(
                power_binding.get("path") if isinstance(power_binding, Mapping) else None
            ),
            sealed_outcome_bundle_dir=(
                outcome_binding.get("path") if isinstance(outcome_binding, Mapping) else None
            ),
            outcome_consumption_path=(
                consumption_binding.get("path")
                if isinstance(consumption_binding, Mapping)
                else None
            ),
            sealed_metric_bundle_dir=(
                metric_binding.get("path")
                if isinstance(metric_binding, Mapping)
                else None
            ),
            confirmatory_multiplicity_bundle_dir=(
                multiplicity_binding.get("path")
                if isinstance(multiplicity_binding, Mapping)
                else None
            ),
        )
        expected = {"decision_sha256": claimed, **_champion_gate_identity(sources)}
        if dict(payload) != expected:
            raise TournamentError("champion gate decision does not rederive from sources")
    return dict(payload)


def _gate_binding(path: str | Path) -> tuple[dict[str, Any], dict[str, Any], Path]:
    root = _safe_resolve(path, "champion gate decision")
    gate = verify_champion_gate_decision(root, reverify_sources=True)
    binding, _ = _frozen_document_binding(
        root,
        filename="champion_gate_decision.json",
        artifact_class="champion_gate_decision",
    )
    return gate, {
        **binding,
        "decision_sha256": gate["decision_sha256"],
        "task_id": gate["task_id"],
    }, root


def _terminal_decision_summary(gate: Mapping[str, Any]) -> dict[str, Any]:
    power = gate.get("power_decision_binding")
    outcome = gate.get("sealed_outcome_binding")
    consumption = gate.get("outcome_consumption_binding")
    predictions = gate.get("prediction_commit_bindings", ())
    return {
        "decision_sha256": gate["decision_sha256"],
        "task_id": gate["task_id"],
        "selected_model_id": gate["selected_model_id"],
        "model_disposition": gate["promotion_mode"],
        "task_disposition": gate["disposition"],
        "sealed": gate["sealed"],
        "open_champion": gate["promotion_mode"] in PROMOTABLE_MODES,
        "gate_passed": gate["passed"],
        "promoted": gate["disposition"] == "promoted",
        "terminal": gate["terminal"],
        "prediction_bindings": [
            {
                "run_id": item["run_id"],
                "prediction_commit_sha256": item["prediction_commit_sha256"],
            }
            for item in predictions
        ],
        "power_decision_sha256": (
            power["power_decision_sha256"] if isinstance(power, Mapping) else None
        ),
        "outcome_bundle_sha256": (
            outcome["outcome_bundle_sha256"] if isinstance(outcome, Mapping) else None
        ),
        "consumption_sha256": (
            consumption["consumption_sha256"]
            if isinstance(consumption, Mapping)
            else None
        ),
    }


def _terminal_authorization_sources(
    *,
    selection_lock_dir: str | Path,
    champion_gate_decision_dirs: Iterable[str | Path],
) -> dict[str, Any]:
    lock, lock_binding, lock_root = _selection_lock_binding(selection_lock_dir)
    gate_paths = tuple(
        _safe_resolve(path, "champion gate decision")
        for path in champion_gate_decision_dirs
    )
    if not gate_paths or len(gate_paths) != len(set(gate_paths)):
        raise TournamentError("terminal authorization requires distinct gate directories")
    gate_records = [_gate_binding(path) for path in gate_paths]
    gate_records.sort(key=lambda item: str(item[0]["task_id"]))
    gates = [item[0] for item in gate_records]
    task_ids = [str(gate["task_id"]) for gate in gates]
    expected_tasks = sorted(str(item["task_id"]) for item in lock.task_decisions)
    if task_ids != expected_tasks:
        raise TournamentError(
            "terminal authorization gates must cover every SelectionLock task exactly"
        )
    for gate in gates:
        if (
            gate["selection_lock_id"] != lock.lock_id
            or gate["selection_lock_binding"] != lock_binding
        ):
            raise TournamentError("terminal gate uses a different SelectionLock")
        if gate.get("terminal") is not True:
            raise TournamentError(
                f"terminal authorization rejects nonterminal gate for {gate['task_id']}: "
                f"{gate.get('disposition')!r}"
            )

    opened_by_marker: dict[str, set[str]] = {}
    for gate in gates:
        marker = gate.get("outcome_consumption_binding")
        if isinstance(marker, Mapping):
            opened_by_marker.setdefault(str(marker["path"]), set()).add(
                str(gate["task_id"])
            )
    for marker_path, opened_tasks in opened_by_marker.items():
        consumption = verify_outcome_consumption(marker_path)
        marker_tasks = {str(item["task_id"]) for item in consumption.power_bindings}
        if opened_tasks != marker_tasks:
            raise TournamentError(
                "every task in one consumed outcome marker needs a recursively verified gate"
            )
        if len(marker_tasks) > 1 and marker_tasks != {CELL_TASK, VARIANT_TASK}:
            raise TournamentError(
                "only the prespecified joint GSE289173 cell/variant pair may share consumption"
            )
    return {
        "lock": lock,
        "lock_binding": lock_binding,
        "lock_root": lock_root,
        "gates": gates,
        "gate_bindings": [item[1] for item in gate_records],
        "gate_roots": [item[2] for item in gate_records],
        "decisions": [_terminal_decision_summary(gate) for gate in gates],
    }


def _terminal_authorization_identity(sources: Mapping[str, Any]) -> dict[str, Any]:
    lock = sources["lock"]
    return {
        "schema_version": TERMINAL_AUTHORIZATION_SCHEMA_VERSION,
        "selection_lock_binding": dict(sources["lock_binding"]),
        "selection_lock_id": lock.lock_id,
        "selection_lock_manifest_sha256": sources["lock_binding"][
            "manifest_sha256"
        ],
        "candidate_manifest_sha256": lock.candidate_manifest_sha256,
        "plan_sha256": lock.plan_sha256,
        "gate_decision_bindings": list(sources["gate_bindings"]),
        "decisions": list(sources["decisions"]),
        "terminal_policy": LOCKED_TERMINAL_POLICY,
    }


def freeze_terminal_evaluation_authorization(
    *,
    selection_lock_dir: str | Path,
    champion_gate_decision_dirs: Iterable[str | Path],
    output_root: str | Path,
) -> Path:
    """Freeze the aggregate task-complete authorization/release decision ledger."""

    sources = _terminal_authorization_sources(
        selection_lock_dir=selection_lock_dir,
        champion_gate_decision_dirs=champion_gate_decision_dirs,
    )
    identity = _terminal_authorization_identity(sources)
    payload = {"authorization_id": canonical_sha256(identity), **identity}
    output = _safe_resolve(output_root, "terminal authorization output root")
    output.mkdir(parents=True, exist_ok=True)
    target = output / f"terminal-authorization--{payload['authorization_id']}"
    if target.exists():
        raise TournamentError(f"terminal authorization already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".terminal-authorization.", dir=output))
    try:
        write_json_exclusive(
            staging / "terminal_evaluation_authorization.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "terminal_evaluation_authorization",
                "authorization_id": payload["authorization_id"],
            },
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_terminal_evaluation_authorization(target, reverify_sources=True)
    return target


def _absolute_source_paths(value: Any) -> set[Path]:
    result: set[Path] = set()
    if isinstance(value, Mapping):
        for key, item in value.items():
            if isinstance(item, str) and (
                key == "path" or key.endswith("_path") or key.endswith("_dir")
            ):
                candidate = Path(item)
                if candidate.is_absolute():
                    result.add(_safe_resolve(candidate, f"protected source {key}"))
            else:
                result.update(_absolute_source_paths(item))
    elif isinstance(value, (list, tuple)):
        for item in value:
            result.update(_absolute_source_paths(item))
    return result


def _terminal_protected_source_roots(sources: Mapping[str, Any]) -> list[str]:
    roots: set[Path] = {sources["lock_root"], *sources["gate_roots"]}
    lock = sources["lock"]
    candidate = lock.metadata.get("candidate_binding")
    if not isinstance(candidate, Mapping):
        raise TournamentError("SelectionLock lacks its protected candidate binding")
    candidate_root = _safe_resolve(candidate.get("path"), "SelectionLock candidate")
    artifact_class = candidate.get("artifact_class")
    if artifact_class == "selection_candidate_ledger":
        ledger = verify_selection_candidate_ledger(candidate_root)
        if (
            sha256_file(candidate_root / "ARTIFACTS.json")
            != candidate.get("manifest_sha256")
            or ledger["ledger_id"] != candidate.get("source_identity_sha256")
        ):
            raise TournamentError("SelectionLock candidate ledger binding changed")
        roots.add(candidate_root)
        roots.update(_absolute_source_paths(ledger))
    elif artifact_class == "candidate_plan":
        _, plan, binding = _candidate_plan_source(candidate_root)
        if (
            binding["manifest_sha256"] != candidate.get("manifest_sha256")
            or plan["plan_sha256"] != candidate.get("source_identity_sha256")
        ):
            raise TournamentError("SelectionLock candidate plan binding changed")
        roots.add(candidate_root)
    else:
        raise TournamentError("SelectionLock candidate binding has unknown artifact class")
    for gate in sources["gates"]:
        roots.update(_absolute_source_paths(gate))
        for binding in gate.get("prediction_commit_bindings", ()):
            commit = verify_prediction_commit(binding["path"], reverify_sources=True)
            roots.add(_safe_resolve(commit.prediction_bundle_path, "PredictionBundle").parent)
        power = gate.get("power_decision_binding")
        if isinstance(power, Mapping):
            verified_power = verify_power_decision(power["path"], reverify_sources=True)
            roots.add(
                _safe_resolve(
                    verified_power.development_power_evidence_dir,
                    "DevelopmentPowerEvidence",
                )
            )
            evidence = verify_development_power_evidence(
                verified_power.development_power_evidence_dir,
                reverify_sources=True,
            )
            roots.update(_absolute_source_paths(evidence.to_dict()))
    return sorted(path.as_posix() for path in roots)


def verify_terminal_evaluation_authorization(
    path: str | Path, *, reverify_sources: bool = True
) -> dict[str, Any]:
    """Verify the aggregate authorization and return non-identity protected roots."""

    root = _safe_resolve(path, "terminal evaluation authorization")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "terminal_evaluation_authorization.json").read_text(
                encoding="utf-8"
            )
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid terminal evaluation authorization: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != TERMINAL_AUTHORIZATION_SCHEMA_VERSION:
        raise TournamentError("unsupported terminal authorization schema")
    claimed = _sha256_identifier(payload.get("authorization_id"), "authorization_id")
    identity = dict(payload)
    identity.pop("authorization_id", None)
    if canonical_sha256(identity) != claimed:
        raise TournamentError("terminal authorization identity hash mismatch")
    if manifest.get("metadata") != {
        "artifact_class": "terminal_evaluation_authorization",
        "authorization_id": claimed,
    }:
        raise TournamentError("terminal authorization manifest metadata mismatch")
    result = dict(payload)
    protected: list[str] = []
    if reverify_sources:
        lock_binding = payload.get("selection_lock_binding")
        gate_bindings = payload.get("gate_decision_bindings")
        if not isinstance(lock_binding, Mapping) or not isinstance(gate_bindings, list):
            raise TournamentError("terminal authorization source bindings are missing")
        lock_root = _verify_bound_frozen_document(
            {
                key: lock_binding.get(key)
                for key in ("path", "manifest_sha256", "document_sha256")
            },
            filename="selection_lock.json",
            artifact_class="selection_lock",
        )
        gate_paths = [
            binding.get("path")
            for binding in gate_bindings
            if isinstance(binding, Mapping)
        ]
        if len(gate_paths) != len(gate_bindings):
            raise TournamentError("terminal gate binding has the wrong schema")
        sources = _terminal_authorization_sources(
            selection_lock_dir=lock_root,
            champion_gate_decision_dirs=gate_paths,
        )
        expected = {
            "authorization_id": claimed,
            **_terminal_authorization_identity(sources),
        }
        if dict(payload) != expected:
            raise TournamentError("terminal authorization does not rederive from gates")
        protected = _terminal_protected_source_roots(sources)
    result["_protected_source_roots"] = protected
    return result


def _complementarity_evidence_payload(
    evidence: ComplementarityEvidence,
) -> dict[str, Any]:
    if not isinstance(evidence, ComplementarityEvidence):
        raise TournamentError("complementarity_evidence must be ComplementarityEvidence")
    relative = {
        _nonempty_identifier(str(key), "relative-deviance endpoint"): _finite_number(
            value, f"relative_deviance_gains[{key}]"
        )
        for key, value in evidence.relative_deviance_gains.items()
    }
    absolute = {
        _nonempty_identifier(str(key), "absolute-gain endpoint"): _finite_number(
            value, f"absolute_correlation_or_f1_gains[{key}]"
        )
        for key, value in evidence.absolute_correlation_or_f1_gains.items()
    }
    return {
        "residual_correlation": _finite_number(
            evidence.residual_correlation, "complementarity_evidence.residual_correlation"
        ),
        "relative_deviance_gains": dict(sorted(relative.items())),
        "absolute_correlation_or_f1_gains": dict(sorted(absolute.items())),
        "qualifying_seed_count": evidence.qualifying_seed_count,
        "evaluated_seed_count": evidence.evaluated_seed_count,
        "study_count": evidence.study_count,
        "cross_fitted": evidence.cross_fitted,
        "nonnegative_stack": evidence.nonnegative_stack,
        "best_open_models_compared": evidence.best_open_models_compared,
        "derivation_source": evidence.derivation_source,
        "stacking_evidence_binding": evidence.stacking_evidence_binding,
    }


def _selected_shortlist_candidate(
    shortlist: Mapping[str, Any], candidate_id: str
) -> Mapping[str, Any]:
    candidate_id = _sha256_identifier(candidate_id, "candidate_id")
    if candidate_id not in shortlist["selected_candidate_ids"]:
        raise TournamentError("residual candidate is not selected by its frozen shortlist")
    matches = [
        item
        for item in shortlist["candidates"]
        if isinstance(item, Mapping) and item.get("candidate_id") == candidate_id
    ]
    if len(matches) != 1:
        raise TournamentError("residual candidate does not resolve uniquely in shortlist")
    return matches[0]


def _residual_contract_from_shortlist(
    shortlist: Mapping[str, Any], candidate_id: str
) -> dict[str, Any]:
    candidate = _selected_shortlist_candidate(shortlist, candidate_id)
    ledger_binding = shortlist["selection_candidate_ledger_binding"]
    ledger = verify_selection_candidate_ledger(ledger_binding["path"])
    model_matches = [
        item
        for item in ledger["model_candidates"]
        if item.get("candidate_id") == candidate_id
        and item.get("task_id") == shortlist["task_id"]
    ]
    if len(model_matches) != 1:
        raise TournamentError("residual candidate does not resolve in its scientific ledger")
    model = model_matches[0]
    model_id = _nonempty_identifier(model.get("model_id"), "residual model_id")
    if model_id in CONDITIONAL_TRIGGER_INELIGIBLE_MODEL_IDS:
        raise TournamentError(
            f"control or ablation {model_id} cannot supply conditional-model residuals"
        )
    disposition_matches = [
        item
        for item in ledger["model_dispositions"]
        if isinstance(item, Mapping) and item.get("model_id") == model_id
    ]
    if (
        len(disposition_matches) != 1
        or disposition_matches[0].get("champion_eligible") is not True
    ):
        raise TournamentError(
            f"conditional-model residual source {model_id} is not an open eligible model"
        )
    run_ids = tuple(model.get("run_ids", ()))
    science_by_run = {
        item["run_id"]: item for item in ledger["scientific_runs"]
    }
    runs_by_id = {item["run_id"]: item for item in ledger["runs"]}
    science = [science_by_run.get(run_id) for run_id in run_ids]
    runs = [runs_by_id.get(run_id) for run_id in run_ids]
    if not run_ids or any(item is None for item in (*science, *runs)):
        raise TournamentError("residual candidate has incomplete scientific run bindings")
    evaluator_ids = {str(item["endpoint_evaluator_id"]) for item in science}
    evaluator_hashes = {str(item["endpoint_evaluator_sha256"]) for item in science}
    folds = {str(item["fold"]) for item in runs}
    immutable_inputs = {canonical_sha256(item["immutable_inputs"]) for item in runs}
    if not (
        len(evaluator_ids) == len(evaluator_hashes) == len(folds) == len(immutable_inputs) == 1
    ):
        raise TournamentError(
            "residual candidate runs do not share evaluator, fold, and dataset identity"
        )
    inputs = dict(runs[0]["immutable_inputs"])
    dataset_locks = ledger.get("dataset_locks")
    if not isinstance(dataset_locks, Mapping) or not set(inputs).issubset(dataset_locks):
        raise TournamentError("residual candidate lacks frozen dataset locks")
    registry_bindings: dict[str, str] = {}
    for dataset_id in sorted(inputs):
        lock = dataset_locks[dataset_id]
        if not isinstance(lock, Mapping):
            raise TournamentError(f"residual dataset lock {dataset_id!r} is invalid")
        if lock.get("activation_sha256") != inputs[dataset_id]:
            raise TournamentError(
                f"residual run changed dataset activation for {dataset_id}"
            )
        registry_bindings[dataset_id] = _sha256_identifier(
            lock.get("registry_sha256"),
            f"residual dataset {dataset_id} registry_sha256",
        )
    task_dispositions = [
        item
        for item in ledger["task_dispositions"]
        if item.get("task_id") == shortlist["task_id"]
    ]
    if len(task_dispositions) != 1:
        raise TournamentError("residual candidate lacks one frozen task disposition")
    resampling_units = task_dispositions[0].get("resampling_units")
    units = _unique_strings(resampling_units, "resampling_units")
    return {
        "task_id": shortlist["task_id"],
        "candidate_id": candidate_id,
        "source_family_id": candidate["family"],
        "endpoint_evaluator_id": next(iter(evaluator_ids)),
        "endpoint_evaluator_sha256": next(iter(evaluator_hashes)),
        "dataset_ids": sorted(inputs),
        "dataset_registry_sha256s": registry_bindings,
        "fold": next(iter(folds)),
        "row_id_field": "row_hash",
        "unit_id_namespace": f"{shortlist['task_id']}:held_development_residual_v1",
        "biological_unit": "+".join(units),
    }


def _residual_rows(path: Path) -> list[dict[str, str]]:
    rows = _read_exact_tsv(path, ("row_hash", "residual"), "development residual table")
    for index, row in enumerate(rows):
        _finite_number(row["residual"], f"residual row {index}")
    return rows


def freeze_development_residual_bundle(
    *,
    development_shortlist_dir: str | Path,
    candidate_id: str,
    residual_table_path: str | Path,
    output_root: str | Path,
) -> Path:
    """Freeze row-level held residuals under the shortlist's scientific requirements."""

    shortlist_root = _safe_resolve(
        development_shortlist_dir, "development shortlist"
    )
    shortlist = verify_development_shortlist(shortlist_root)
    contract = _residual_contract_from_shortlist(shortlist, candidate_id)
    rows = _residual_rows(_safe_resolve(residual_table_path, "development residual table"))
    shortlist_binding, _ = _frozen_document_binding(
        shortlist_root,
        filename="development_shortlist.json",
        artifact_class="development_shortlist",
    )
    shortlist_binding = {**shortlist_binding, "shortlist_id": shortlist["shortlist_id"]}
    output = _safe_resolve(output_root, "development residual output root")
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".development-residuals.", dir=output))
    try:
        table_path = staging / "residuals.tsv"
        write_text_exclusive(
            table_path,
            _tsv_text(("row_hash", "residual"), rows),
            mode=0o440,
        )
        table_ref = ArtifactRef.from_path(
            table_path,
            relative_to=staging,
            media_type="text/tab-separated-values",
            role=f"held_development_residuals:{contract['task_id']}",
        )
        identity = {
            "schema_version": DEVELOPMENT_RESIDUAL_BUNDLE_SCHEMA_VERSION,
            "development_shortlist_binding": shortlist_binding,
            **contract,
            "residual_table": table_ref.to_dict(),
            "n_rows": len(rows),
            "row_set_sha256": canonical_sha256([row["row_hash"] for row in rows]),
            "created_before_outcome_unblind": True,
            "sealed_results_used": False,
        }
        payload = {"residual_bundle_id": canonical_sha256(identity), **identity}
        write_json_exclusive(
            staging / "development_residual_bundle.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "development_residual_bundle",
                "residual_bundle_id": payload["residual_bundle_id"],
            },
        )
        target = output / f"development-residuals--{payload['residual_bundle_id']}"
        if target.exists():
            raise TournamentError(f"development residual bundle already exists: {target}")
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_development_residual_bundle(target)
    return target


def verify_development_residual_bundle(path: str | Path) -> dict[str, Any]:
    """Recursively verify a row-level held residual bundle."""

    root = _safe_resolve(path, "development residual bundle")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "development_residual_bundle.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid frozen development residual bundle: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != DEVELOPMENT_RESIDUAL_BUNDLE_SCHEMA_VERSION:
        raise TournamentError("unsupported development residual bundle schema")
    binding = payload.get("development_shortlist_binding")
    if not isinstance(binding, Mapping) or set(binding) != {
        "path", "manifest_sha256", "document_sha256", "shortlist_id"
    }:
        raise TournamentError("development residual shortlist binding has the wrong schema")
    shortlist_root = _verify_bound_frozen_document(
        {key: binding[key] for key in ("path", "manifest_sha256", "document_sha256")},
        filename="development_shortlist.json",
        artifact_class="development_shortlist",
    )
    shortlist = verify_development_shortlist(shortlist_root)
    if shortlist["shortlist_id"] != binding["shortlist_id"]:
        raise TournamentError("development residual shortlist binding changed")
    contract = _residual_contract_from_shortlist(shortlist, payload.get("candidate_id"))
    for field, expected in contract.items():
        if payload.get(field) != expected:
            raise TournamentError(f"development residual {field} changed from shortlist")
    table_ref = ArtifactRef.from_dict(payload.get("residual_table"))
    if table_ref.role != f"held_development_residuals:{contract['task_id']}":
        raise TournamentError("development residual table role is invalid")
    rows = _residual_rows(table_ref.validate(root, require_relative=True))
    identity = dict(payload)
    claimed = _sha256_identifier(identity.pop("residual_bundle_id", None), "residual_bundle_id")
    if (
        payload.get("n_rows") != len(rows)
        or payload.get("row_set_sha256")
        != canonical_sha256([row["row_hash"] for row in rows])
        or payload.get("created_before_outcome_unblind") is not True
        or payload.get("sealed_results_used") is not False
        or canonical_sha256(identity) != claimed
    ):
        raise TournamentError("development residual bundle identity does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": "development_residual_bundle",
        "residual_bundle_id": claimed,
    }:
        raise TournamentError("development residual bundle manifest metadata mismatch")
    return dict(payload)


def _pearson_residual_correlation(
    first_root: Path,
    first: Mapping[str, Any],
    second_root: Path,
    second: Mapping[str, Any],
) -> float:
    identity_fields = (
        "task_id",
        "endpoint_evaluator_id",
        "endpoint_evaluator_sha256",
        "dataset_ids",
        "dataset_registry_sha256s",
        "fold",
        "row_id_field",
        "unit_id_namespace",
        "biological_unit",
        "n_rows",
        "row_set_sha256",
    )
    for field in identity_fields:
        if first.get(field) != second.get(field):
            raise TournamentError(
                f"conditional residual bundles differ in frozen {field} identity"
            )
    first_rows = _residual_rows(
        ArtifactRef.from_dict(first["residual_table"]).validate(
            first_root, require_relative=True
        )
    )
    second_rows = _residual_rows(
        ArtifactRef.from_dict(second["residual_table"]).validate(
            second_root, require_relative=True
        )
    )
    if [row["row_hash"] for row in first_rows] != [
        row["row_hash"] for row in second_rows
    ]:
        raise TournamentError("conditional residual row-ID sets are not exactly aligned")
    if len(first_rows) < 3:
        raise TournamentError("Pearson residual correlation requires at least three rows")
    x = [_finite_number(row["residual"], "first residual") for row in first_rows]
    y = [_finite_number(row["residual"], "second residual") for row in second_rows]
    x_mean = fsum(x) / len(x)
    y_mean = fsum(y) / len(y)
    x_centered = [value - x_mean for value in x]
    y_centered = [value - y_mean for value in y]
    denominator = sqrt(
        fsum(value * value for value in x_centered)
        * fsum(value * value for value in y_centered)
    )
    if denominator == 0.0:
        raise TournamentError("Pearson residual correlation is undefined for zero variance")
    correlation = fsum(
        first_value * second_value
        for first_value, second_value in zip(x_centered, y_centered, strict=True)
    ) / denominator
    if not isfinite(correlation) or not -1.0 <= correlation <= 1.0:
        raise TournamentError("Pearson residual correlation is not finite in [-1, 1]")
    return correlation


def _complementarity_decision_identity(
    *,
    triggered: bool,
    evidence_receipt_sha256: str,
    complementarity_evidence_sha256: str,
    pairwise_correlations_sha256: str,
    shortlist_ids: Sequence[str],
    passing_tracks: Sequence[str],
    source_families: Sequence[str],
    selected_candidate_ids: Sequence[str],
    complementary_pairs: Sequence[Sequence[str]],
    created_before_outcome_unblind: bool,
    reasons: Sequence[str],
) -> dict[str, Any]:
    return {
        "schema_version": CONDITIONAL_DECISION_SCHEMA_VERSION,
        "triggered": triggered,
        "evidence_receipt_sha256": evidence_receipt_sha256,
        "complementarity_evidence_sha256": complementarity_evidence_sha256,
        "pairwise_correlations_sha256": pairwise_correlations_sha256,
        "shortlist_ids": list(shortlist_ids),
        "passing_tracks": list(passing_tracks),
        "source_families": list(source_families),
        "selected_candidate_ids": list(selected_candidate_ids),
        "complementary_pairs": [list(pair) for pair in complementary_pairs],
        "created_before_outcome_unblind": created_before_outcome_unblind,
        "terminal_track": CONDITIONAL_CONTEXT_TASK,
        "reasons": list(reasons),
    }


def conditional_new_model_trigger(
    development_shortlist_paths: Iterable[str | Path],
    *,
    development_residual_bundle_paths: Iterable[str | Path],
    complementarity_evidence: ComplementarityEvidence,
    evidence_receipt_sha256: str,
    comparison_task_ids: Sequence[str] = (VARIANT_TASK, RNA_ATAC_TASK),
    terminal_history: Iterable[Any] = (),
    new_model_version: str | None = None,
    holdout_id: str | None = None,
) -> ComplementarityDecision:
    """Evaluate the one-model trigger from frozen pre-selection development output files."""

    evidence_receipt_sha256 = _sha256_identifier(
        evidence_receipt_sha256, "evidence_receipt_sha256"
    )
    try:
        complementarity_evidence.require_verified_artifact()
    except ConditionalModelError as exc:
        raise TournamentError(
            f"invalid frozen complementarity evidence: {exc}"
        ) from exc
    if (
        evidence_receipt_sha256
        != complementarity_evidence.stacking_evidence_binding["manifest_sha256"]
    ):
        raise TournamentError(
            "evidence_receipt_sha256 must equal the verified stacking-evidence "
            "manifest SHA-256"
        )
    comparison = _unique_strings(
        comparison_task_ids, "comparison_task_ids", sorted_required=False
    )
    required_comparison = tuple(sorted((VARIANT_TASK, RNA_ATAC_TASK)))
    if tuple(sorted(comparison)) != required_comparison:
        raise TournamentError(
            "comparison_task_ids must be exactly variant_to_regulation and "
            "rna_conditioned_atac"
        )
    comparison = required_comparison
    if (new_model_version is None) != (holdout_id is None):
        raise TournamentError("new_model_version and holdout_id must be supplied together")
    if new_model_version is not None and holdout_id is not None:
        contextual_history = tuple(
            record
            for record in terminal_history
            if str(_read_first(record, ("task_id", "track"), "")).strip()
            == CONDITIONAL_CONTEXT_TASK
        )
        enforce_terminal_failure_policy(
            contextual_history,
            model_version=new_model_version,
            task_id=CONDITIONAL_CONTEXT_TASK,
            sealed_outcome_bundle_id=holdout_id,
        )

    paths = tuple(
        _safe_resolve(path, "development shortlist")
        for path in development_shortlist_paths
    )
    if len(paths) != 2 or len(set(paths)) != 2:
        raise TournamentError("conditional trigger requires exactly two distinct shortlists")
    shortlists = tuple(verify_development_shortlist(path) for path in paths)
    by_task: dict[str, Mapping[str, Any]] = {}
    path_by_shortlist_id: dict[str, Path] = {}
    reasons: list[str] = []
    shortlist_ids: list[str] = []
    selected_candidate_ids: set[str] = set()
    for path, shortlist in zip(paths, shortlists, strict=True):
        task_id = str(shortlist["task_id"])
        if task_id in by_task:
            raise TournamentError(f"duplicate development shortlist for {task_id!r}")
        by_task[task_id] = shortlist
        shortlist_id = str(shortlist["shortlist_id"])
        shortlist_ids.append(shortlist_id)
        path_by_shortlist_id[shortlist_id] = path
        selected_candidate_ids.update(str(item) for item in shortlist["selected_candidate_ids"])
    if set(by_task) != set(comparison):
        raise TournamentError(
            "development shortlists must exactly match variant_to_regulation and "
            "rna_conditioned_atac"
        )
    residual_paths = tuple(
        _safe_resolve(path, "development residual bundle")
        for path in development_residual_bundle_paths
    )
    if len(residual_paths) != 2 or len(set(residual_paths)) != 2:
        raise TournamentError("conditional trigger requires exactly two residual bundles")
    residuals = tuple(verify_development_residual_bundle(path) for path in residual_paths)
    residual_bundle_ids = tuple(str(item["residual_bundle_id"]) for item in residuals)
    if len(set(residual_bundle_ids)) != 2:
        raise TournamentError("conditional residual bundles must have distinct identities")
    residual_task_ids = {str(item["task_id"]) for item in residuals}
    if len(residual_task_ids) != 1:
        raise TournamentError(
            "conditional residual correlation requires two candidates from the same "
            "comparison task; cross-task residual correlation is prohibited"
        )
    residual_task_id = next(iter(residual_task_ids))
    if residual_task_id not in comparison:
        raise TournamentError("residual bundles do not belong to a comparison task")
    first_binding = residuals[0].get("development_shortlist_binding")
    second_binding = residuals[1].get("development_shortlist_binding")
    if not isinstance(first_binding, Mapping) or not isinstance(second_binding, Mapping):
        raise TournamentError("residual bundles lack frozen shortlist bindings")
    if dict(first_binding) != dict(second_binding):
        raise TournamentError(
            "conditional residual bundles must bind the same frozen development shortlist"
        )
    residual_shortlist_id = str(first_binding.get("shortlist_id", ""))
    bound_shortlist = by_task[residual_task_id]
    if (
        residual_shortlist_id != bound_shortlist["shortlist_id"]
        or residual_shortlist_id not in path_by_shortlist_id
        or _safe_resolve(first_binding.get("path"), "residual shortlist binding")
        != path_by_shortlist_id[residual_shortlist_id]
    ):
        raise TournamentError(
            "residual shortlist does not match either bound comparison-task shortlist"
        )
    residual_candidate_ids = tuple(str(item["candidate_id"]) for item in residuals)
    if len(set(residual_candidate_ids)) != 2:
        raise TournamentError(
            "conditional residual bundles must represent two distinct selected candidates"
        )
    residual_candidates = tuple(
        _selected_shortlist_candidate(bound_shortlist, candidate_id)
        for candidate_id in residual_candidate_ids
    )
    source_families = {
        _nonempty_identifier(str(candidate["family"]), "residual source family")
        for candidate in residual_candidates
    }
    for residual, candidate in zip(residuals, residual_candidates, strict=True):
        if residual.get("source_family_id") != candidate["family"]:
            raise TournamentError(
                "residual candidate source family changed from its bound shortlist"
            )
    if len(source_families) != 2:
        reasons.append(
            "comparison residuals must represent distinct model/source families"
        )
    correlation = _pearson_residual_correlation(
        residual_paths[0],
        residuals[0],
        residual_paths[1],
        residuals[1],
    )
    evidence_payload = _complementarity_evidence_payload(complementarity_evidence)
    qualifying_endpoint_keys = {
        endpoint
        for endpoint, gain in evidence_payload["relative_deviance_gains"].items()
        if gain >= 0.05
    } | {
        endpoint
        for endpoint, gain in evidence_payload[
            "absolute_correlation_or_f1_gains"
        ].items()
        if gain >= 0.02
    }
    missing_comparison_evidence = sorted(
        set(comparison).difference(qualifying_endpoint_keys)
    )
    if missing_comparison_evidence:
        reasons.append(
            "qualifying evidence endpoint keys must include both comparison task IDs; "
            "missing " + ", ".join(missing_comparison_evidence)
        )
    if abs(correlation - float(evidence_payload["residual_correlation"])) > 1e-12:
        reasons.append("pairwise residual correlation does not match frozen evidence")
    elif correlation >= 0.80:
        reasons.append("held residual correlation is not strictly below 0.80")
    try:
        complementarity_evidence.require_trigger()
    except ConditionalModelError as exc:
        reasons.append(str(exc))

    passing_tracks = tuple(sorted(set(comparison).intersection(qualifying_endpoint_keys)))
    frozen_shortlist_ids = tuple(sorted(shortlist_ids))
    frozen_source_families = tuple(sorted(source_families))
    frozen_candidate_ids = tuple(sorted(selected_candidate_ids))
    pair = tuple(sorted(residual_candidate_ids))
    residual_identity_fields = (
        "task_id",
        "endpoint_evaluator_id",
        "endpoint_evaluator_sha256",
        "dataset_ids",
        "dataset_registry_sha256s",
        "fold",
        "row_id_field",
        "unit_id_namespace",
        "biological_unit",
        "n_rows",
        "row_set_sha256",
    )
    residual_evidence = {
        "comparison_task_ids": list(comparison),
        "residual_task_id": residual_task_id,
        "shortlist_id": residual_shortlist_id,
        "candidate_ids": list(pair),
        "source_family_ids": list(frozen_source_families),
        "residual_bundle_ids": sorted(residual_bundle_ids),
        "row_set_sha256": residuals[0]["row_set_sha256"],
        "frozen_identity_sha256": canonical_sha256(
            {field: residuals[0][field] for field in residual_identity_fields}
        ),
        "pearson_residual_correlation": correlation,
    }
    complementary_pairs = (pair,) if correlation < 0.80 and not reasons else ()
    identity = _complementarity_decision_identity(
        triggered=not reasons,
        evidence_receipt_sha256=evidence_receipt_sha256,
        complementarity_evidence_sha256=canonical_sha256(evidence_payload),
        pairwise_correlations_sha256=canonical_sha256(residual_evidence),
        shortlist_ids=frozen_shortlist_ids,
        passing_tracks=passing_tracks,
        source_families=frozen_source_families,
        selected_candidate_ids=frozen_candidate_ids,
        complementary_pairs=complementary_pairs,
        created_before_outcome_unblind=True,
        reasons=tuple(reasons),
    )
    return ComplementarityDecision(
        decision_id=canonical_sha256(identity),
        triggered=not reasons,
        evidence_receipt_sha256=evidence_receipt_sha256,
        complementarity_evidence_sha256=canonical_sha256(evidence_payload),
        pairwise_correlations_sha256=canonical_sha256(residual_evidence),
        shortlist_ids=frozen_shortlist_ids,
        passing_tracks=passing_tracks,
        source_families=frozen_source_families,
        selected_candidate_ids=frozen_candidate_ids,
        complementary_pairs=complementary_pairs,
        created_before_outcome_unblind=True,
        reasons=tuple(reasons),
    )


def freeze_conditional_model_decision(
    *,
    development_shortlist_paths: Iterable[str | Path],
    development_residual_bundle_paths: Iterable[str | Path],
    complementarity_evidence: ComplementarityEvidence,
    evidence_receipt_sha256: str,
    output_root: str | Path,
    comparison_task_ids: Sequence[str] = (VARIANT_TASK, RNA_ATAC_TASK),
) -> Path:
    """Freeze the pre-selection conditional trigger and all shortlist bindings."""

    paths = tuple(
        _safe_resolve(path, "development shortlist")
        for path in development_shortlist_paths
    )
    residual_paths = tuple(
        _safe_resolve(path, "development residual bundle")
        for path in development_residual_bundle_paths
    )
    decision = conditional_new_model_trigger(
        paths,
        development_residual_bundle_paths=residual_paths,
        complementarity_evidence=complementarity_evidence,
        evidence_receipt_sha256=evidence_receipt_sha256,
        comparison_task_ids=comparison_task_ids,
    )
    bindings = []
    for path in sorted(paths):
        shortlist = verify_development_shortlist(path)
        bindings.append(
            {
                "path": path.as_posix(),
                "manifest_sha256": sha256_file(path / "ARTIFACTS.json"),
                "document_sha256": sha256_file(path / "development_shortlist.json"),
                "shortlist_id": shortlist["shortlist_id"],
            }
        )
    residual_bindings = []
    for path in sorted(residual_paths):
        residual = verify_development_residual_bundle(path)
        binding, _ = _frozen_document_binding(
            path,
            filename="development_residual_bundle.json",
            artifact_class="development_residual_bundle",
        )
        residual_bindings.append(
            {**binding, "residual_bundle_id": residual["residual_bundle_id"]}
        )
    payload = {
        **decision.to_dict(),
        "schema_version": CONDITIONAL_DECISION_SCHEMA_VERSION,
        "terminal_track": CONDITIONAL_CONTEXT_TASK,
        "comparison_task_ids": sorted(comparison_task_ids),
        "development_shortlist_bindings": bindings,
        "complementarity_evidence": _complementarity_evidence_payload(
            complementarity_evidence
        ),
        "development_residual_bindings": residual_bindings,
    }
    root = _safe_resolve(output_root, "conditional decision output root")
    root.mkdir(parents=True, exist_ok=True)
    target = root / f"conditional-decision--{decision.decision_id}"
    if target.exists():
        raise TournamentError(f"conditional decision already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=".conditional-decision.", dir=root))
    try:
        write_json_exclusive(staging / "conditional_decision.json", payload, mode=0o440)
        freeze_tree(
            staging,
            {"artifact_class": "conditional_decision", "decision_id": decision.decision_id},
        )
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_conditional_model_decision(target)
    return target


def verify_conditional_model_decision(path: str | Path) -> dict[str, Any]:
    """Recursively rederive a frozen pre-selection conditional trigger."""

    root = _safe_resolve(path, "conditional decision")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads((root / "conditional_decision.json").read_text(encoding="utf-8"))
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TournamentError(f"invalid frozen conditional decision: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise TournamentError("conditional decision must be a JSON object")
    if payload.get("schema_version") != CONDITIONAL_DECISION_SCHEMA_VERSION:
        raise TournamentError("unsupported conditional decision schema_version")
    bindings = payload.get("development_shortlist_bindings")
    if not isinstance(bindings, list) or len(bindings) != 2:
        raise TournamentError("conditional decision must bind exactly two shortlists")
    paths: list[Path] = []
    for index, raw in enumerate(bindings):
        if not isinstance(raw, Mapping) or set(raw) != {
            "path",
            "manifest_sha256",
            "document_sha256",
            "shortlist_id",
        }:
            raise TournamentError(f"conditional shortlist binding {index} has the wrong schema")
        source = Path(_nonempty_identifier(raw["path"], f"shortlist binding {index}.path"))
        if not source.is_absolute() or source.is_symlink() or not source.is_dir():
            raise TournamentError("conditional shortlist bindings require absolute frozen directories")
        source = _safe_resolve(source, "conditional shortlist binding")
        shortlist = verify_development_shortlist(source)
        if (
            sha256_file(source / "ARTIFACTS.json")
            != _sha256_identifier(raw["manifest_sha256"], "manifest_sha256")
            or sha256_file(source / "development_shortlist.json")
            != _sha256_identifier(raw["document_sha256"], "document_sha256")
            or shortlist["shortlist_id"] != raw["shortlist_id"]
        ):
            raise TournamentError("conditional shortlist binding changed")
        paths.append(source)
    evidence_raw = payload.get("complementarity_evidence")
    if not isinstance(evidence_raw, Mapping):
        raise TournamentError("conditional complementarity evidence must be an object")
    evidence_binding = evidence_raw.get("stacking_evidence_binding")
    if not isinstance(evidence_binding, Mapping) or set(evidence_binding) != {
        "path",
        "manifest_sha256",
        "document_sha256",
    }:
        raise TournamentError(
            "conditional complementarity evidence lacks its frozen binding"
        )
    try:
        evidence = ComplementarityEvidence.from_stacking_evidence(
            evidence_binding["path"]
        )
    except (ConditionalModelError, KeyError, TypeError) as exc:
        raise TournamentError(
            f"conditional complementarity evidence is invalid: {exc}"
        ) from exc
    if evidence.stacking_evidence_binding != dict(evidence_binding):
        raise TournamentError("conditional stacking-evidence binding changed")
    residual_bindings = payload.get("development_residual_bindings")
    if not isinstance(residual_bindings, list) or len(residual_bindings) != 2:
        raise TournamentError("conditional decision must bind exactly two residual bundles")
    residual_paths: list[Path] = []
    for index, raw in enumerate(residual_bindings):
        if not isinstance(raw, Mapping) or set(raw) != {
            "path", "manifest_sha256", "document_sha256", "residual_bundle_id"
        }:
            raise TournamentError(f"conditional residual binding {index} has the wrong schema")
        source = _verify_bound_frozen_document(
            {key: raw[key] for key in ("path", "manifest_sha256", "document_sha256")},
            filename="development_residual_bundle.json",
            artifact_class="development_residual_bundle",
        )
        residual = verify_development_residual_bundle(source)
        if residual["residual_bundle_id"] != raw["residual_bundle_id"]:
            raise TournamentError("conditional residual binding changed")
        residual_paths.append(source)
    derived = conditional_new_model_trigger(
        paths,
        development_residual_bundle_paths=residual_paths,
        complementarity_evidence=evidence,
        evidence_receipt_sha256=payload.get("evidence_receipt_sha256"),
        comparison_task_ids=payload.get("comparison_task_ids", ()),
    )
    expected = {
        **derived.to_dict(),
        "schema_version": CONDITIONAL_DECISION_SCHEMA_VERSION,
        "terminal_track": CONDITIONAL_CONTEXT_TASK,
        "comparison_task_ids": sorted(payload["comparison_task_ids"]),
        "development_shortlist_bindings": bindings,
        "complementarity_evidence": _complementarity_evidence_payload(evidence),
        "development_residual_bindings": residual_bindings,
    }
    if dict(payload) != expected:
        raise TournamentError("conditional decision does not rederive from frozen evidence")
    if manifest.get("metadata") != {
        "artifact_class": "conditional_decision",
        "decision_id": derived.decision_id,
    }:
        raise TournamentError("conditional decision manifest metadata mismatch")
    return dict(payload)


__all__ = [
    "BLOCKED_MISSING_EXTERNAL_FAMILY",
    "BLOCKED_MISSING_SEALED_SOURCE",
    "BULK_TASK",
    "CELL_TASK",
    "CHAMPION_GATE_THRESHOLDS",
    "CONDITIONAL_CONTEXT_TASK",
    "CONDITIONAL_SEALED",
    "ComplementarityDecision",
    "DEVELOPMENT_ENSEMBLE_POLICY_ID",
    "DEVELOPMENT_ONLY",
    "DevelopmentGateDecision",
    "ENSEMBLE_SCORE_SCALE_POLICY",
    "EXPLORATORY_ONLY",
    "GRAPH_TASK",
    "GateCheck",
    "LEDGER_AUTHORITY_FINALIST",
    "LEDGER_AUTHORITY_SCREENING",
    "LOCKED_MULTIPLICITY_PLAN",
    "Objective",
    "PERTURBATION_TASK",
    "PROMOTION_GATE_POLICY",
    "RNA_ATAC_TASK",
    "SCREENING_CAMPAIGN_UNIVERSE_SCHEMA_VERSION",
    "SEALED_CONFIRMATORY",
    "SHORTLIST_SOURCE_WAVE_CONTRACTS",
    "SelectionLockError",
    "TerminalFailureError",
    "TournamentError",
    "UNIFIED_TASK",
    "UNIVERSE_KIND_FINALIST",
    "UNIVERSE_KIND_SCREENING",
    "VARIANT_TASK",
    "campaign_universe_sha256",
    "conditional_new_model_trigger",
    "enforce_terminal_failure_policy",
    "family_model_slots",
    "freeze_campaign_universe",
    "freeze_champion_gate_decision",
    "freeze_conditional_model_decision",
    "freeze_confirmatory_multiplicity_bundle",
    "freeze_development_outcome_bundle",
    "freeze_development_residual_bundle",
    "freeze_development_shortlist",
    "freeze_finalist_campaign_universe",
    "freeze_finalist_candidate_ledger",
    "freeze_finalist_development_metric_bundle",
    "freeze_scientific_run_receipt",
    "freeze_screening_campaign_universe",
    "freeze_screening_candidate_ledger",
    "freeze_sealed_metric_bundle",
    "freeze_terminal_evaluation_authorization",
    "freeze_variant_secondary_run_receipt",
    "one_standard_error_candidates",
    "pareto_frontier",
    "require_selection_lock",
    "select_family_top_two",
    "verify_campaign_universe",
    "verify_champion_gate_decision",
    "verify_conditional_model_decision",
    "verify_confirmatory_multiplicity_bundle",
    "verify_development_outcome_bundle",
    "verify_development_residual_bundle",
    "verify_development_shortlist",
    "verify_finalist_campaign_universe",
    "verify_finalist_development_metric_bundle",
    "verify_scientific_run_receipt",
    "verify_screening_campaign_universe",
    "verify_sealed_metric_bundle",
    "verify_selection_candidate_ledger",
    "verify_terminal_evaluation_authorization",
    "verify_variant_secondary_run_receipt",
]
