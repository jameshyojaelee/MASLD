"""Task-native sealed metrics rederived from five-seed frozen predictions.

Seeds define one ensemble and one stability diagnostic.  Inferential intervals
resample donors or LD blocks, never seeds or rows as biological replicates.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from math import floor, isfinite
from random import Random
from typing import Any, Mapping, Sequence

from ..hashing import canonical_sha256
from .metrics import (
    average_precision,
    donor_class_balanced_weights,
    fisher_z_mean,
    multiclass_brier_score,
    spearman_correlation,
)


CELL_TASK = "cell_state_mapping"
VARIANT_TASK = "variant_to_regulation"
SEALED_METRIC_EVALUATOR_ID = "five_seed_task_native_sealed_metrics_v1"
_VARIANT_CAPABILITY_BINDING_SCHEMA_VERSION = (
    "masld-bench-locked-variant-primary-capability-v1"
)
_VARIANT_PRIMARY_ENDPOINT_ID = "signed_cell_type_eqtl_effect"
_VARIANT_PRIMARY_EVALUATOR_ID = (
    "variant_ld_block_fisher_z_spearman_gain_v1"
)
_VARIANT_SECONDARY_BINDING_SCHEMA_VERSION = (
    "masld-bench-locked-variant-secondary-evaluation-v1"
)
_VARIANT_RETRIEVAL_ENDPOINT_ID = "eqtl_retrieval"
_VARIANT_RETRIEVAL_EVALUATOR_ID = "variant_ld_block_eqtl_retrieval_auprc_v1"
_ABSOLUTE_SIGNED_EFFECT_TRANSFORM_ID = "absolute_signed_effect_v1"
_IDENTITY_LINK_SCORE_TRANSFORM_ID = "identity_link_score_v1"
_MANDATORY_RETRIEVAL_MODELS = ("abc", "nearest_gene", "re2g")


class SealedMetricError(ValueError):
    """Raised when sealed sources cannot support an exact registered metric."""


def _locked_string_array(value: Any) -> list[str] | None:
    """Read a locked string array that may be a list or a contract tuple.

    These capability records come out of a ``SelectionLock`` task decision,
    and ``contracts._as_metadata`` freezes every nested array there into a
    tuple.  Testing ``isinstance(value, list)`` or comparing against a list
    literal would therefore reject a perfectly valid locked record.
    """

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        return None
    items = list(value)
    if any(not isinstance(item, str) or not item for item in items):
        return None
    return items


@dataclass(frozen=True, slots=True)
class SealedMetricResult:
    task_id: str
    evaluator_id: str
    metrics: Mapping[str, Any]
    ensemble_rows: tuple[Mapping[str, str], ...]
    n_resamples: int
    bootstrap_seed: int
    independent_unit: str
    independent_unit_count: int
    row_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "task_id": self.task_id,
            "evaluator_id": self.evaluator_id,
            "metrics": dict(self.metrics),
            "ensemble_rows": [dict(row) for row in self.ensemble_rows],
            "n_resamples": self.n_resamples,
            "bootstrap_seed": self.bootstrap_seed,
            "independent_unit": self.independent_unit,
            "independent_unit_count": self.independent_unit_count,
            "row_count": self.row_count,
        }


def _validate_variant_primary_capability(
    value: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise SealedMetricError(
            "variant sealed metrics require a locked primary-capability binding"
        )
    root_fields = {
        "schema_version",
        "task_id",
        "primary_endpoint_id",
        "primary_evaluator_id",
        "capability_registry_sha256",
        "selected",
        "baseline",
    }
    if set(value) != root_fields:
        raise SealedMetricError(
            "variant primary-capability binding has missing or unknown fields"
        )
    expected = {
        "schema_version": _VARIANT_CAPABILITY_BINDING_SCHEMA_VERSION,
        "task_id": VARIANT_TASK,
        "primary_endpoint_id": _VARIANT_PRIMARY_ENDPOINT_ID,
        "primary_evaluator_id": _VARIANT_PRIMARY_EVALUATOR_ID,
    }
    if any(
        value.get(field) != expected_value
        for field, expected_value in expected.items()
    ):
        raise SealedMetricError("variant primary-capability identity is not canonical")
    registry_sha256 = value.get("capability_registry_sha256")
    if (
        not isinstance(registry_sha256, str)
        or len(registry_sha256) != 64
        or any(character not in "0123456789abcdef" for character in registry_sha256)
    ):
        raise SealedMetricError("variant capability registry SHA-256 is invalid")
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
    for role in ("selected", "baseline"):
        role_binding = value.get(role)
        if not isinstance(role_binding, Mapping) or set(role_binding) != role_fields:
            raise SealedMetricError(
                f"variant {role} capability binding has missing or unknown fields"
            )
        model_id = role_binding.get("model_id")
        capability = role_binding.get("capability")
        if not isinstance(model_id, str) or not model_id:
            raise SealedMetricError(f"variant {role} model_id is invalid")
        if not isinstance(capability, Mapping) or set(capability) != capability_fields:
            raise SealedMetricError(
                f"variant {role} capability record has missing or unknown fields"
            )
        if capability.get("model_id") != model_id:
            raise SealedMetricError(
                f"variant {role} capability record names the wrong model"
            )
        if capability.get("role") not in {
            "baseline",
            "conditional_only",
            "primary_candidate",
        }:
            raise SealedMetricError(
                f"variant {role} capability has an invalid primary role"
            )
        for field_name in ("native_outputs", "allowed_endpoints"):
            values = _locked_string_array(capability.get(field_name))
            if (
                not values
                or values != sorted(set(values))
            ):
                raise SealedMetricError(
                    f"variant {role} capability {field_name} is invalid"
                )
        for field_name in (
            "primary_eligible",
            "requires_fitted_head",
            "requires_observed_target_context",
            "is_mandatory_baseline",
        ):
            if not isinstance(capability.get(field_name), bool):
                raise SealedMetricError(
                    f"variant {role} capability {field_name} must be boolean"
                )
        if (
            capability.get("primary_eligible") is not True
            or capability.get("requires_observed_target_context") is not False
            or _VARIANT_PRIMARY_ENDPOINT_ID
            not in capability.get("allowed_endpoints", ())
        ):
            raise SealedMetricError(
                f"variant {role} model is ineligible for signed-effect scoring"
            )
        if not set(capability.get("native_outputs", ())).intersection(
            {"gene_expression_delta", "rna_coverage_delta"}
        ):
            raise SealedMetricError(
                f"variant {role} model lacks a signed gene-level output"
            )
        capability_sha256 = role_binding.get("capability_sha256")
        if capability_sha256 != canonical_sha256(capability):
            raise SealedMetricError(
                f"variant {role} capability SHA-256 does not match its record"
            )
    return dict(value)


def _validate_variant_secondary_evaluation(
    value: Mapping[str, Any] | None,
    *,
    primary_binding: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise SealedMetricError(
            "variant sealed metrics require a locked secondary-evaluation binding"
        )
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
    if set(value) != root_fields:
        raise SealedMetricError(
            "variant secondary-evaluation binding has missing or unknown fields"
        )
    expected = {
        "schema_version": _VARIANT_SECONDARY_BINDING_SCHEMA_VERSION,
        "task_id": VARIANT_TASK,
        "endpoint_id": _VARIANT_RETRIEVAL_ENDPOINT_ID,
        "evaluator_id": _VARIANT_RETRIEVAL_EVALUATOR_ID,
        "candidate_transform_id": _ABSOLUTE_SIGNED_EFFECT_TRANSFORM_ID,
        "primary_baseline_transform_id": _ABSOLUTE_SIGNED_EFFECT_TRANSFORM_ID,
    }
    if any(value.get(field) != expected_value for field, expected_value in expected.items()):
        raise SealedMetricError(
            "variant secondary-evaluation identity or score transform is not canonical"
        )
    registry_sha256 = value.get("capability_registry_sha256")
    if registry_sha256 != primary_binding.get("capability_registry_sha256"):
        raise SealedMetricError(
            "variant primary and secondary capability registries differ"
        )
    mandatory = _locked_string_array(value.get("mandatory_model_ids"))
    if mandatory != list(_MANDATORY_RETRIEVAL_MODELS):
        raise SealedMetricError(
            "variant secondary evaluation must name the exact mandatory comparator roster"
        )
    comparators = value.get("comparators")
    if (
        isinstance(comparators, (str, bytes))
        or not isinstance(comparators, Sequence)
        or len(comparators) != len(mandatory)
    ):
        raise SealedMetricError(
            "variant secondary evaluation must bind every mandatory comparator once"
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
    expected_roles = {
        "abc": "link_only",
        "nearest_gene": "baseline",
        "re2g": "link_only",
    }
    comparator_ids: list[str] = []
    for index, comparator in enumerate(comparators):
        if not isinstance(comparator, Mapping) or set(comparator) != comparator_fields:
            raise SealedMetricError(
                f"variant secondary comparator {index} has an invalid schema"
            )
        model_id = comparator.get("model_id")
        if not isinstance(model_id, str) or not model_id:
            raise SealedMetricError(
                f"variant secondary comparator {index} has an invalid model_id"
            )
        comparator_ids.append(model_id)
        run_id = comparator.get("run_id")
        if (
            not isinstance(run_id, str)
            or len(run_id) != 64
            or any(character not in "0123456789abcdef" for character in run_id)
        ):
            raise SealedMetricError(
                f"variant secondary comparator {model_id} has an invalid run_id"
            )
        if comparator.get("score_transform_id") != _IDENTITY_LINK_SCORE_TRANSFORM_ID:
            raise SealedMetricError(
                f"variant secondary comparator {model_id} has an invalid score transform"
            )
        capability = comparator.get("capability")
        if not isinstance(capability, Mapping) or set(capability) != capability_fields:
            raise SealedMetricError(
                f"variant secondary comparator {model_id} capability is invalid"
            )
        for field_name in ("native_outputs", "allowed_endpoints"):
            identifiers = _locked_string_array(capability.get(field_name))
            if (
                not identifiers
                or identifiers != sorted(set(identifiers))
            ):
                raise SealedMetricError(
                    f"variant secondary comparator {model_id} {field_name} is invalid"
                )
        for field_name in (
            "primary_eligible",
            "requires_fitted_head",
            "requires_observed_target_context",
            "is_mandatory_baseline",
        ):
            if not isinstance(capability.get(field_name), bool):
                raise SealedMetricError(
                    f"variant secondary comparator {model_id} {field_name} must be boolean"
                )
        if (
            capability.get("model_id") != model_id
            or capability.get("role") != expected_roles.get(model_id)
            or _locked_string_array(capability.get("native_outputs"))
            != ["enhancer_gene_link_score"]
            or _locked_string_array(capability.get("allowed_endpoints"))
            != ["enhancer_gene_link", _VARIANT_RETRIEVAL_ENDPOINT_ID]
            or capability.get("primary_eligible") is not False
            or capability.get("requires_fitted_head") is not False
            or capability.get("requires_observed_target_context") is not False
            or capability.get("is_mandatory_baseline") is not True
        ):
            raise SealedMetricError(
                f"variant secondary comparator {model_id} is not an eligible link baseline"
            )
        if comparator.get("capability_sha256") != canonical_sha256(capability):
            raise SealedMetricError(
                f"variant secondary comparator {model_id} capability SHA-256 changed"
            )
    if comparator_ids != list(_MANDATORY_RETRIEVAL_MODELS):
        raise SealedMetricError(
            "variant secondary comparators are missing, duplicated, or not canonical"
        )
    primary_model_ids = {
        str(primary_binding[role]["model_id"]) for role in ("selected", "baseline")
    }
    if primary_model_ids.intersection(comparator_ids):
        raise SealedMetricError(
            "variant signed models and secondary comparators must be distinct"
        )
    return dict(value)


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise SealedMetricError(f"{label} must be numeric, not boolean")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise SealedMetricError(f"{label} must be numeric") from error
    if not isfinite(result):
        raise SealedMetricError(f"{label} must be finite")
    return result


def _index_rows(
    rows: Sequence[Mapping[str, Any]], *, label: str
) -> dict[str, Mapping[str, Any]]:
    if not rows:
        raise SealedMetricError(f"{label} contains no rows")
    result: dict[str, Mapping[str, Any]] = {}
    for row in rows:
        row_hash = row.get("row_hash")
        if not isinstance(row_hash, str) or len(row_hash) != 64:
            raise SealedMetricError(f"{label} has an invalid row_hash")
        if row_hash in result:
            raise SealedMetricError(f"{label} repeats row_hash {row_hash}")
        result[row_hash] = row
    return result


def _aligned_seed_indices(
    rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]], *, label: str
) -> tuple[tuple[int, ...], dict[int, dict[str, Mapping[str, Any]]], tuple[str, ...]]:
    seeds = tuple(sorted(rows_by_seed))
    if len(seeds) != 5 or len(set(seeds)) != 5:
        raise SealedMetricError(f"{label} must contain exactly five unique seeds")
    indices = {
        seed: _index_rows(rows_by_seed[seed], label=f"{label} seed {seed}")
        for seed in seeds
    }
    row_sets = {tuple(sorted(index)) for index in indices.values()}
    if len(row_sets) != 1:
        raise SealedMetricError(f"{label} seed row sets differ")
    return seeds, indices, next(iter(row_sets))


def _percentile(values: Sequence[float], probability: float) -> float:
    if not values:
        raise SealedMetricError("cannot take a percentile of an empty bootstrap")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = floor(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _one_sided_bootstrap_p(values: Sequence[float]) -> float:
    return (1.0 + sum(value <= 0.0 for value in values)) / (len(values) + 1.0)


def _probability_index(
    rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    *,
    class_roster: Sequence[str],
    label: str,
) -> dict[int, dict[str, dict[str, float]]]:
    seeds, indices, _ = _aligned_seed_indices(rows_by_seed, label=label)
    roster = tuple(class_roster)
    result: dict[int, dict[str, dict[str, float]]] = {}
    for seed in seeds:
        result[seed] = {}
        for row_hash, row in indices[seed].items():
            raw = row.get("probabilities")
            if not isinstance(raw, Mapping) or set(raw) != set(roster):
                raise SealedMetricError(
                    f"{label} seed {seed} probabilities must exactly cover class roster"
                )
            probabilities = {
                class_id: _finite(raw[class_id], f"{label} probability {class_id}")
                for class_id in roster
            }
            if any(value < 0.0 for value in probabilities.values()) or abs(
                sum(probabilities.values()) - 1.0
            ) > 1e-6:
                raise SealedMetricError(
                    f"{label} seed {seed} probabilities must be nonnegative and sum to one"
                )
            result[seed][row_hash] = probabilities
    return result


def _macro_f1(
    observed: Sequence[str],
    predicted: Sequence[str],
    donors: Sequence[str],
    roster: Sequence[str],
    *,
    require_full_roster: bool,
) -> tuple[float, dict[str, float]]:
    if not observed or not (len(observed) == len(predicted) == len(donors)):
        raise SealedMetricError("cell classification vectors must be aligned")
    if not set(observed).issubset(roster) or not set(predicted).issubset(roster):
        raise SealedMetricError("cell labels fall outside the frozen class roster")
    if require_full_roster and set(observed) != set(roster):
        raise SealedMetricError("cell observations do not cover the frozen class roster")
    weights = donor_class_balanced_weights(donors, observed)
    per_class: dict[str, float] = {}
    for class_id in roster:
        true_positive = sum(
            weight
            for truth, guess, weight in zip(
                observed, predicted, weights, strict=True
            )
            if truth == class_id and guess == class_id
        )
        false_positive = sum(
            weight
            for truth, guess, weight in zip(
                observed, predicted, weights, strict=True
            )
            if truth != class_id and guess == class_id
        )
        false_negative = sum(
            weight
            for truth, guess, weight in zip(
                observed, predicted, weights, strict=True
            )
            if truth == class_id and guess != class_id
        )
        denominator = 2.0 * true_positive + false_positive + false_negative
        per_class[class_id] = (
            0.0 if denominator == 0.0 else 2.0 * true_positive / denominator
        )
    return sum(per_class.values()) / len(roster), per_class


def _cell_metrics(
    *,
    selected_rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    baseline_rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    selected_probabilities_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    baseline_probabilities_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    outcome_rows: Sequence[Mapping[str, Any]],
    class_roster: Sequence[str],
    n_resamples: int,
    bootstrap_seed: int,
) -> SealedMetricResult:
    roster = tuple(class_roster)
    if not roster or tuple(sorted(set(roster))) != roster:
        raise SealedMetricError("class_roster must be non-empty, unique, and sorted")
    selected_seeds, selected, row_ids = _aligned_seed_indices(
        selected_rows_by_seed, label="selected cell predictions"
    )
    baseline_seeds, baseline, baseline_ids = _aligned_seed_indices(
        baseline_rows_by_seed, label="baseline cell predictions"
    )
    if selected_seeds != baseline_seeds or row_ids != baseline_ids:
        raise SealedMetricError("selected and baseline cell seed predictions do not align")
    selected_probabilities = _probability_index(
        selected_probabilities_by_seed,
        class_roster=roster,
        label="selected cell probabilities",
    )
    baseline_probabilities = _probability_index(
        baseline_probabilities_by_seed,
        class_roster=roster,
        label="baseline cell probabilities",
    )
    outcomes = _index_rows(outcome_rows, label="sealed cell outcomes")
    if set(outcomes) != set(row_ids):
        raise SealedMetricError("sealed cell outcome rows differ from predictions")

    observed: list[str] = []
    donors: list[str] = []
    selected_classes: list[str] = []
    baseline_classes: list[str] = []
    selected_ensemble_probabilities: list[dict[str, float]] = []
    baseline_ensemble_probabilities: list[dict[str, float]] = []
    ensemble_rows: list[dict[str, str]] = []
    for row_hash in row_ids:
        outcome = outcomes[row_hash]
        observed_class = str(outcome.get("observed_class", ""))
        donor = str(outcome.get("unit_hash", ""))
        if observed_class not in roster or len(donor) != 64:
            raise SealedMetricError("sealed cell outcome has an invalid class or donor")
        candidate_probs = {
            class_id: sum(
                selected_probabilities[seed][row_hash][class_id]
                for seed in selected_seeds
            )
            / len(selected_seeds)
            for class_id in roster
        }
        baseline_probs = {
            class_id: sum(
                baseline_probabilities[seed][row_hash][class_id]
                for seed in baseline_seeds
            )
            / len(baseline_seeds)
            for class_id in roster
        }
        candidate_class = min(roster, key=lambda item: (-candidate_probs[item], item))
        baseline_class = min(roster, key=lambda item: (-baseline_probs[item], item))
        for seed in selected_seeds:
            if selected[seed][row_hash].get("predicted_class") != min(
                roster,
                key=lambda item: (-selected_probabilities[seed][row_hash][item], item),
            ):
                raise SealedMetricError(
                    "selected hard class differs from its probability argmax"
                )
            if baseline[seed][row_hash].get("predicted_class") != min(
                roster,
                key=lambda item: (-baseline_probabilities[seed][row_hash][item], item),
            ):
                raise SealedMetricError(
                    "baseline hard class differs from its probability argmax"
                )
        observed.append(observed_class)
        donors.append(donor)
        selected_classes.append(candidate_class)
        baseline_classes.append(baseline_class)
        selected_ensemble_probabilities.append(candidate_probs)
        baseline_ensemble_probabilities.append(baseline_probs)
        ensemble_rows.append(
            {
                "row_hash": row_hash,
                "unit_hash": donor,
                "observed_class": observed_class,
                "candidate_class": candidate_class,
                "baseline_class": baseline_class,
                "candidate_probabilities_json": canonical_sha256(candidate_probs),
                "baseline_probabilities_json": canonical_sha256(baseline_probs),
            }
        )

    selected_primary, per_class = _macro_f1(
        observed, selected_classes, donors, roster, require_full_roster=True
    )
    baseline_primary, _ = _macro_f1(
        observed, baseline_classes, donors, roster, require_full_roster=True
    )
    weights = donor_class_balanced_weights(donors, observed)
    selected_brier = multiclass_brier_score(
        observed, selected_ensemble_probabilities, weights
    )
    baseline_brier = multiclass_brier_score(
        observed, baseline_ensemble_probabilities, weights
    )

    seed_deltas: dict[str, float] = {}
    for seed in selected_seeds:
        candidate = [str(selected[seed][row_id]["predicted_class"]) for row_id in row_ids]
        comparator = [str(baseline[seed][row_id]["predicted_class"]) for row_id in row_ids]
        candidate_score, _ = _macro_f1(
            observed, candidate, donors, roster, require_full_roster=True
        )
        baseline_score, _ = _macro_f1(
            observed, comparator, donors, roster, require_full_roster=True
        )
        seed_deltas[str(seed)] = candidate_score - baseline_score

    unit_roster = sorted(set(donors))
    by_unit = {
        unit: [index for index, donor in enumerate(donors) if donor == unit]
        for unit in unit_roster
    }
    rng = Random(bootstrap_seed)
    delta_bootstrap: list[float] = []
    brier_delta_bootstrap: list[float] = []
    for _ in range(n_resamples):
        sampled_indices: list[int] = []
        sampled_donors: list[str] = []
        for draw_index in range(len(unit_roster)):
            unit = unit_roster[rng.randrange(len(unit_roster))]
            sampled_indices.extend(by_unit[unit])
            sampled_donors.extend(
                [f"{draw_index}:{unit}"] * len(by_unit[unit])
            )
        sampled_observed = [observed[index] for index in sampled_indices]
        candidate_score, _ = _macro_f1(
            sampled_observed,
            [selected_classes[index] for index in sampled_indices],
            sampled_donors,
            roster,
            require_full_roster=False,
        )
        comparator_score, _ = _macro_f1(
            sampled_observed,
            [baseline_classes[index] for index in sampled_indices],
            sampled_donors,
            roster,
            require_full_roster=False,
        )
        sampled_weights = donor_class_balanced_weights(
            sampled_donors, sampled_observed
        )
        candidate_brier = multiclass_brier_score(
            sampled_observed,
            [selected_ensemble_probabilities[index] for index in sampled_indices],
            sampled_weights,
        )
        comparator_brier = multiclass_brier_score(
            sampled_observed,
            [baseline_ensemble_probabilities[index] for index in sampled_indices],
            sampled_weights,
        )
        delta_bootstrap.append(candidate_score - comparator_score)
        brier_delta_bootstrap.append(candidate_brier - comparator_brier)

    metrics = {
        "donor_class_balanced_macro_f1": selected_primary,
        "baseline_donor_class_balanced_macro_f1": baseline_primary,
        "delta_donor_class_balanced_macro_f1": selected_primary - baseline_primary,
        "delta_donor_class_balanced_macro_f1_ci_low": _percentile(
            delta_bootstrap, 0.025
        ),
        "per_class_f1": per_class,
        "multiclass_brier_score": selected_brier,
        "baseline_multiclass_brier_score": baseline_brier,
        "delta_multiclass_brier_score": selected_brier - baseline_brier,
        "raw_confirmatory_p": _one_sided_bootstrap_p(delta_bootstrap),
        "positive_seed_count": sum(value > 0.0 for value in seed_deltas.values()),
        "seed_deltas": seed_deltas,
        "seeds_are_biological_replicates": False,
        "delta_bootstrap_sha256": canonical_sha256(delta_bootstrap),
        "brier_delta_bootstrap_sha256": canonical_sha256(brier_delta_bootstrap),
    }
    return SealedMetricResult(
        task_id=CELL_TASK,
        evaluator_id=SEALED_METRIC_EVALUATOR_ID,
        metrics=metrics,
        ensemble_rows=tuple(ensemble_rows),
        n_resamples=n_resamples,
        bootstrap_seed=bootstrap_seed,
        independent_unit="donor",
        independent_unit_count=len(unit_roster),
        row_count=len(row_ids),
    )


def _variant_primary(
    rows: Sequence[Mapping[str, Any]], strata: Sequence[str]
) -> tuple[float, dict[str, float], dict[str, float]]:
    candidate_by_stratum: dict[str, float] = {}
    baseline_by_stratum: dict[str, float] = {}
    for stratum in strata:
        subset = [row for row in rows if row["stratum"] == stratum]
        if len(subset) < 2:
            raise SealedMetricError(f"variant stratum {stratum} has fewer than two rows")
        observed = [_finite(row["observed"], "observed variant effect") for row in subset]
        candidate_by_stratum[stratum] = spearman_correlation(
            observed,
            [_finite(row["candidate"], "candidate variant score") for row in subset],
        )
        baseline_by_stratum[stratum] = spearman_correlation(
            observed,
            [_finite(row["baseline"], "baseline variant score") for row in subset],
        )
    return (
        fisher_z_mean([candidate_by_stratum[item] for item in strata]),
        candidate_by_stratum,
        baseline_by_stratum,
    )


def _conservative_average_precision(labels: Sequence[int], scores: Sequence[float]) -> float:
    if sum(labels) == 0:
        return 0.0
    return average_precision(labels, scores)


def _variant_metrics(
    *,
    selected_rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    baseline_rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    secondary_rows_by_model: Mapping[str, Sequence[Mapping[str, Any]]],
    outcome_rows: Sequence[Mapping[str, Any]],
    secondary_outcome_rows: Sequence[Mapping[str, Any]],
    primary_binding: Mapping[str, Any],
    secondary_binding: Mapping[str, Any],
    strata: Sequence[str],
    n_resamples: int,
    bootstrap_seed: int,
) -> SealedMetricResult:
    roster = tuple(strata)
    if not roster or tuple(sorted(set(roster))) != roster:
        raise SealedMetricError("variant strata must be non-empty, unique, and sorted")
    selected_seeds, selected, row_ids = _aligned_seed_indices(
        selected_rows_by_seed, label="selected variant predictions"
    )
    baseline_seeds, baseline, baseline_ids = _aligned_seed_indices(
        baseline_rows_by_seed, label="baseline variant predictions"
    )
    if selected_seeds != baseline_seeds or row_ids != baseline_ids:
        raise SealedMetricError("selected and baseline variant predictions do not align")
    comparator_models = tuple(
        str(comparator["model_id"])
        for comparator in secondary_binding["comparators"]
    )
    if set(secondary_rows_by_model) != set(comparator_models):
        raise SealedMetricError(
            "variant secondary prediction tables do not exactly cover the locked comparators"
        )
    secondary_indices: dict[str, dict[str, Mapping[str, Any]]] = {}
    for model_id in comparator_models:
        index = _index_rows(
            secondary_rows_by_model[model_id],
            label=f"secondary variant predictions for {model_id}",
        )
        if tuple(sorted(index)) != row_ids:
            raise SealedMetricError(
                f"secondary variant predictions for {model_id} do not align"
            )
        secondary_indices[model_id] = index
    outcomes = _index_rows(outcome_rows, label="sealed variant outcomes")
    secondary = _index_rows(
        secondary_outcome_rows, label="sealed variant secondary outcomes"
    )
    if set(outcomes) != set(row_ids) or set(secondary) != set(row_ids):
        raise SealedMetricError("sealed variant outcome rows differ from predictions")

    ensemble_rows: list[dict[str, str]] = []
    for row_hash in row_ids:
        outcome = outcomes[row_hash]
        extra = secondary[row_hash]
        metadata = {
            field: str(outcome.get(field, ""))
            for field in ("unit_hash", "block_hash", "stratum")
        }
        if any(str(extra.get(field, "")) != value for field, value in metadata.items()):
            raise SealedMetricError("variant primary and secondary outcome metadata differ")
        for seed in selected_seeds:
            if any(
                str(selected[seed][row_hash].get(field, "")) != value
                or str(baseline[seed][row_hash].get(field, "")) != value
                for field, value in metadata.items()
            ):
                raise SealedMetricError(
                    "variant prediction metadata differs from the sealed outcome"
                )
        for model_id in comparator_models:
            if any(
                str(secondary_indices[model_id][row_hash].get(field, "")) != value
                for field, value in metadata.items()
            ):
                raise SealedMetricError(
                    f"secondary variant metadata for {model_id} differs from the sealed outcome"
                )
        if metadata["stratum"] not in roster:
            raise SealedMetricError("variant outcome stratum is outside the frozen roster")
        if str(extra.get("observed_binary", "")) not in {"0", "1"}:
            raise SealedMetricError("variant observed_binary must be 0 or 1")
        if str(extra.get("summary_statistics_complete", "")).lower() != "true":
            complete = "false"
        else:
            complete = "true"
        candidate = sum(
            _finite(selected[seed][row_hash].get("predicted"), "selected prediction")
            for seed in selected_seeds
        ) / len(selected_seeds)
        comparator = sum(
            _finite(baseline[seed][row_hash].get("predicted"), "baseline prediction")
            for seed in baseline_seeds
        ) / len(baseline_seeds)
        retrieval_scores = {
            "candidate": abs(candidate),
            "primary_baseline": abs(comparator),
            **{
                model_id: _finite(
                    secondary_indices[model_id][row_hash].get("predicted"),
                    f"secondary prediction for {model_id}",
                )
                for model_id in comparator_models
            },
        }
        ensemble_rows.append(
            {
                "row_hash": row_hash,
                **metadata,
                "observed": str(_finite(outcome.get("observed"), "observed effect")),
                "observed_binary": str(extra["observed_binary"]),
                "summary_statistics_complete": complete,
                "candidate": str(candidate),
                "baseline": str(comparator),
                "retrieval_scores_json": json.dumps(
                    retrieval_scores, sort_keys=True, separators=(",", ":")
                ),
            }
        )

    candidate_primary, candidate_strata, baseline_strata = _variant_primary(
        ensemble_rows, roster
    )
    baseline_primary = fisher_z_mean([baseline_strata[item] for item in roster])
    stratum_deltas = {
        item: candidate_strata[item] - baseline_strata[item] for item in roster
    }
    selected_model_id = str(primary_binding["selected"]["model_id"])
    baseline_model_id = str(primary_binding["baseline"]["model_id"])
    retrieval_model_ids = (selected_model_id, baseline_model_id, *comparator_models)
    if len(set(retrieval_model_ids)) != len(retrieval_model_ids):
        raise SealedMetricError("variant retrieval model identifiers must be distinct")
    retrieval_scores_by_model: dict[str, dict[str, float]] = {
        model_id: {} for model_id in retrieval_model_ids
    }
    for row in ensemble_rows:
        row_hash = row["row_hash"]
        scores = json.loads(row["retrieval_scores_json"])
        retrieval_scores_by_model[selected_model_id][row_hash] = _finite(
            scores["candidate"], "candidate retrieval score"
        )
        retrieval_scores_by_model[baseline_model_id][row_hash] = _finite(
            scores["primary_baseline"], "primary-baseline retrieval score"
        )
        for model_id in comparator_models:
            retrieval_scores_by_model[model_id][row_hash] = _finite(
                scores[model_id], f"retrieval score for {model_id}"
            )
    labels = [int(row["observed_binary"]) for row in ensemble_rows]
    retrieval_auprc_by_model = {
        model_id: _conservative_average_precision(
            labels,
            [
                retrieval_scores_by_model[model_id][row["row_hash"]]
                for row in ensemble_rows
            ],
        )
        for model_id in retrieval_model_ids
    }
    candidate_auprc = retrieval_auprc_by_model[selected_model_id]
    primary_baseline_auprc = retrieval_auprc_by_model[baseline_model_id]
    retrieval_comparator_ids = (baseline_model_id, *comparator_models)
    best_retrieval_comparator_id = min(
        retrieval_comparator_ids,
        key=lambda model_id: (-retrieval_auprc_by_model[model_id], model_id),
    )
    best_retrieval_comparator_auprc = retrieval_auprc_by_model[
        best_retrieval_comparator_id
    ]
    seed_deltas: dict[str, float] = {}
    for seed in selected_seeds:
        seed_rows = [
            {
                **row,
                "candidate": str(selected[seed][row["row_hash"]]["predicted"]),
                "baseline": str(baseline[seed][row["row_hash"]]["predicted"]),
            }
            for row in ensemble_rows
        ]
        seed_candidate, _, seed_baseline_strata = _variant_primary(seed_rows, roster)
        seed_baseline = fisher_z_mean(
            [seed_baseline_strata[item] for item in roster]
        )
        seed_deltas[str(seed)] = seed_candidate - seed_baseline

    blocks = sorted({row["block_hash"] for row in ensemble_rows})
    if len(blocks) < 2:
        raise SealedMetricError("variant sealed metric requires at least two LD blocks")
    by_block = {
        block: [row for row in ensemble_rows if row["block_hash"] == block]
        for block in blocks
    }
    rng = Random(bootstrap_seed)
    primary_delta_bootstrap: list[float] = []
    retrieval_delta_bootstrap: list[float] = []
    for _ in range(n_resamples):
        sampled = [
            row
            for _draw in range(len(blocks))
            for row in by_block[blocks[rng.randrange(len(blocks))]]
        ]
        sampled_candidate, _, sampled_baseline_strata = _variant_primary(
            sampled, roster
        )
        sampled_baseline = fisher_z_mean(
            [sampled_baseline_strata[item] for item in roster]
        )
        primary_delta_bootstrap.append(sampled_candidate - sampled_baseline)
        sampled_labels = [int(row["observed_binary"]) for row in sampled]
        sampled_auprcs = {
            model_id: _conservative_average_precision(
                sampled_labels,
                [
                    retrieval_scores_by_model[model_id][row["row_hash"]]
                    for row in sampled
                ],
            )
            for model_id in retrieval_model_ids
        }
        retrieval_delta_bootstrap.append(
            sampled_auprcs[selected_model_id]
            - max(sampled_auprcs[model_id] for model_id in retrieval_comparator_ids)
        )

    metrics = {
        "complete_summary_statistics_available": all(
            row["summary_statistics_complete"] == "true" for row in ensemble_rows
        ),
        "mean_fisher_z_celltype_spearman": candidate_primary,
        "baseline_mean_fisher_z_celltype_spearman": baseline_primary,
        "delta_mean_fisher_z_celltype_spearman": candidate_primary - baseline_primary,
        "delta_mean_fisher_z_celltype_spearman_ci_low": _percentile(
            primary_delta_bootstrap, 0.025
        ),
        "major_lineage_signed_effect_spearman_deltas": stratum_deltas,
        "eqtl_retrieval_auprc": candidate_auprc,
        "primary_baseline_eqtl_retrieval_auprc": primary_baseline_auprc,
        "retrieval_auprc_by_model": retrieval_auprc_by_model,
        "retrieval_score_transform_by_model": {
            selected_model_id: _ABSOLUTE_SIGNED_EFFECT_TRANSFORM_ID,
            baseline_model_id: _ABSOLUTE_SIGNED_EFFECT_TRANSFORM_ID,
            **{
                model_id: _IDENTITY_LINK_SCORE_TRANSFORM_ID
                for model_id in comparator_models
            },
        },
        "best_eqtl_retrieval_comparator_model_id": best_retrieval_comparator_id,
        "best_eqtl_retrieval_comparator_auprc": best_retrieval_comparator_auprc,
        "delta_vs_best_retrieval_comparator_auprc": (
            candidate_auprc - best_retrieval_comparator_auprc
        ),
        "delta_vs_best_retrieval_comparator_auprc_ci_low": _percentile(
            retrieval_delta_bootstrap, 0.025
        ),
        "retrieval_endpoint_in_confirmatory_holm_family": False,
        "raw_confirmatory_p": _one_sided_bootstrap_p(primary_delta_bootstrap),
        "positive_seed_count": sum(value > 0.0 for value in seed_deltas.values()),
        "seed_deltas": seed_deltas,
        "seeds_are_biological_replicates": False,
        "delta_bootstrap_sha256": canonical_sha256(primary_delta_bootstrap),
        "retrieval_delta_bootstrap_sha256": canonical_sha256(
            retrieval_delta_bootstrap
        ),
    }
    return SealedMetricResult(
        task_id=VARIANT_TASK,
        evaluator_id=SEALED_METRIC_EVALUATOR_ID,
        metrics=metrics,
        ensemble_rows=tuple(ensemble_rows),
        n_resamples=n_resamples,
        bootstrap_seed=bootstrap_seed,
        independent_unit="ld_block",
        independent_unit_count=len(blocks),
        row_count=len(row_ids),
    )


def recompute_sealed_metrics(
    *,
    task_id: str,
    selected_rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    baseline_rows_by_seed: Mapping[int, Sequence[Mapping[str, Any]]],
    outcome_rows: Sequence[Mapping[str, Any]],
    endpoint_parameters: Mapping[str, Any],
    n_resamples: int,
    bootstrap_seed: int,
    selected_probabilities_by_seed: Mapping[
        int, Sequence[Mapping[str, Any]]
    ] | None = None,
    baseline_probabilities_by_seed: Mapping[
        int, Sequence[Mapping[str, Any]]
    ] | None = None,
    secondary_outcome_rows: Sequence[Mapping[str, Any]] | None = None,
    variant_primary_capability: Mapping[str, Any] | None = None,
    variant_secondary_evaluation: Mapping[str, Any] | None = None,
    variant_secondary_rows_by_model: Mapping[
        str, Sequence[Mapping[str, Any]]
    ] | None = None,
) -> SealedMetricResult:
    """Dispatch one exact externally opened task to its registered evaluator."""

    if task_id == VARIANT_TASK:
        primary_binding = _validate_variant_primary_capability(
            variant_primary_capability
        )
        secondary_binding = _validate_variant_secondary_evaluation(
            variant_secondary_evaluation,
            primary_binding=primary_binding,
        )
    elif (
        variant_primary_capability is not None
        or variant_secondary_evaluation is not None
        or variant_secondary_rows_by_model is not None
    ):
        raise SealedMetricError(
            "variant capability, secondary evaluation, or rows are invalid for a nonvariant task"
        )
    if n_resamples != 10_000:
        raise SealedMetricError("sealed metrics require exactly 10,000 resamples")
    if not isinstance(bootstrap_seed, int) or isinstance(bootstrap_seed, bool):
        raise SealedMetricError("bootstrap_seed must be an integer")
    if task_id == CELL_TASK:
        if set(endpoint_parameters) != {"class_roster"}:
            raise SealedMetricError("cell endpoint parameters require class_roster only")
        if selected_probabilities_by_seed is None or baseline_probabilities_by_seed is None:
            raise SealedMetricError(
                "cell sealed metrics require selected and baseline probability artifacts"
            )
        return _cell_metrics(
            selected_rows_by_seed=selected_rows_by_seed,
            baseline_rows_by_seed=baseline_rows_by_seed,
            selected_probabilities_by_seed=selected_probabilities_by_seed,
            baseline_probabilities_by_seed=baseline_probabilities_by_seed,
            outcome_rows=outcome_rows,
            class_roster=endpoint_parameters["class_roster"],
            n_resamples=n_resamples,
            bootstrap_seed=bootstrap_seed,
        )
    if task_id == VARIANT_TASK:
        if set(endpoint_parameters) != {"strata"}:
            raise SealedMetricError("variant endpoint parameters require strata only")
        if secondary_outcome_rows is None:
            raise SealedMetricError(
                "variant sealed metrics require a secondary eQTL-status outcome artifact"
            )
        if variant_secondary_rows_by_model is None:
            raise SealedMetricError(
                "variant sealed metrics require locked secondary comparator predictions"
            )
        return _variant_metrics(
            selected_rows_by_seed=selected_rows_by_seed,
            baseline_rows_by_seed=baseline_rows_by_seed,
            secondary_rows_by_model=variant_secondary_rows_by_model,
            outcome_rows=outcome_rows,
            secondary_outcome_rows=secondary_outcome_rows,
            primary_binding=primary_binding,
            secondary_binding=secondary_binding,
            strata=endpoint_parameters["strata"],
            n_resamples=n_resamples,
            bootstrap_seed=bootstrap_seed,
        )
    raise SealedMetricError(
        f"task {task_id!r} has no implemented positive sealed-metric evaluator"
    )


__all__ = [
    "SEALED_METRIC_EVALUATOR_ID",
    "SealedMetricError",
    "SealedMetricResult",
    "recompute_sealed_metrics",
]
