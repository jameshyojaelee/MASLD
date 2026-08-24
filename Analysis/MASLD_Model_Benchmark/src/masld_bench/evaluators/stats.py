"""Small, deterministic statistical helpers for the MASLD benchmark.

The benchmark deliberately keeps its inferential unit explicit.  In particular,
``paired_cluster_bootstrap`` first reduces observations within a biological
cluster and then resamples clusters with equal weight.  Cells, spots, guides,
or sequencing runs therefore cannot silently inflate the sample size.
"""

from __future__ import annotations

from collections.abc import Callable, Hashable, Mapping, Sequence
from dataclasses import asdict, dataclass
from math import floor, isfinite, sqrt
from random import Random
from statistics import NormalDist, fmean
from typing import Any, TypeAlias


PValueInput: TypeAlias = Sequence[float] | Mapping[Hashable, float]
PValueOutput: TypeAlias = list[float] | dict[Hashable, float]
DEFAULT_BOOTSTRAP_SEED = 17
DEFAULT_TRAINING_SEEDS = (1103, 2909, 4721, 6673, 8111)


@dataclass(frozen=True, slots=True)
class BootstrapInterval:
    """Paired cluster-bootstrap estimate and percentile interval."""

    estimate: float
    lower: float
    upper: float
    confidence_level: float
    n_resamples: int
    n_clusters: int
    n_observations: int
    seed: int
    bootstrap_p_value: float
    probability_improvement: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class TwoWayBootstrapInterval:
    """Paired interval from independent donor and genomic-block resampling."""

    estimate: float
    lower: float
    upper: float
    confidence_level: float
    n_resamples: int
    n_donors: int
    n_genomic_blocks: int
    n_cells: int
    n_observations: int
    seed: int
    bootstrap_p_value: float
    probability_improvement: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class SeedDirectionStabilityResult:
    """Five-seed direction diagnostic, never an inferential sample size."""

    passed: bool
    qualifying_seed_count: int
    evaluated_seed_count: int
    required_seed_count: int
    seed_ids: tuple[int, ...]
    oriented_differences: tuple[float, ...]
    seeds_are_biological_replicates: bool
    supports_inferential_test: bool
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PowerGateResult:
    """Prospective power check expressed through a minimum detectable effect."""

    passed: bool
    n_units: int
    paired_sd: float
    minimum_effect: float
    minimum_detectable_effect: float
    alpha: float
    effective_alpha: float
    power: float
    n_primary_claims: int
    min_units: int
    two_sided: bool
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class EmpiricalPowerGateResult:
    """Power derived from a frozen endpoint-recomputed paired bootstrap.

    The bootstrap distribution must contain the task's actual nonlinear
    endpoint recomputed from resampled biological units.  This class therefore
    supports macro-F1, AUPRC, rank correlation, and deviance without pretending
    that their sampling errors are normal or that cells are replicates.
    """

    method_id: str
    passed: bool
    n_units: int
    min_units: int
    observed_effect: float
    minimum_effect: float
    empirical_standard_error: float
    critical_value: float
    estimated_power: float
    power_target: float
    alpha: float
    effective_alpha: float
    n_primary_claims: int
    two_sided: bool
    n_resamples: int
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _validated_p_values(p_values: PValueInput) -> tuple[list[Hashable], list[float], bool]:
    is_mapping = isinstance(p_values, Mapping)
    if is_mapping:
        keys = list(p_values)
        values = [p_values[key] for key in keys]
    else:
        if isinstance(p_values, (str, bytes)):
            raise TypeError("p_values must be a sequence of numbers, not text")
        keys = list(range(len(p_values)))
        values = list(p_values)

    checked: list[float] = []
    for value in values:
        if isinstance(value, bool):
            raise TypeError("p-values must be numeric probabilities, not booleans")
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"invalid p-value: {value!r}") from exc
        if not isfinite(number) or not 0.0 <= number <= 1.0:
            raise ValueError(f"p-value must be finite and in [0, 1], got {number!r}")
        checked.append(number)
    return keys, checked, is_mapping


def _restore_p_value_shape(
    keys: Sequence[Hashable], values: Sequence[float], is_mapping: bool
) -> PValueOutput:
    if is_mapping:
        return {key: value for key, value in zip(keys, values, strict=True)}
    return list(values)


def holm_correction(p_values: PValueInput) -> PValueOutput:
    """Return Holm family-wise-error adjusted p-values.

    Mapping inputs retain their insertion order and keys.  Sequence inputs
    return a list in the original order.
    """

    keys, values, is_mapping = _validated_p_values(p_values)
    count = len(values)
    if count == 0:
        return _restore_p_value_shape(keys, [], is_mapping)

    order = sorted(range(count), key=lambda index: (values[index], index))
    adjusted = [0.0] * count
    running_max = 0.0
    for rank, index in enumerate(order):
        candidate = min(1.0, (count - rank) * values[index])
        running_max = max(running_max, candidate)
        adjusted[index] = running_max
    return _restore_p_value_shape(keys, adjusted, is_mapping)


def benjamini_hochberg(p_values: PValueInput) -> PValueOutput:
    """Return Benjamini-Hochberg false-discovery-rate adjusted p-values."""

    keys, values, is_mapping = _validated_p_values(p_values)
    count = len(values)
    if count == 0:
        return _restore_p_value_shape(keys, [], is_mapping)

    order = sorted(range(count), key=lambda index: (values[index], index))
    adjusted = [0.0] * count
    running_min = 1.0
    for reverse_rank in range(count - 1, -1, -1):
        index = order[reverse_rank]
        rank = reverse_rank + 1
        candidate = min(1.0, values[index] * count / rank)
        running_min = min(running_min, candidate)
        adjusted[index] = running_min
    return _restore_p_value_shape(keys, adjusted, is_mapping)


# Short aliases are convenient in result-table code and retain explicit names
# above for public documentation.
holm_adjust = holm_correction
bh_adjust = benjamini_hochberg


def _coerce_finite_series(values: Sequence[float], name: str) -> list[float]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be a numeric sequence")
    result: list[float] = []
    for value in values:
        if isinstance(value, bool):
            raise TypeError(f"{name} cannot contain booleans")
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{name} contains a non-numeric value: {value!r}") from exc
        if not isfinite(number):
            raise ValueError(f"{name} must contain only finite values")
        result.append(number)
    return result


def _validated_unit_ids(
    values: Sequence[Hashable], name: str
) -> list[Hashable]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be a sequence of independent-unit identifiers")
    result = list(values)
    for value in result:
        if value is None or isinstance(value, bool):
            raise TypeError(f"{name} must contain non-null, non-boolean identifiers")
        if isinstance(value, str) and not value.strip():
            raise ValueError(f"{name} cannot contain an empty identifier")
        try:
            hash(value)
        except TypeError as exc:
            raise TypeError(f"{name} identifiers must be hashable") from exc
    return result


def _validate_independent_unit_count(
    observed: int,
    *,
    name: str,
    minimum: int,
    expected: int | None,
) -> None:
    if not isinstance(minimum, int) or isinstance(minimum, bool) or minimum < 2:
        raise ValueError(f"minimum {name} count must be an integer of at least two")
    if observed < minimum:
        raise ValueError(
            f"observed {observed} independent {name} is below the required {minimum}"
        )
    if expected is None:
        return
    if not isinstance(expected, int) or isinstance(expected, bool) or expected < 2:
        raise ValueError(f"expected {name} count must be an integer of at least two")
    if observed != expected:
        raise ValueError(
            f"observed {observed} independent {name} does not match expected {expected}"
        )


def _validate_bootstrap_controls(
    n_resamples: int, confidence_level: float, seed: int
) -> tuple[float, int]:
    if not isinstance(n_resamples, int) or isinstance(n_resamples, bool) or n_resamples < 1:
        raise ValueError("n_resamples must be a positive integer")
    if isinstance(confidence_level, bool):
        raise TypeError("confidence_level must be numeric, not boolean")
    try:
        checked_confidence = float(confidence_level)
    except (TypeError, ValueError) as exc:
        raise TypeError("confidence_level must be numeric") from exc
    if not isfinite(checked_confidence) or not 0.0 < checked_confidence < 1.0:
        raise ValueError("confidence_level must be finite and between zero and one")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise TypeError("seed must be an integer")
    return checked_confidence, seed


def _percentile(sorted_values: Sequence[float], probability: float) -> float:
    if not sorted_values:
        raise ValueError("cannot take a percentile of an empty sequence")
    position = (len(sorted_values) - 1) * probability
    lower_index = floor(position)
    upper_index = min(lower_index + 1, len(sorted_values) - 1)
    fraction = position - lower_index
    return (
        sorted_values[lower_index] * (1.0 - fraction)
        + sorted_values[upper_index] * fraction
    )


def empirical_bootstrap_power_gate(
    *,
    paired_difference_bootstrap: Sequence[float],
    observed_effect: float,
    n_units: int,
    minimum_effect: float,
    alpha: float = 0.05,
    power_target: float = 0.80,
    n_primary_claims: int = 1,
    min_units: int = 10,
    two_sided: bool = True,
) -> EmpiricalPowerGateResult:
    """Estimate prospective detection power from endpoint-recomputed errors.

    ``paired_difference_bootstrap`` is not a vector of row-level model
    differences.  It must be the frozen distribution obtained by resampling the
    registered independent units and recomputing the complete paired endpoint
    for candidate and baseline on every resample.  The function centers that
    distribution at ``observed_effect`` to obtain empirical sampling errors,
    then evaluates the prespecified ``minimum_effect`` against the first-step
    Holm threshold.  The production evidence contract separately requires
    10,000 resamples and binds the source predictions, outcomes, evaluator, and
    unit set; this pure function also permits smaller deterministic test
    fixtures.
    """

    distribution = _coerce_finite_series(
        paired_difference_bootstrap, "paired_difference_bootstrap"
    )
    if len(distribution) < 100:
        raise ValueError(
            "paired_difference_bootstrap requires at least 100 endpoint-recomputed resamples"
        )
    if not isinstance(n_units, int) or isinstance(n_units, bool) or n_units < 2:
        raise ValueError("n_units must be an integer of at least two")
    if not isinstance(min_units, int) or isinstance(min_units, bool) or min_units < 2:
        raise ValueError("min_units must be an integer of at least two")
    if (
        not isinstance(n_primary_claims, int)
        or isinstance(n_primary_claims, bool)
        or n_primary_claims < 1
    ):
        raise ValueError("n_primary_claims must be a positive integer")
    checked: dict[str, float] = {}
    for name, value in (
        ("observed_effect", observed_effect),
        ("minimum_effect", minimum_effect),
        ("alpha", alpha),
        ("power_target", power_target),
    ):
        if isinstance(value, bool):
            raise TypeError(f"{name} must be numeric, not boolean")
        try:
            number = float(value)
        except (TypeError, ValueError) as error:
            raise TypeError(f"{name} must be numeric") from error
        if not isfinite(number):
            raise ValueError(f"{name} must be finite")
        checked[name] = number
    observed = checked["observed_effect"]
    minimum = checked["minimum_effect"]
    alpha_value = checked["alpha"]
    target = checked["power_target"]
    if minimum <= 0.0:
        raise ValueError("minimum_effect must be positive")
    if not 0.0 < alpha_value < 1.0:
        raise ValueError("alpha must be between zero and one")
    if not 0.5 < target < 1.0:
        raise ValueError("power_target must be greater than 0.5 and less than one")
    if not isinstance(two_sided, bool):
        raise TypeError("two_sided must be boolean")

    effective_alpha = alpha_value / n_primary_claims
    errors = [value - observed for value in distribution]
    empirical_standard_error = sqrt(fmean(value * value for value in errors))
    if two_sided:
        null_statistics = sorted(abs(value) for value in errors)
        critical_value = _percentile(null_statistics, 1.0 - effective_alpha)
        detected = sum(
            abs(minimum + error) > critical_value for error in errors
        )
    else:
        null_statistics = sorted(errors)
        critical_value = _percentile(null_statistics, 1.0 - effective_alpha)
        detected = sum(minimum + error > critical_value for error in errors)
    estimated_power = (detected + 1.0) / (len(errors) + 1.0)
    enough_units = n_units >= min_units
    enough_power = estimated_power >= target
    passed = enough_units and enough_power
    if not enough_units:
        reason = f"{n_units} independent units is below the locked floor of {min_units}"
    elif not enough_power:
        reason = (
            f"empirical power {estimated_power:.6g} is below the locked target "
            f"{target:.6g} at minimum effect {minimum:.6g}"
        )
    else:
        reason = (
            f"empirical power {estimated_power:.6g} meets the locked target "
            f"{target:.6g} at minimum effect {minimum:.6g}"
        )
    return EmpiricalPowerGateResult(
        method_id="empirical_paired_endpoint_bootstrap_v1",
        passed=passed,
        n_units=n_units,
        min_units=min_units,
        observed_effect=observed,
        minimum_effect=minimum,
        empirical_standard_error=empirical_standard_error,
        critical_value=critical_value,
        estimated_power=estimated_power,
        power_target=target,
        alpha=alpha_value,
        effective_alpha=effective_alpha,
        n_primary_claims=n_primary_claims,
        two_sided=two_sided,
        n_resamples=len(distribution),
        reason=reason,
    )


def paired_cluster_bootstrap(
    candidate: Sequence[float],
    baseline: Sequence[float],
    clusters: Sequence[Hashable],
    *,
    n_resamples: int = 10_000,
    confidence_level: float = 0.95,
    seed: int = 17,
    higher_is_better: bool = True,
    within_cluster_statistic: Callable[[Sequence[float]], float] = fmean,
    between_cluster_statistic: Callable[[Sequence[float]], float] = fmean,
) -> BootstrapInterval:
    """Bootstrap a paired model difference at the independent-cluster level.

    The candidate and baseline must have aligned observations.  Each cluster is
    reduced independently, producing one paired difference per biological unit.
    Clusters are then sampled with replacement and receive equal weight.  For a
    non-additive metric such as macro-F1, callers should pass one precomputed
    metric per donor and a unique cluster ID for each donor.

    Positive estimates always favor ``candidate``.  Set ``higher_is_better`` to
    false for losses such as deviance or Brier score.
    """

    candidate_values = _coerce_finite_series(candidate, "candidate")
    baseline_values = _coerce_finite_series(baseline, "baseline")
    cluster_values = list(clusters)
    if not candidate_values:
        raise ValueError("at least one paired observation is required")
    if not (
        len(candidate_values) == len(baseline_values) == len(cluster_values)
    ):
        raise ValueError("candidate, baseline, and clusters must have equal lengths")
    if not isinstance(n_resamples, int) or isinstance(n_resamples, bool) or n_resamples < 1:
        raise ValueError("n_resamples must be a positive integer")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between zero and one")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise TypeError("seed must be an integer")

    grouped_indices: dict[Hashable, list[int]] = {}
    for index, cluster in enumerate(cluster_values):
        try:
            grouped_indices.setdefault(cluster, []).append(index)
        except TypeError as exc:
            raise TypeError("cluster identifiers must be hashable") from exc
    if len(grouped_indices) < 2:
        raise ValueError("paired cluster bootstrap requires at least two clusters")

    direction = 1.0 if higher_is_better else -1.0
    cluster_differences: list[float] = []
    for indices in grouped_indices.values():
        candidate_cluster = [candidate_values[index] for index in indices]
        baseline_cluster = [baseline_values[index] for index in indices]
        candidate_summary = float(within_cluster_statistic(candidate_cluster))
        baseline_summary = float(within_cluster_statistic(baseline_cluster))
        if not isfinite(candidate_summary) or not isfinite(baseline_summary):
            raise ValueError("within_cluster_statistic returned a non-finite value")
        cluster_differences.append(direction * (candidate_summary - baseline_summary))

    estimate = float(between_cluster_statistic(cluster_differences))
    if not isfinite(estimate):
        raise ValueError("between_cluster_statistic returned a non-finite value")

    rng = Random(seed)
    cluster_count = len(cluster_differences)
    bootstrap_estimates: list[float] = []
    null_differences = [value - estimate for value in cluster_differences]
    null_estimates: list[float] = []
    for _ in range(n_resamples):
        indices = [rng.randrange(cluster_count) for _ in range(cluster_count)]
        sample = [cluster_differences[index] for index in indices]
        bootstrap_estimate = float(between_cluster_statistic(sample))
        if not isfinite(bootstrap_estimate):
            raise ValueError("between_cluster_statistic returned a non-finite bootstrap value")
        bootstrap_estimates.append(bootstrap_estimate)
        null_estimate = float(
            between_cluster_statistic([null_differences[index] for index in indices])
        )
        if not isfinite(null_estimate):
            raise ValueError("between_cluster_statistic returned a non-finite null value")
        null_estimates.append(null_estimate)

    bootstrap_estimates.sort()
    tail = (1.0 - confidence_level) / 2.0
    lower = _percentile(bootstrap_estimates, tail)
    upper = _percentile(bootstrap_estimates, 1.0 - tail)
    bootstrap_p_value = (
        sum(abs(value) >= abs(estimate) for value in null_estimates) + 1
    ) / (n_resamples + 1)
    probability_improvement = sum(
        value > 0.0 for value in bootstrap_estimates
    ) / n_resamples

    return BootstrapInterval(
        estimate=estimate,
        lower=lower,
        upper=upper,
        confidence_level=confidence_level,
        n_resamples=n_resamples,
        n_clusters=cluster_count,
        n_observations=len(candidate_values),
        seed=seed,
        bootstrap_p_value=bootstrap_p_value,
        probability_improvement=probability_improvement,
    )


def two_way_donor_block_bootstrap(
    candidate: Sequence[float],
    baseline: Sequence[float],
    donor_ids: Sequence[Hashable],
    genomic_block_ids: Sequence[Hashable],
    *,
    n_resamples: int = 10_000,
    confidence_level: float = 0.95,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
    higher_is_better: bool = True,
    minimum_donors: int = 2,
    minimum_genomic_blocks: int = 2,
    expected_n_donors: int | None = None,
    expected_n_genomic_blocks: int | None = None,
    within_cell_statistic: Callable[[Sequence[float]], float] = fmean,
    between_cell_statistic: Callable[[Sequence[float]], float] = fmean,
) -> TwoWayBootstrapInterval:
    """Bootstrap a paired difference across crossed donors and genomic blocks.

    Observations are first reduced within each donor-by-block cell. Donors and
    genomic blocks are then sampled independently with replacement, and every
    sampled donor-by-block cell receives equal weight. A complete crossed grid
    is required so missing cells cannot silently change the weighting or the
    independent-unit counts.

    Donors and genomic blocks are the independent units. Repeated observations
    within a cell only improve that cell's summary and never increase either
    inferential count. Positive estimates always favor ``candidate``.
    """

    candidate_values = _coerce_finite_series(candidate, "candidate")
    baseline_values = _coerce_finite_series(baseline, "baseline")
    donors = _validated_unit_ids(donor_ids, "donor_ids")
    blocks = _validated_unit_ids(genomic_block_ids, "genomic_block_ids")
    if not candidate_values:
        raise ValueError("at least one paired observation is required")
    if not (
        len(candidate_values)
        == len(baseline_values)
        == len(donors)
        == len(blocks)
    ):
        raise ValueError(
            "candidate, baseline, donor_ids, and genomic_block_ids must have equal lengths"
        )
    confidence_level, seed = _validate_bootstrap_controls(
        n_resamples, confidence_level, seed
    )
    if not isinstance(higher_is_better, bool):
        raise TypeError("higher_is_better must be boolean")
    if not callable(within_cell_statistic) or not callable(between_cell_statistic):
        raise TypeError("bootstrap statistics must be callable")

    donor_order = list(dict.fromkeys(donors))
    block_order = list(dict.fromkeys(blocks))
    _validate_independent_unit_count(
        len(donor_order),
        name="donors",
        minimum=minimum_donors,
        expected=expected_n_donors,
    )
    _validate_independent_unit_count(
        len(block_order),
        name="genomic blocks",
        minimum=minimum_genomic_blocks,
        expected=expected_n_genomic_blocks,
    )

    cell_indices: dict[tuple[Hashable, Hashable], list[int]] = {}
    for index, (donor, block) in enumerate(zip(donors, blocks, strict=True)):
        cell_indices.setdefault((donor, block), []).append(index)
    missing_cell_count = sum(
        (donor, block) not in cell_indices
        for donor in donor_order
        for block in block_order
    )
    if missing_cell_count:
        raise ValueError(
            "two-way bootstrap requires a complete donor-by-genomic-block grid; "
            f"{missing_cell_count} cells are missing"
        )

    direction = 1.0 if higher_is_better else -1.0
    cell_differences: dict[tuple[Hashable, Hashable], float] = {}
    for cell, indices in cell_indices.items():
        candidate_summary = float(
            within_cell_statistic([candidate_values[index] for index in indices])
        )
        baseline_summary = float(
            within_cell_statistic([baseline_values[index] for index in indices])
        )
        if not isfinite(candidate_summary) or not isfinite(baseline_summary):
            raise ValueError("within_cell_statistic returned a non-finite value")
        cell_differences[cell] = direction * (candidate_summary - baseline_summary)

    ordered_differences = [
        cell_differences[(donor, block)]
        for donor in donor_order
        for block in block_order
    ]
    estimate = float(between_cell_statistic(ordered_differences))
    if not isfinite(estimate):
        raise ValueError("between_cell_statistic returned a non-finite value")
    null_differences = {
        cell: difference - estimate for cell, difference in cell_differences.items()
    }

    rng = Random(seed)
    bootstrap_estimates: list[float] = []
    null_estimates: list[float] = []
    for _ in range(n_resamples):
        sampled_donors = [
            donor_order[rng.randrange(len(donor_order))]
            for _ in range(len(donor_order))
        ]
        sampled_blocks = [
            block_order[rng.randrange(len(block_order))]
            for _ in range(len(block_order))
        ]
        sample = [
            cell_differences[(donor, block)]
            for donor in sampled_donors
            for block in sampled_blocks
        ]
        bootstrap_estimate = float(between_cell_statistic(sample))
        if not isfinite(bootstrap_estimate):
            raise ValueError(
                "between_cell_statistic returned a non-finite bootstrap value"
            )
        bootstrap_estimates.append(bootstrap_estimate)
        null_estimate = float(
            between_cell_statistic(
                [
                    null_differences[(donor, block)]
                    for donor in sampled_donors
                    for block in sampled_blocks
                ]
            )
        )
        if not isfinite(null_estimate):
            raise ValueError(
                "between_cell_statistic returned a non-finite null value"
            )
        null_estimates.append(null_estimate)

    bootstrap_estimates.sort()
    tail = (1.0 - confidence_level) / 2.0
    lower = _percentile(bootstrap_estimates, tail)
    upper = _percentile(bootstrap_estimates, 1.0 - tail)
    bootstrap_p_value = (
        sum(abs(value) >= abs(estimate) for value in null_estimates) + 1
    ) / (n_resamples + 1)
    probability_improvement = sum(
        value > 0.0 for value in bootstrap_estimates
    ) / n_resamples

    return TwoWayBootstrapInterval(
        estimate=estimate,
        lower=lower,
        upper=upper,
        confidence_level=confidence_level,
        n_resamples=n_resamples,
        n_donors=len(donor_order),
        n_genomic_blocks=len(block_order),
        n_cells=len(cell_indices),
        n_observations=len(candidate_values),
        seed=seed,
        bootstrap_p_value=bootstrap_p_value,
        probability_improvement=probability_improvement,
    )


def ld_block_paired_bootstrap(
    candidate: Sequence[float],
    baseline: Sequence[float],
    ld_block_ids: Sequence[Hashable],
    *,
    n_resamples: int = 10_000,
    confidence_level: float = 0.95,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
    higher_is_better: bool = True,
    minimum_ld_blocks: int = 2,
    expected_n_ld_blocks: int | None = None,
    within_block_statistic: Callable[[Sequence[float]], float] = fmean,
    between_block_statistic: Callable[[Sequence[float]], float] = fmean,
) -> BootstrapInterval:
    """Strict variant-only specialization of the paired cluster bootstrap.

    Variants or loci within an LD block are reduced before resampling. The LD
    blocks, not the variant rows, are the independent units.
    """

    blocks = _validated_unit_ids(ld_block_ids, "ld_block_ids")
    unique_blocks = len(dict.fromkeys(blocks))
    _validate_independent_unit_count(
        unique_blocks,
        name="LD blocks",
        minimum=minimum_ld_blocks,
        expected=expected_n_ld_blocks,
    )
    return paired_cluster_bootstrap(
        candidate,
        baseline,
        blocks,
        n_resamples=n_resamples,
        confidence_level=confidence_level,
        seed=seed,
        higher_is_better=higher_is_better,
        within_cluster_statistic=within_block_statistic,
        between_cluster_statistic=between_block_statistic,
    )


# Match the frozen inference-policy identifier while retaining the explicit
# paired name above.
ld_block_bootstrap = ld_block_paired_bootstrap


def five_seed_direction_stability(
    candidate: Sequence[float],
    baseline: Sequence[float],
    seed_ids: Sequence[int] = DEFAULT_TRAINING_SEEDS,
    *,
    higher_is_better: bool = True,
) -> SeedDirectionStabilityResult:
    """Require candidate improvement in at least four of five training seeds.

    This is a deterministic optimization-stability diagnostic. Training seeds
    are repeated fits of the same biological data, not biological replicates.
    The result must not be used to calculate a p-value, confidence interval, or
    inferential sample size; those require donor- or genomic-block-level data.
    """

    candidate_values = _coerce_finite_series(candidate, "candidate")
    baseline_values = _coerce_finite_series(baseline, "baseline")
    seeds = list(seed_ids)
    if len(candidate_values) != len(baseline_values):
        raise ValueError("candidate and baseline must have equal lengths")
    if len(candidate_values) != 5 or len(seeds) != 5:
        raise ValueError("direction stability requires exactly five aligned seeds")
    if any(not isinstance(seed, int) or isinstance(seed, bool) for seed in seeds):
        raise TypeError("seed_ids must contain integers")
    if len(set(seeds)) != 5:
        raise ValueError("seed_ids must contain five unique seeds")
    if not isinstance(higher_is_better, bool):
        raise TypeError("higher_is_better must be boolean")

    direction = 1.0 if higher_is_better else -1.0
    differences = tuple(
        direction * (candidate_value - baseline_value)
        for candidate_value, baseline_value in zip(
            candidate_values, baseline_values, strict=True
        )
    )
    qualifying_count = sum(difference > 0.0 for difference in differences)
    required_count = 4
    passed = qualifying_count >= required_count
    reason = (
        f"candidate improves in {qualifying_count}/5 training seeds; "
        "seeds are optimization repeats, not biological replicates"
    )
    return SeedDirectionStabilityResult(
        passed=passed,
        qualifying_seed_count=qualifying_count,
        evaluated_seed_count=5,
        required_seed_count=required_count,
        seed_ids=tuple(seeds),
        oriented_differences=differences,
        seeds_are_biological_replicates=False,
        supports_inferential_test=False,
        reason=reason,
    )


def prospective_power_gate(
    *,
    n_units: int,
    paired_sd: float,
    minimum_effect: float,
    alpha: float = 0.05,
    power: float = 0.80,
    n_primary_claims: int = 1,
    min_units: int = 10,
    two_sided: bool = True,
) -> PowerGateResult:
    """Check whether a sealed endpoint can resolve its practical effect floor.

    ``effective_alpha`` is the first-step Holm/Bonferroni threshold
    ``alpha / n_primary_claims``.  This is deliberately conservative for
    prospective admission: later Holm steps can only be less stringent.
    """

    if not isinstance(n_units, int) or isinstance(n_units, bool) or n_units < 2:
        raise ValueError("n_units must be an integer of at least two")
    if not isinstance(min_units, int) or isinstance(min_units, bool) or min_units < 2:
        raise ValueError("min_units must be an integer of at least two")
    if (
        not isinstance(n_primary_claims, int)
        or isinstance(n_primary_claims, bool)
        or n_primary_claims < 1
    ):
        raise ValueError("n_primary_claims must be a positive integer")
    for name, value in (
        ("paired_sd", paired_sd),
        ("minimum_effect", minimum_effect),
        ("alpha", alpha),
        ("power", power),
    ):
        if isinstance(value, bool):
            raise TypeError(f"{name} must be numeric, not boolean")
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{name} must be numeric") from exc
        if not isfinite(number):
            raise ValueError(f"{name} must be finite")
    paired_sd = float(paired_sd)
    minimum_effect = float(minimum_effect)
    alpha = float(alpha)
    power = float(power)
    if paired_sd < 0.0:
        raise ValueError("paired_sd cannot be negative")
    if minimum_effect <= 0.0:
        raise ValueError("minimum_effect must be positive")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be between zero and one")
    if not 0.5 < power < 1.0:
        raise ValueError("power must be greater than 0.5 and less than one")
    if not isinstance(two_sided, bool):
        raise TypeError("two_sided must be boolean")

    effective_alpha = alpha / n_primary_claims
    tail_divisor = 2.0 if two_sided else 1.0
    normal = NormalDist()
    critical_alpha = normal.inv_cdf(1.0 - effective_alpha / tail_divisor)
    critical_power = normal.inv_cdf(power)
    detectable_effect = (critical_alpha + critical_power) * paired_sd / sqrt(n_units)
    enough_units = n_units >= min_units
    enough_power = detectable_effect <= minimum_effect
    passed = enough_units and enough_power
    if not enough_units:
        reason = f"{n_units} independent units is below the locked floor of {min_units}"
    elif not enough_power:
        reason = (
            f"minimum detectable effect {detectable_effect:.6g} exceeds the "
            f"locked practical effect {minimum_effect:.6g}"
        )
    else:
        reason = (
            f"minimum detectable effect {detectable_effect:.6g} is within the "
            f"locked practical effect {minimum_effect:.6g}"
        )

    return PowerGateResult(
        passed=passed,
        n_units=n_units,
        paired_sd=paired_sd,
        minimum_effect=minimum_effect,
        minimum_detectable_effect=detectable_effect,
        alpha=alpha,
        effective_alpha=effective_alpha,
        power=power,
        n_primary_claims=n_primary_claims,
        min_units=min_units,
        two_sided=two_sided,
        reason=reason,
    )
