"""Answer the frozen aspect-separability question, with the prespecification proved first.

Stage 0 of the MASLD showcase model, job 2 of 2. The wrapper verifies the
frozen prespecification tree with ``sha256sum --check --strict`` before this
runs; this script then re-reads the criteria and their thresholds *out of that
frozen file* rather than restating them, so a threshold cannot drift between
the freeze and the answer without aborting the run.

Three criteria decide whether steatosis, hepatocyte ballooning and lobular
inflammation are three outcome axes or one severity axis:

C1  every pairwise Spearman among the three aspects has ``|rho| < 0.80``
C2  the participation ratio of the 3x3 Spearman eigenvalues is ``>= 2.0``
C3  at least one of the three per-gene aspect-association vector pairs has
    ``|rho| < 0.80`` -- the decisive one, because two correlated labels can
    still carry distinguishable gene signatures and two weakly correlated
    labels can still rank genes identically

Two diagnostics are reported beside the criteria and never inside them. Both
exist because a C3 pass is the easiest of the three to manufacture:

*What the C3 statistic reads with no shared aspect signal.*  The two
association vectors in any C3 pair are estimated from the *same* 99
participants and the *same* expression matrix, so their sampling errors are
shared rather than independent. A gene that is noisy in this particular sample
is noisy for both aspects at once, which pushes the C3 statistic up even when
the two aspects have nothing in common. The floor is therefore measured, by
permuting the two aspect vectors independently and correlating the association
vectors that result. Because the errors are shared, a Spearman-Brown
attenuation correction does not apply here and is deliberately not computed;
split-half reliability is still reported, but as a statement about how
reproducible one aspect's gene ordering is across independent participants,
never as a ceiling to divide by.

*What the C3 statistic reads for axes that really are different.*  Fibrosis,
lobular necrosis, the deposited component sum and recorded sex are put through
the identical pipeline. Fibrosis is the informative one: it is a genuinely
different histological axis, so its association-vector correlation with the
three aspects says whether a C3 value near 0.87 means "these three are alike"
or merely "any two liver axes read that high on this substrate", which would
make the 0.80 threshold unreachable by construction.

*Severity residualization.*  Each aspect rank vector residualized on the rank
of the deposited NASH-CRN component sum, then re-associated. This asks whether
any difference in gene ordering survives removing the shared severity axis.
Residualizing three parts on their own total induces negative dependence among
the residuals by construction, so negative correlations here are expected
arithmetic and are reported with that caveat attached rather than as a
finding.

Nulls are measured, never assumed. Every null is drawn from the realised
vectors, so the heavy tie structure of these scores -- 24 participants with an
all-zero component sum, three- and four-level ordinals -- is carried into the
null rather than idealised away. The multiplicity-corrected detectable ``|rho|``
is a max-statistic permutation floor over the realised gene universe, which
accounts for gene-gene correlation; the analytic ``1/sqrt(n-1)`` reference is
reported beside it, not instead of it.

``lobular_necrosis`` is deposited in the endpoint table but is not a NASH-CRN
NAS component and is excluded from the three aspects. Fibrosis is a fourth axis
reported alongside and never folded in.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import rankdata

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.evaluators.auprc_reference import (
    spearman_reference,
    stratified_spearman_reference,
)
from masld_bench.evaluators.metrics import _average_ranks, spearman_correlation
from masld_bench.evaluators.stats import benjamini_hochberg


class SeparabilityError(RuntimeError):
    """Raised when the substrate cannot support the frozen question."""


#: The three NASH-CRN NAS components. lobular_necrosis is deliberately absent.
ASPECTS = ("steatosis", "ballooning", "lobular_inflammation")
PAIRS = (
    ("steatosis", "ballooning"),
    ("steatosis", "lobular_inflammation"),
    ("ballooning", "lobular_inflammation"),
)
#: Reported alongside the three, never folded into them.
CONTEXT_AXES = ("fibrosis", "lobular_necrosis")

C1_THRESHOLD = 0.80
C2_THRESHOLD = 2.0
C3_THRESHOLD = 0.80

#: Literal threshold prose that must still be in the frozen prespecification.
FROZEN_THRESHOLDS = {
    "c1_labels_are_not_redundant": "every pair |rho| < 0.80",
    "c2_effective_dimensionality": ">= 2.0 of a possible 3.0",
    "c3_aspects_rank_genes_differently": "at least one pair with |rho| < 0.80",
}


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def load_prespecification(directory: Path) -> dict:
    """Read the frozen criteria and refuse to run against a drifted freeze."""

    verify_frozen_tree(directory)
    payload = json.loads(
        (directory / "aspect_separability_prespecification.json").read_text(
            encoding="utf-8"
        )
    )
    criteria = payload["criteria"]
    if sorted(criteria) != sorted(FROZEN_THRESHOLDS):
        raise SeparabilityError("the frozen criteria are not the three expected ones")
    for name, expected in FROZEN_THRESHOLDS.items():
        if criteria[name]["threshold"] != expected:
            raise SeparabilityError(
                f"{name} threshold drifted since the freeze: "
                f"{criteria[name]['threshold']!r} is not {expected!r}"
            )
    if criteria["c3_aspects_rank_genes_differently"].get("decisive") is not True:
        raise SeparabilityError("c3 lost its decisive mark since the freeze")
    return payload


def verify_recorded_marginals(
    rows: Sequence[dict[str, str]], prespec: dict
) -> dict[str, object]:
    """The prespecification's recorded marginals must be the deposited ones."""

    inspected = prespec["what_was_inspected_before_freezing"]
    checked: dict[str, object] = {}
    for column in ASPECTS + ("fibrosis",):
        observed = Counter(row[column] for row in rows)
        recorded = {key: int(value) for key, value in inspected[column].items()}
        if dict(observed) != recorded:
            raise SeparabilityError(
                f"{column} marginal differs from the frozen prespecification"
            )
        checked[column] = dict(sorted(observed.items()))
    zeros = sum(1 for row in rows if row["nash_crn_component_sum"] == "0")
    if zeros != inspected["nash_crn_component_sum_zero_participants"]:
        raise SeparabilityError("the all-zero tied block differs from the freeze")
    checked["nash_crn_component_sum_zero_participants"] = zeros
    return checked


def verify_nas_identity(rows: Sequence[dict[str, str]]) -> dict[str, object]:
    """The component sum is the three aspects, and necrosis is not one of them."""

    three = 0
    with_necrosis = 0
    for row in rows:
        deposited = int(row["nash_crn_component_sum"])
        total = sum(int(row[aspect]) for aspect in ASPECTS)
        three += total == deposited
        with_necrosis += total + int(row["lobular_necrosis"]) == deposited
    if three != len(rows):
        raise SeparabilityError("the component sum is not the three aspects")
    return {
        "component_sum_equals_three_aspects": three,
        "component_sum_equals_three_aspects_plus_necrosis": with_necrosis,
        "adding_necrosis_breaks_the_identity_for": len(rows) - with_necrosis,
        "lobular_necrosis_is_a_nas_component": False,
    }


# --------------------------------------------------------------------------
# tie structure
# --------------------------------------------------------------------------


def tie_structure(values: Sequence[float], name: str) -> dict[str, object]:
    """Describe how much of the ranking is decided by ties, not by order."""

    counts = Counter(values)
    n = len(values)
    total_pairs = n * (n - 1) // 2
    tied_pairs = sum(count * (count - 1) // 2 for count in counts.values())
    return {
        "axis": name,
        "n_participants": n,
        "distinct_values": len(counts),
        "value_counts": {str(key): counts[key] for key in sorted(counts)},
        "largest_tied_block": max(counts.values()),
        "largest_tied_block_fraction": max(counts.values()) / n,
        "tied_pairs": tied_pairs,
        "total_pairs": total_pairs,
        "fraction_of_pairs_tied": tied_pairs / total_pairs,
    }


def joint_tie_structure(
    rows: Sequence[dict[str, str]],
) -> dict[str, object]:
    """Ties in the three-aspect profile, which is what a model must separate."""

    triples = Counter(tuple(int(row[aspect]) for aspect in ASPECTS) for row in rows)
    n = len(rows)
    tied_pairs = sum(count * (count - 1) // 2 for count in triples.values())
    ordered = sorted(triples.items(), key=lambda item: (-item[1], item[0]))
    return {
        "profile": "(steatosis, ballooning, lobular_inflammation)",
        "n_participants": n,
        "distinct_profiles": len(triples),
        "largest_tied_block": ordered[0][1],
        "largest_tied_profile": list(ordered[0][0]),
        "all_zero_profile_count": triples.get((0, 0, 0), 0),
        "singleton_profiles": sum(1 for count in triples.values() if count == 1),
        "fraction_of_pairs_with_an_identical_profile": tied_pairs / (n * (n - 1) // 2),
        "most_common_profiles": [
            {"profile": list(profile), "n": count} for profile, count in ordered[:8]
        ],
    }


# --------------------------------------------------------------------------
# ranks
# --------------------------------------------------------------------------


def average_ranks(matrix: np.ndarray) -> np.ndarray:
    """Average ranks down each column, matching the package's tie handling."""

    return rankdata(matrix, method="average", axis=0)


def validate_rank_agreement(
    matrix: np.ndarray, columns: Sequence[int]
) -> dict[str, object]:
    """Prove the vectorized ranks equal the package's own, on real columns."""

    ranked = average_ranks(matrix[:, list(columns)])
    worst = 0.0
    for offset, column in enumerate(columns):
        reference = np.asarray(_average_ranks(list(matrix[:, column])), dtype=float)
        worst = max(worst, float(np.max(np.abs(ranked[:, offset] - reference))))
    return {
        "columns_checked": len(columns),
        "max_absolute_rank_difference": worst,
        "agrees_with_masld_bench_average_ranks": worst == 0.0,
    }


def standardize(ranked: np.ndarray) -> np.ndarray:
    """Centre and scale ranked columns so a dot product is a correlation."""

    centred = ranked - ranked.mean(axis=0, keepdims=True)
    norms = np.sqrt((centred**2).sum(axis=0, keepdims=True))
    if np.any(norms == 0):
        raise SeparabilityError("a constant column reached the correlation step")
    return centred / norms


def association_vector(
    standardized_genes: np.ndarray, outcome: np.ndarray
) -> np.ndarray:
    """Per-gene Spearman between expression and one outcome vector."""

    ranked = rankdata(outcome, method="average").astype(float)
    centred = ranked - ranked.mean()
    norm = math.sqrt(float((centred**2).sum()))
    if norm == 0.0:
        raise SeparabilityError("a constant outcome reached the association step")
    return standardized_genes.T @ (centred / norm)


# --------------------------------------------------------------------------
# criteria
# --------------------------------------------------------------------------


def permutation_pair_null(
    left: Sequence[float],
    right: Sequence[float],
    *,
    n_draws: int,
    seed: int,
) -> dict[str, object]:
    """Null for two heavily tied vectors: permute one, keep both tie structures.

    ``spearman_reference`` answers a different question -- what a *continuous*
    random scorer reaches against this outcome -- and is reported beside this,
    because the two references differ exactly by the tie structure of the
    permuted vector and the gap is the thing worth seeing.
    """

    generator = np.random.default_rng(seed)
    left_ranked = rankdata(np.asarray(left, dtype=float), method="average")
    right_ranked = rankdata(np.asarray(right, dtype=float), method="average")
    left_standard = (left_ranked - left_ranked.mean()) / np.sqrt(
        ((left_ranked - left_ranked.mean()) ** 2).sum()
    )
    right_centred = right_ranked - right_ranked.mean()
    right_norm = np.sqrt((right_centred**2).sum())
    observed = float(left_standard @ (right_centred / right_norm))
    draws = np.empty(n_draws, dtype=float)
    working = right_centred / right_norm
    for index in range(n_draws):
        draws[index] = left_standard @ generator.permutation(working)
    absolute = np.abs(draws)
    extreme = int(np.sum(absolute >= abs(observed)))
    return {
        "null": "permutation_of_one_realised_aspect_vector_ties_preserved",
        "observed_spearman": observed,
        "null_mean_absolute": float(absolute.mean()),
        "null_sd": float(draws.std(ddof=1)),
        "null_absolute_percentile_95": float(np.percentile(absolute, 95)),
        "null_absolute_percentile_99": float(np.percentile(absolute, 99)),
        "n_draws": n_draws,
        "seed": seed,
        "p_value_two_sided": (extreme + 1) / (n_draws + 1),
    }


def participation_ratio(matrix: np.ndarray) -> tuple[float, list[float]]:
    eigenvalues = np.linalg.eigvalsh(matrix)
    ratio = float(eigenvalues.sum() ** 2 / (eigenvalues**2).sum())
    return ratio, [float(value) for value in sorted(eigenvalues, reverse=True)]


def split_half_reliability(
    raw_genes: np.ndarray,
    outcomes: dict[str, np.ndarray],
    *,
    n_splits: int,
    seed: int,
) -> dict[str, dict[str, float]]:
    """How well an aspect's own gene ordering reproduces on independent halves.

    This is the ceiling every cross-aspect correlation is measured against. It
    is a diagnostic, not one of the frozen criteria.

    Ranks are recomputed inside each half, because a Spearman estimated on 50
    participants is the thing whose reproducibility is in question. A gene that
    varies across all 99 can still be constant inside a half; those genes are
    dropped for that split only, on the intersection of the two halves, and the
    surviving count is reported rather than assumed to be the whole universe.
    """

    generator = np.random.default_rng(seed)
    n = raw_genes.shape[0]
    raw: dict[str, list[float]] = {aspect: [] for aspect in outcomes}
    surviving: list[int] = []
    for _ in range(n_splits):
        order = generator.permutation(n)
        halves = (order[: n // 2], order[n // 2 :])
        blocks = []
        keep = np.ones(raw_genes.shape[1], dtype=bool)
        for indices in halves:
            half_values = raw_genes[indices, :]
            keep &= half_values.min(axis=0) != half_values.max(axis=0)
            blocks.append((indices, half_values))
        if not keep.any():
            continue
        surviving.append(int(keep.sum()))
        standardized_blocks = [
            (indices, standardize(average_ranks(half_values[:, keep])))
            for indices, half_values in blocks
        ]
        for aspect, outcome in outcomes.items():
            vectors = []
            usable = True
            for indices, block in standardized_blocks:
                values = outcome[indices]
                if len(set(values.tolist())) < 2:
                    usable = False
                    break
                vectors.append(association_vector(block, values))
            if not usable:
                continue
            raw[aspect].append(spearman_of_vectors(vectors[0], vectors[1]))
    summary: dict[str, dict[str, float]] = {}
    median_surviving = float(np.median(surviving)) if surviving else float("nan")
    for aspect, values in raw.items():
        array = np.asarray(values, dtype=float)
        if array.size == 0:
            summary[aspect] = {"n_usable_splits": 0}
            continue
        half = float(np.median(array))
        # Spearman-Brown correction from half length back to the full cohort.
        full = 2.0 * half / (1.0 + half) if half > -1.0 else float("nan")
        summary[aspect] = {
            "n_usable_splits": int(array.size),
            "median_genes_non_constant_in_both_halves": median_surviving,
            "median_half_length_reliability": half,
            "percentile_5": float(np.percentile(array, 5)),
            "percentile_95": float(np.percentile(array, 95)),
            "spearman_brown_full_length_reliability": full,
        }
    return summary


def spearman_of_vectors(left: np.ndarray, right: np.ndarray) -> float:
    """Spearman between two long vectors, via ranks and a dot product."""

    left_ranked = rankdata(left, method="average")
    right_ranked = rankdata(right, method="average")
    left_centred = left_ranked - left_ranked.mean()
    right_centred = right_ranked - right_ranked.mean()
    denominator = math.sqrt(
        float((left_centred**2).sum()) * float((right_centred**2).sum())
    )
    if denominator == 0.0:
        raise SeparabilityError("a constant association vector reached the comparison")
    return float(left_centred @ right_centred / denominator)


def bootstrap_label_spearman(
    left: np.ndarray,
    right: np.ndarray,
    *,
    n_resamples: int,
    seed: int,
) -> dict[str, object]:
    """Participant-level percentile interval for one pairwise label Spearman.

    The unit resampled is the participant, which is the unit of inference here.
    A resample whose labels go constant carries no correlation and is counted
    rather than dropped silently.
    """

    generator = np.random.default_rng(seed)
    n = left.size
    draws: list[float] = []
    degenerate = 0
    for _ in range(n_resamples):
        indices = generator.integers(0, n, size=n)
        sample_left = left[indices]
        sample_right = right[indices]
        if sample_left.min() == sample_left.max() or (
            sample_right.min() == sample_right.max()
        ):
            degenerate += 1
            continue
        draws.append(
            spearman_correlation(sample_left.tolist(), sample_right.tolist())
        )
    array = np.abs(np.asarray(draws, dtype=float))
    return {
        "resampling_unit": "participant",
        "n_resamples": n_resamples,
        "n_degenerate_resamples": degenerate,
        "seed": seed,
        "abs_spearman_percentile_2_5": float(np.percentile(array, 2.5)),
        "abs_spearman_percentile_97_5": float(np.percentile(array, 97.5)),
        "fraction_of_resamples_below_the_threshold": float(
            np.mean(array < C1_THRESHOLD)
        ),
    }


def bootstrap_association_vector_spearman(
    raw_genes: np.ndarray,
    outcomes: dict[str, np.ndarray],
    *,
    n_resamples: int,
    seed: int,
) -> dict[str, dict[str, object]]:
    """Participant-level percentile interval for every C3 pair at once.

    One resample of participants produces all three association vectors, so the
    three intervals share their resamples and stay comparable. Ranks are
    recomputed inside each resample; genes that go constant there are dropped
    for that resample only, and the surviving count is reported.
    """

    generator = np.random.default_rng(seed)
    n = raw_genes.shape[0]
    draws: dict[tuple[str, str], list[float]] = {pair: [] for pair in PAIRS}
    surviving: list[int] = []
    skipped = 0
    for _ in range(n_resamples):
        indices = generator.integers(0, n, size=n)
        block = raw_genes[indices, :]
        keep = block.min(axis=0) != block.max(axis=0)
        usable = {
            aspect: values[indices]
            for aspect, values in outcomes.items()
        }
        if not keep.any() or any(
            values.min() == values.max() for values in usable.values()
        ):
            skipped += 1
            continue
        surviving.append(int(keep.sum()))
        standardized = standardize(average_ranks(block[:, keep]))
        vectors = {
            aspect: association_vector(standardized, values)
            for aspect, values in usable.items()
        }
        for pair in PAIRS:
            draws[pair].append(
                spearman_of_vectors(vectors[pair[0]], vectors[pair[1]])
            )
    summary: dict[str, dict[str, object]] = {}
    for pair, values in draws.items():
        array = np.abs(np.asarray(values, dtype=float))
        summary["|".join(pair)] = {
            "resampling_unit": "participant",
            "n_resamples": n_resamples,
            "n_usable_resamples": int(array.size),
            "n_skipped_resamples": skipped,
            "median_genes_non_constant_in_the_resample": (
                float(np.median(surviving)) if surviving else float("nan")
            ),
            "seed": seed,
            "abs_spearman_percentile_2_5": float(np.percentile(array, 2.5)),
            "abs_spearman_percentile_97_5": float(np.percentile(array, 97.5)),
            "fraction_of_resamples_below_the_threshold": float(
                np.mean(array < C3_THRESHOLD)
            ),
        }
    return summary


def shared_sample_c3_floor(
    standardized_genes: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    *,
    n_draws: int,
    seed: int,
) -> dict[str, object]:
    """What C3 reads when two outcomes share the sample but share no signal.

    Both association vectors in a C3 pair come from the same participants and
    the same expression matrix, so a gene that is noisy in this sample is noisy
    for both aspects at once. That shared error inflates C3 on its own. Drawing
    the floor by permuting each aspect independently -- keeping each one's own
    tie structure -- measures how much, instead of assuming it away.
    """

    generator = np.random.default_rng(seed)
    draws = np.empty(n_draws, dtype=float)
    for index in range(n_draws):
        draws[index] = spearman_of_vectors(
            association_vector(standardized_genes, generator.permutation(left)),
            association_vector(standardized_genes, generator.permutation(right)),
        )
    absolute = np.abs(draws)
    return {
        "null": "independent_permutation_of_both_aspects_same_sample_same_matrix",
        "n_draws": n_draws,
        "seed": seed,
        "mean_abs": float(absolute.mean()),
        "sd": float(draws.std(ddof=1)),
        "abs_percentile_95": float(np.percentile(absolute, 95)),
        "abs_percentile_99": float(np.percentile(absolute, 99)),
        "abs_maximum": float(absolute.max()),
    }


def max_statistic_floor(
    standardized_genes: np.ndarray,
    outcome: np.ndarray,
    observed: np.ndarray,
    *,
    n_permutations: int,
    block: int,
    seed: int,
) -> dict[str, object]:
    """Family-wise detectable ``|rho|`` over the realised gene universe.

    Permuting the *outcome* and taking the maximum absolute association across
    every realised gene gives a floor that already accounts for gene-gene
    correlation, which a per-gene analytic threshold cannot. Per-gene
    exceedance counts from the same draws give the permutation p-values.
    """

    generator = np.random.default_rng(seed)
    ranked = rankdata(outcome, method="average").astype(float)
    centred = ranked - ranked.mean()
    unit = centred / math.sqrt(float((centred**2).sum()))
    maxima: list[float] = []
    exceedances = np.zeros(standardized_genes.shape[1], dtype=np.int64)
    absolute_observed = np.abs(observed)
    drawn = 0
    while drawn < n_permutations:
        size = min(block, n_permutations - drawn)
        permuted = np.empty((size, unit.size), dtype=float)
        for index in range(size):
            permuted[index] = generator.permutation(unit)
        null_block = np.abs(permuted @ standardized_genes)
        maxima.extend(null_block.max(axis=1).tolist())
        exceedances += (null_block >= absolute_observed).sum(axis=0)
        drawn += size
    maxima_array = np.asarray(maxima, dtype=float)
    p_values = (exceedances + 1) / (n_permutations + 1)
    adjusted = np.asarray(benjamini_hochberg(p_values.tolist()), dtype=float)
    discovered = adjusted < 0.05
    return {
        "null": "outcome_permutation_max_statistic_over_realised_gene_universe",
        "n_permutations": n_permutations,
        "seed": seed,
        "n_genes": int(standardized_genes.shape[1]),
        "familywise_detectable_abs_rho_p95": float(np.percentile(maxima_array, 95)),
        "familywise_detectable_abs_rho_p99": float(np.percentile(maxima_array, 99)),
        "max_statistic_null_mean": float(maxima_array.mean()),
        "analytic_rank_null_sd_1_over_sqrt_n_minus_1": 1.0 / math.sqrt(
            outcome.size - 1
        ),
        "observed_max_abs_rho": float(absolute_observed.max()),
        "genes_bh_below_0_05": int(discovered.sum()),
        "smallest_abs_rho_reaching_bh_0_05": (
            float(absolute_observed[discovered].min()) if discovered.any() else None
        ),
        "genes_above_the_familywise_floor": int(
            (absolute_observed > np.percentile(maxima_array, 95)).sum()
        ),
        "permutation_p_value_resolution": 1.0 / (n_permutations + 1),
        "bh_is_reachable_only_if_at_least_k_genes_sit_at_that_floor": math.ceil(
            standardized_genes.shape[1] / (0.05 * (n_permutations + 1))
        ),
        "genes_at_the_permutation_p_floor": int(
            np.sum(p_values <= 1.0 / (n_permutations + 1))
        ),
    }


# --------------------------------------------------------------------------
# strata
# --------------------------------------------------------------------------


def assess_strata(
    rows: Sequence[dict[str, str]],
    aspect_values: dict[str, np.ndarray],
    *,
    n_draws: int,
    seed: int,
) -> dict[str, object]:
    """Decide whether a genuine nuisance stratum exists, and say so either way."""

    folds = [row["outer_fold"] for row in rows]
    sexes = [row["recorded_sex"] for row in rows]
    fold_by_stage: dict[str, Counter] = {}
    for row in rows:
        fold_by_stage.setdefault(row["outer_fold"], Counter())[row["stage5"]] += 1

    sex_binary = np.asarray([1.0 if value == "M" else 0.0 for value in sexes])
    sex_association = {}
    for aspect, values in aspect_values.items():
        sex_association[aspect] = permutation_pair_null(
            sex_binary.tolist(), values.tolist(), n_draws=n_draws, seed=seed
        )

    return {
        "outer_fold": {
            "counts": dict(sorted(Counter(folds).items())),
            "stage5_composition_by_fold": {
                fold: dict(sorted(counter.items()))
                for fold, counter in sorted(fold_by_stage.items())
            },
            "is_a_genuine_nuisance_stratum": False,
            "why": (
                "outer_fold is a cross-validation partition this campaign "
                "assigned, not an external source of structure. There is one "
                "study, one platform and one processing run here. The folds are "
                "balanced on histology stage by construction, so permuting "
                "within a fold would hold fixed part of the very outcome "
                "variation the criteria are asking about, and would remove "
                "signal rather than nuisance."
            ),
        },
        "recorded_sex": {
            "counts": dict(sorted(Counter(sexes).items())),
            "association_with_each_aspect": sex_association,
            "why": (
                "recorded_sex is the only candidate biological stratum in this "
                "cohort. Whether it is a genuine one is decided by whether it "
                "carries between-stratum structure in the aspects, which the "
                "measured association above answers rather than assumes."
            ),
        },
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prespecification", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--molecular", type=Path, required=True)
    parser.add_argument("--folds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--permutation-block", type=int, default=250)
    parser.add_argument("--label-draws", type=int, default=20000)
    parser.add_argument("--splits", type=int, default=200)
    parser.add_argument("--c3-floor-draws", type=int, default=500)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--c3-bootstrap", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise SeparabilityError("refusing to overwrite a separability result")

    prespec = load_prespecification(arguments.prespecification)
    print("verified the frozen prespecification", flush=True)

    rows = read_table(arguments.endpoints)
    participants = [row["participant_id"] for row in rows]
    if len(set(participants)) != len(participants):
        raise SeparabilityError("participant_id is not unique in the endpoint table")
    recorded = verify_recorded_marginals(rows, prespec)
    identity = verify_nas_identity(rows)
    print(f"endpoints: {len(rows)} rows, {len(set(participants))} participants",
          flush=True)

    folds = {row["participant_id"]: row["outer_fold"] for row in read_table(arguments.folds)}
    if folds != {row["participant_id"]: row["outer_fold"] for row in rows}:
        raise SeparabilityError("the fold table disagrees with the endpoint table")

    aspect_values = {
        aspect: np.asarray([int(row[aspect]) for row in rows], dtype=float)
        for aspect in ASPECTS
    }
    context_values = {
        axis: np.asarray([int(row[axis]) for row in rows], dtype=float)
        for axis in CONTEXT_AXES
    }
    severity = np.asarray(
        [int(row["nash_crn_component_sum"]) for row in rows], dtype=float
    )

    ties = {
        name: tie_structure(values.tolist(), name)
        for name, values in {**aspect_values, **context_values,
                             "nash_crn_component_sum": severity}.items()
    }
    joint = joint_tie_structure(rows)

    strata = assess_strata(
        rows, aspect_values, n_draws=arguments.label_draws, seed=arguments.seed
    )
    sex_carries_structure = any(
        entry["p_value_two_sided"] < 0.05
        for entry in strata["recorded_sex"]["association_with_each_aspect"].values()
    )
    strata["recorded_sex"]["is_a_genuine_nuisance_stratum"] = bool(
        sex_carries_structure
    )
    strata["conclusion"] = (
        "recorded_sex carries measurable between-stratum structure in at least "
        "one aspect, so the stratified null is reported beside the global one "
        "and the stratified null is the primary reference."
        if sex_carries_structure
        else "No genuine nuisance stratum exists in this cohort. outer_fold is an "
        "analysis partition and recorded_sex shows no measurable association "
        "with any aspect. The global permutation null is the primary "
        "reference; the sex-stratified null is reported beside it only to show "
        "the two agree."
    )
    print(f"strata assessed: sex_is_a_stratum={sex_carries_structure}", flush=True)

    # ---------------- C1 ----------------
    sexes = [row["recorded_sex"] for row in rows]
    c1_pairs = []
    for left, right in PAIRS:
        rho = spearman_correlation(
            aspect_values[left].tolist(), aspect_values[right].tolist()
        )
        pair_null = permutation_pair_null(
            aspect_values[left].tolist(),
            aspect_values[right].tolist(),
            n_draws=arguments.label_draws,
            seed=arguments.seed,
        )
        continuous = spearman_reference(
            aspect_values[right].tolist(), n_draws=arguments.label_draws,
            seed=arguments.seed,
        ).to_dict()
        stratified_reference, stratified_observed, stratified_p = (
            stratified_spearman_reference(
                aspect_values[right].tolist(),
                aspect_values[left].tolist(),
                sexes,
                n_draws=arguments.label_draws,
                seed=arguments.seed,
            )
        )
        c1_pairs.append({
            "pair": [left, right],
            "spearman": rho,
            "abs_spearman": abs(rho),
            "meets_c1": abs(rho) < C1_THRESHOLD,
            "participant_bootstrap": bootstrap_label_spearman(
                aspect_values[left],
                aspect_values[right],
                n_resamples=arguments.bootstrap,
                seed=arguments.seed,
            ),
            "global_permutation_null": pair_null,
            "continuous_random_scorer_reference": continuous,
            "sex_stratified_null": {
                **stratified_reference.to_dict(),
                "observed_spearman": stratified_observed,
                "p_value_two_sided": stratified_p,
            },
        })
    c1_pass = all(entry["meets_c1"] for entry in c1_pairs)
    print(f"C1 computed: pass={c1_pass}", flush=True)

    context_correlations = []
    for axis, values in context_values.items():
        for aspect in ASPECTS:
            context_correlations.append({
                "context_axis": axis,
                "aspect": aspect,
                "spearman": spearman_correlation(
                    values.tolist(), aspect_values[aspect].tolist()
                ),
                "inside_the_criteria": False,
            })

    # ---------------- C2 ----------------
    matrix = np.eye(3)
    lookup = {tuple(sorted(entry["pair"])): entry["spearman"] for entry in c1_pairs}
    for i, left in enumerate(ASPECTS):
        for j, right in enumerate(ASPECTS):
            if i < j:
                value = lookup[tuple(sorted((left, right)))]
                matrix[i, j] = matrix[j, i] = value
    ratio, eigenvalues = participation_ratio(matrix)
    sum_of_squares = sum(
        entry["spearman"] ** 2 for entry in c1_pairs
    )
    closed_form = 9.0 / (3.0 + 2.0 * sum_of_squares)
    if abs(closed_form - ratio) > 1e-9:
        raise SeparabilityError("the participation ratio disagrees with its closed form")
    four_axis = np.eye(4)
    axes4 = ASPECTS + ("fibrosis",)
    all_values = {**aspect_values, "fibrosis": context_values["fibrosis"]}
    for i, left in enumerate(axes4):
        for j, right in enumerate(axes4):
            if i < j:
                value = spearman_correlation(
                    all_values[left].tolist(), all_values[right].tolist()
                )
                four_axis[i, j] = four_axis[j, i] = value
    ratio4, eigenvalues4 = participation_ratio(four_axis)
    c2_pass = ratio >= C2_THRESHOLD
    print(f"C2 computed: participation_ratio={ratio:.4f} pass={c2_pass}", flush=True)

    # ---------------- gene association vectors ----------------
    axis_rows = read_table(arguments.molecular / "participant_axis.tsv")
    if [row["participant_id"] for row in axis_rows] != participants:
        raise SeparabilityError(
            "the molecular participant axis is not in the endpoint table's order"
        )
    features = read_table(arguments.molecular / "rna_feature_axis.tsv")
    gene_ids = [row["stable_gene_id"] for row in features]
    if len(set(gene_ids)) != len(gene_ids):
        raise SeparabilityError("stable_gene_id is not unique on the feature axis")
    values = np.load(arguments.molecular / "rna_values.npy")
    if values.shape != (len(participants), len(gene_ids)):
        raise SeparabilityError("the expression matrix does not match its axes")
    print(f"loaded expression {values.shape}", flush=True)

    ordered_values = np.sort(values, axis=0)
    distinct = (ordered_values[1:] != ordered_values[:-1]).sum(axis=0) + 1
    realised = distinct >= 2
    del ordered_values
    universe = np.flatnonzero(realised)
    if universe.size == 0:
        raise SeparabilityError("no gene varies across participants")
    realised_values = values[:, universe]
    realised_ids = [gene_ids[index] for index in universe]
    print(
        f"realised gene universe: {universe.size} of {len(gene_ids)} "
        f"({len(gene_ids) - universe.size} constant across all 99 participants)",
        flush=True,
    )

    generator = np.random.default_rng(arguments.seed)
    sample = sorted(
        generator.choice(universe.size, size=min(300, universe.size), replace=False)
        .tolist()
    )
    rank_agreement = validate_rank_agreement(realised_values, sample)
    if not rank_agreement["agrees_with_masld_bench_average_ranks"]:
        raise SeparabilityError("the vectorized ranks disagree with the package")

    standardized = standardize(average_ranks(realised_values))
    print("ranked and standardized the realised universe", flush=True)

    vectors = {
        aspect: association_vector(standardized, aspect_values[aspect])
        for aspect in ASPECTS
    }
    reference_check = []
    for aspect in ASPECTS:
        worst = 0.0
        for column in sample[:100]:
            direct = spearman_correlation(
                realised_values[:, column].tolist(), aspect_values[aspect].tolist()
            )
            worst = max(worst, abs(direct - float(vectors[aspect][column])))
        reference_check.append({
            "aspect": aspect,
            "genes_checked": 100,
            "max_absolute_difference_vs_masld_bench_spearman": worst,
        })
    if max(entry["max_absolute_difference_vs_masld_bench_spearman"]
           for entry in reference_check) > 1e-10:
        raise SeparabilityError(
            "the vectorized association vector disagrees with the package metric"
        )
    print("association vectors validated against masld_bench.spearman_correlation",
          flush=True)

    # ---------------- C3 ----------------
    c3_pairs = []
    for left, right in PAIRS:
        rho = spearman_of_vectors(vectors[left], vectors[right])
        c3_pairs.append({
            "pair": [left, right],
            "spearman_of_association_vectors": rho,
            "abs_spearman": abs(rho),
            "meets_c3_for_this_pair": abs(rho) < C3_THRESHOLD,
            "pearson_of_association_vectors": float(
                np.corrcoef(vectors[left], vectors[right])[0, 1]
            ),
        })
    c3_pass = any(entry["meets_c3_for_this_pair"] for entry in c3_pairs)
    print(f"C3 computed: pass={c3_pass}", flush=True)

    c3_bootstrap = bootstrap_association_vector_spearman(
        realised_values,
        aspect_values,
        n_resamples=arguments.c3_bootstrap,
        seed=arguments.seed + 2,
    )
    for entry, pair in zip(c3_pairs, PAIRS):
        entry["participant_bootstrap"] = c3_bootstrap["|".join(pair)]
    print("C3 participant bootstrap done", flush=True)

    # ---------------- detectable effect floor ----------------
    floors = {}
    for aspect in ASPECTS:
        floors[aspect] = max_statistic_floor(
            standardized,
            aspect_values[aspect],
            vectors[aspect],
            n_permutations=arguments.permutations,
            block=arguments.permutation_block,
            seed=arguments.seed,
        )
        print(f"detectable-effect floor done for {aspect}", flush=True)

    # ---------------- diagnostics ----------------
    # What C3 reads with no shared aspect signal, on this sample and matrix.
    noise_floor = []
    for entry, (left, right) in zip(c3_pairs, PAIRS):
        floor = shared_sample_c3_floor(
            standardized,
            aspect_values[left],
            aspect_values[right],
            n_draws=arguments.c3_floor_draws,
            seed=arguments.seed,
        )
        noise_floor.append({
            "pair": [left, right],
            "observed_abs_spearman": entry["abs_spearman"],
            **floor,
            "observed_over_the_floor_p95": (
                entry["abs_spearman"] / floor["abs_percentile_95"]
                if floor["abs_percentile_95"] > 0
                else None
            ),
        })
    print("shared-sample C3 floor done", flush=True)

    # What C3 reads for axes that really are different, same pipeline.
    contrast_axes = {
        "fibrosis": context_values["fibrosis"],
        "lobular_necrosis": context_values["lobular_necrosis"],
        "nash_crn_component_sum": severity,
        "recorded_sex": np.asarray(
            [1.0 if row["recorded_sex"] == "M" else 0.0 for row in rows]
        ),
    }
    contrast = []
    for name, values_vector in contrast_axes.items():
        contrast_vector = association_vector(standardized, values_vector)
        for aspect in ASPECTS:
            contrast.append({
                "contrast_axis": name,
                "aspect": aspect,
                "spearman_of_association_vectors": spearman_of_vectors(
                    contrast_vector, vectors[aspect]
                ),
                "inside_the_criteria": False,
            })
    print("contrast-axis association vectors done", flush=True)

    reliability = split_half_reliability(
        realised_values,
        aspect_values,
        n_splits=arguments.splits,
        seed=arguments.seed + 1,
    )
    print("split-half reliability done", flush=True)

    severity_ranked = rankdata(severity, method="average").astype(float)
    design = np.column_stack([np.ones_like(severity_ranked), severity_ranked])
    residual_vectors = {}
    for aspect in ASPECTS:
        ranked = rankdata(aspect_values[aspect], method="average").astype(float)
        coefficients, *_ = np.linalg.lstsq(design, ranked, rcond=None)
        residual_vectors[aspect] = association_vector(
            standardized, ranked - design @ coefficients
        )
    residualized = []
    for left, right in PAIRS:
        residualized.append({
            "pair": [left, right],
            "spearman_of_severity_residualized_association_vectors": (
                spearman_of_vectors(residual_vectors[left], residual_vectors[right])
            ),
        })
    print("severity-residualized sensitivity done", flush=True)

    go = bool(c1_pass and c2_pass and c3_pass)
    if go:
        decision = "GO"
        finding = (
            "The three aspects are separable in this cohort: they are not "
            "redundant, they span more than one effective dimension, and at "
            "least one pair orders genes differently."
        )
    else:
        failed = [
            name
            for name, passed in (
                ("c1_labels_are_not_redundant", c1_pass),
                ("c2_effective_dimensionality", c2_pass),
                ("c3_aspects_rank_genes_differently", c3_pass),
            )
            if not passed
        ]
        decision = "STOP"
        if "c3_aspects_rank_genes_differently" in failed and len(failed) == 1:
            finding = (
                "The aspects are distinct scores that rank genes identically. "
                "An aspect-resolved model cannot beat a severity model, so "
                "disease_aspect is indeterminate."
            )
        else:
            finding = (
                "The three aspects collapse to one severity axis. An "
                "aspect-resolved model is undefined and every disease_aspect "
                "field is indeterminate."
            )

    payload = {
        "schema_version": "masld-bench-aspect-separability-result-v1",
        "prespec_id": prespec["prespec_id"],
        "decision": decision,
        "finding": finding,
        "a_stop_is_a_result_not_a_failure": True,
        "cohort": {
            "dataset_id": prespec["cohort"]["dataset_id"],
            "endpoint_rows": len(rows),
            "distinct_participants": len(set(participants)),
            "unit_the_arithmetic_uses": "participant",
            "a_row_count_is_not_a_unit_count": (
                "Here they agree: 99 rows are 99 distinct participants, one "
                "biopsy each. Every correlation below is over 99 participants."
            ),
        },
        "recorded_marginals_reproduced": recorded,
        "nas_component_identity": identity,
        "tie_structure": {
            "per_axis": ties,
            "joint_profile": joint,
            "why_it_matters": (
                "Every null below is drawn from these realised vectors, so the "
                "ties are in the null. A null assumed from n alone would be "
                "wrong for scores this coarse."
            ),
        },
        "strata": strata,
        "criteria": {
            "c1_labels_are_not_redundant": {
                "threshold": FROZEN_THRESHOLDS["c1_labels_are_not_redundant"],
                "pairs": c1_pairs,
                "passed": c1_pass,
            },
            "c2_effective_dimensionality": {
                "threshold": FROZEN_THRESHOLDS["c2_effective_dimensionality"],
                "spearman_matrix": matrix.tolist(),
                "eigenvalues": eigenvalues,
                "participation_ratio": ratio,
                "participation_ratio_closed_form": closed_form,
                "sum_of_squared_pairwise_spearman": sum_of_squares,
                "boundary": (
                    "PR >= 2.0 is exactly r12^2 + r13^2 + r23^2 <= 0.75, so for "
                    "homogeneous correlations C2 binds at |rho| = 0.5 and is "
                    "stricter than C1's 0.80."
                ),
                "passed": c2_pass,
                "four_axis_including_fibrosis": {
                    "axes": list(axes4),
                    "spearman_matrix": four_axis.tolist(),
                    "eigenvalues": eigenvalues4,
                    "participation_ratio": ratio4,
                    "inside_the_criteria": False,
                },
            },
            "c3_aspects_rank_genes_differently": {
                "threshold": FROZEN_THRESHOLDS["c3_aspects_rank_genes_differently"],
                "gene_universe": {
                    "features_on_axis": len(gene_ids),
                    "realised_non_constant": int(universe.size),
                    "constant_across_all_participants": len(gene_ids)
                    - int(universe.size),
                    "join_key": "stable_gene_id",
                },
                "association_definition": (
                    "per-gene Spearman between expression and the aspect score "
                    "across all 99 participants, marginal"
                ),
                "pairs": c3_pairs,
                "passed": c3_pass,
                "decisive": True,
            },
        },
        "detectable_effect_floor": {
            "per_aspect": floors,
            "analytic_reference": {
                "one_over_sqrt_n_minus_1": 1.0 / math.sqrt(len(rows) - 1),
                "note": (
                    "The analytic figure is a per-gene rank-null sd, not a "
                    "multiplicity-corrected floor. The max-statistic figure "
                    "above is the one to quote over this gene universe."
                ),
            },
        },
        "diagnostics_not_criteria": {
            "why": (
                "A C3 pass is the easiest of the three to manufacture, because "
                "a noisy association vector correlates poorly with anything, "
                "including a second measurement of itself. These bound how "
                "much of any observed difference in gene ordering could be "
                "sampling noise."
            ),
            "split_half_reliability": {
                "per_aspect": reliability,
                "what_it_is": (
                    "how well one aspect's own gene ordering reproduces on two "
                    "independent halves of the cohort, Spearman-Brown corrected "
                    "back to full length"
                ),
                "what_it_is_not": (
                    "an attenuation ceiling. The two association vectors in a "
                    "C3 pair come from the same participants and the same "
                    "expression matrix, so their errors are shared, not "
                    "independent, and the classical disattenuation formula does "
                    "not apply. Dividing by this reliability returns values "
                    "above 1, which is the arithmetic saying so. The measured "
                    "shared-sample floor below is the correct reference."
                ),
            },
            "shared_sample_c3_floor": {
                "what_it_answers": (
                    "how high the C3 statistic reads when two outcomes share "
                    "only this sample and this expression matrix and share no "
                    "aspect signal at all"
                ),
                "pairs": noise_floor,
            },
            "contrast_axes_through_the_same_pipeline": {
                "what_it_answers": (
                    "what the C3 statistic reads for axes that really are "
                    "different. Fibrosis is the informative one: it is a "
                    "separate histological axis, so its value says whether an "
                    "observed C3 near the threshold means the three aspects are "
                    "alike or merely that any two liver axes read that high "
                    "here, which would make 0.80 unreachable by construction. "
                    "recorded_sex is the orthogonal reference."
                ),
                "pairs": contrast,
            },
            "severity_residualized_sensitivity": {
                "definition": (
                    "each aspect rank residualized on the rank of "
                    "nash_crn_component_sum, then re-associated per gene"
                ),
                "caveat": (
                    "The component sum is the sum of these three aspects, so "
                    "residualizing the parts on their own total forces the "
                    "residuals to be negatively dependent. Negative values here "
                    "are expected arithmetic, not evidence of opposing biology."
                ),
                "pairs": residualized,
            },
        },
        "context_axes_outside_the_criteria": context_correlations,
        "reuse_validation": {
            "rank_agreement_with_masld_bench": rank_agreement,
            "association_agreement_with_masld_bench": reference_check,
        },
        "controls": {
            "lobular_necrosis_excluded_from_the_three_aspects": True,
            "fibrosis_reported_as_a_fourth_axis_never_folded_in": True,
            "nulls_measured_from_realised_vectors": True,
            "seed": arguments.seed,
            "permutations_per_aspect": arguments.permutations,
            "label_permutation_draws": arguments.label_draws,
            "split_half_repeats": arguments.splits,
            "shared_sample_c3_floor_draws": arguments.c3_floor_draws,
            "label_bootstrap_resamples": arguments.bootstrap,
            "association_vector_bootstrap_resamples": arguments.c3_bootstrap,
            "bootstrap_resampling_unit": "participant",
        },
        "claim_boundary": prespec["claim_boundary"],
    }

    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "aspect_separability.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    header = "stable_gene_id\t" + "\t".join(ASPECTS) + "\n"
    lines = [header]
    for index, gene in enumerate(realised_ids):
        lines.append(
            gene
            + "\t"
            + "\t".join(f"{float(vectors[aspect][index]):.10f}" for aspect in ASPECTS)
            + "\n"
        )
    (arguments.output / "gene_association_vectors.tsv").write_text(
        "".join(lines), encoding="utf-8"
    )
    freeze_tree(arguments.output, {
        "artifact_class": "gse267145_aspect_separability_result",
        "prespec_id": prespec["prespec_id"],
        "decision": decision,
        "c1_passed": c1_pass,
        "c2_passed": c2_pass,
        "c3_passed": c3_pass,
        "model_fitted": False,
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)

    print(json.dumps({
        "decision": decision,
        "finding": finding,
        "c1": {"passed": c1_pass,
               "abs_spearman": {f"{a}|{b}": round(entry["abs_spearman"], 4)
                                for entry, (a, b) in zip(c1_pairs, PAIRS)}},
        "c2": {"passed": c2_pass, "participation_ratio": round(ratio, 4)},
        "c3": {"passed": c3_pass,
               "abs_spearman": {f"{a}|{b}": round(entry["abs_spearman"], 4)
                                for entry, (a, b) in zip(c3_pairs, PAIRS)}},
        "realised_gene_universe": int(universe.size),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
