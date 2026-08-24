"""Registered endpoint-recomputed bootstrap evaluators for prospective power.

The power firewall never accepts caller-supplied donor summaries.  It binds one
strict joined development table, dispatches by a frozen evaluator identifier,
and recomputes candidate, baseline, and paired-difference endpoints on the same
biological-unit draws. Unknown evaluators fail closed.
"""

from __future__ import annotations

import csv
from collections import OrderedDict
from dataclasses import dataclass, replace
import json
from math import isfinite
from pathlib import Path
from random import Random
import re
from threading import RLock
from typing import Callable, Mapping, Sequence

from ..artifacts import ArtifactError, reject_symlink_components
from ..hashing import (
    HashingError,
    canonical_json,
    canonical_sha256,
    sha256_bytes,
    sha256_file,
)
from .metrics import (
    MetricError,
    average_precision,
    donor_class_balanced_weights,
    fisher_z_mean,
    multinomial_deviance,
    spearman_correlation,
)


class EndpointPowerError(ValueError):
    """Raised when endpoint evidence cannot reproduce the registered metric."""


@dataclass(frozen=True, slots=True)
class EndpointBootstrapResult:
    evaluator_id: str
    observed_effect: float
    candidate_primary_metric: float
    baseline_primary_metric: float
    paired_difference_bootstrap: tuple[float, ...]
    candidate_primary_bootstrap: tuple[float, ...]
    baseline_primary_bootstrap: tuple[float, ...]
    n_units: int
    independent_unit_counts: Mapping[str, int]
    unit_set_sha256: str
    n_rows: int
    row_set_sha256: str
    strata: tuple[str, ...]
    seed: int
    n_resamples: int

    def to_dict(self) -> dict[str, object]:
        return {
            "evaluator_id": self.evaluator_id,
            "observed_effect": self.observed_effect,
            "candidate_primary_metric": self.candidate_primary_metric,
            "baseline_primary_metric": self.baseline_primary_metric,
            "paired_difference_bootstrap": list(self.paired_difference_bootstrap),
            "candidate_primary_bootstrap": list(self.candidate_primary_bootstrap),
            "baseline_primary_bootstrap": list(self.baseline_primary_bootstrap),
            "n_units": self.n_units,
            "independent_unit_counts": dict(self.independent_unit_counts),
            "unit_set_sha256": self.unit_set_sha256,
            "n_rows": self.n_rows,
            "row_set_sha256": self.row_set_sha256,
            "strata": list(self.strata),
            "seed": self.seed,
            "n_resamples": self.n_resamples,
        }


_HASH = re.compile(r"^[0-9a-f]{64}$")


def _strict_parameters(
    parameters: Mapping[str, object], expected: frozenset[str]
) -> dict[str, object]:
    if not isinstance(parameters, Mapping):
        raise EndpointPowerError("evaluator parameters must be an object")
    if set(parameters) != expected:
        missing = sorted(expected.difference(parameters))
        unknown = sorted(set(parameters).difference(expected))
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if unknown:
            details.append("unknown=" + ",".join(unknown))
        raise EndpointPowerError(
            "evaluator parameters do not match the registered schema: "
            + "; ".join(details)
        )
    return dict(parameters)


def _string_roster(value: object, name: str) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)) or not value:
        raise EndpointPowerError(f"{name} must be a non-empty array")
    if any(not isinstance(item, str) for item in value):
        raise EndpointPowerError(f"{name} must contain only strings")
    result = tuple(value)
    if any(not item or "\x00" in item for item in result):
        raise EndpointPowerError(f"{name} contains an invalid identifier")
    if len(set(result)) != len(result) or tuple(sorted(result)) != result:
        raise EndpointPowerError(f"{name} must be unique and canonically sorted")
    return result


def _read_table(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise EndpointPowerError(f"endpoint table is not a regular file: {path}")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or tuple(reader.fieldnames) != fields:
                raise EndpointPowerError(
                    "endpoint table header must be exactly: " + "\t".join(fields)
                )
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise EndpointPowerError(f"cannot read endpoint table: {error}") from error
    if not rows:
        raise EndpointPowerError("endpoint table must contain at least one row")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise EndpointPowerError(
            "endpoint table rows must contain exactly the registered columns"
        )
    row_hashes = [row["row_hash"] for row in rows]
    if any(not _HASH.fullmatch(value) for value in row_hashes):
        raise EndpointPowerError("row_hash values must be lowercase SHA-256 identifiers")
    if len(set(row_hashes)) != len(row_hashes):
        raise EndpointPowerError("endpoint table row_hash values must be unique")
    return rows


def _hash_ids(kind: str, values: Sequence[str]) -> str:
    return canonical_sha256({kind: sorted(values)})


def _finite(row: Mapping[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise EndpointPowerError(f"{field} must contain finite numeric values") from error
    if not isfinite(value):
        raise EndpointPowerError(f"{field} must contain finite numeric values")
    return value


def _hashed(row: Mapping[str, str], field: str) -> str:
    value = row.get(field, "")
    if not _HASH.fullmatch(value):
        raise EndpointPowerError(f"{field} values must be lowercase SHA-256 identifiers")
    return value


def _resample_primary_metrics(
    *,
    n_resamples: int,
    draw: Callable[[], tuple[float, float]],
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    if not isinstance(n_resamples, int) or isinstance(n_resamples, bool) or n_resamples < 1:
        raise EndpointPowerError("n_resamples must be a positive integer")
    candidate_values: list[float] = []
    baseline_values: list[float] = []
    differences: list[float] = []
    for replicate_index in range(n_resamples):
        try:
            candidate_raw, baseline_raw = draw()
            candidate = float(candidate_raw)
            baseline = float(baseline_raw)
        except (EndpointPowerError, MetricError, ValueError, ZeroDivisionError) as error:
            raise EndpointPowerError(
                "registered endpoint is undefined for fixed bootstrap replicate "
                f"{replicate_index}; draws may not be retried or discarded"
            ) from error
        difference = candidate - baseline
        if not (
            isfinite(candidate) and isfinite(baseline) and isfinite(difference)
        ):
            raise EndpointPowerError(
                "registered endpoint returned a non-finite value for fixed bootstrap "
                f"replicate {replicate_index}"
            )
        candidate_values.append(candidate)
        baseline_values.append(baseline)
        differences.append(difference)
    return tuple(candidate_values), tuple(baseline_values), tuple(differences)


def _fixed_macro_f1(
    observed: Sequence[str],
    predicted: Sequence[str],
    donor_ids: Sequence[str],
    class_roster: Sequence[str],
) -> float:
    weights = donor_class_balanced_weights(donor_ids, observed)
    scores = []
    for label in class_roster:
        true_positive = sum(
            weight
            for truth, guess, weight in zip(observed, predicted, weights, strict=True)
            if truth == label and guess == label
        )
        false_positive = sum(
            weight
            for truth, guess, weight in zip(observed, predicted, weights, strict=True)
            if truth != label and guess == label
        )
        false_negative = sum(
            weight
            for truth, guess, weight in zip(observed, predicted, weights, strict=True)
            if truth == label and guess != label
        )
        denominator = 2.0 * true_positive + false_positive + false_negative
        scores.append(0.0 if denominator == 0.0 else 2.0 * true_positive / denominator)
    return sum(scores) / len(scores)


def _cell_macro_f1(
    path: Path, n_resamples: int, seed: int, parameters: Mapping[str, object]
) -> EndpointBootstrapResult:
    raw = _strict_parameters(parameters, frozenset({"class_roster"}))
    roster = _string_roster(raw["class_roster"], "class_roster")
    rows = _read_table(
        path,
        (
            "row_hash",
            "unit_hash",
            "observed_class",
            "candidate_class",
            "baseline_class",
        ),
    )
    for row in rows:
        _hashed(row, "unit_hash")
        for field in ("observed_class", "candidate_class", "baseline_class"):
            if row[field] not in roster:
                raise EndpointPowerError(f"{field} is outside the frozen class_roster")
    units = sorted({row["unit_hash"] for row in rows})
    if len(units) < 2:
        raise EndpointPowerError("cell endpoint requires at least two donors")
    by_unit: dict[str, list[dict[str, str]]] = {unit: [] for unit in units}
    for row in rows:
        by_unit[row["unit_hash"]].append(row)

    def metric(selected: Sequence[tuple[int, str]]) -> tuple[float, float]:
        sample_rows: list[dict[str, str]] = []
        sample_units: list[str] = []
        for draw_index, unit in selected:
            members = by_unit[unit]
            sample_rows.extend(members)
            sample_units.extend([f"{draw_index}:{unit}"] * len(members))
        observed = [row["observed_class"] for row in sample_rows]
        candidate = [row["candidate_class"] for row in sample_rows]
        baseline = [row["baseline_class"] for row in sample_rows]
        return (
            _fixed_macro_f1(observed, candidate, sample_units, roster),
            _fixed_macro_f1(observed, baseline, sample_units, roster),
        )

    candidate_primary, baseline_primary = metric(tuple(enumerate(units)))
    observed_effect = candidate_primary - baseline_primary
    rng = Random(seed)

    def draw() -> tuple[float, float]:
        selected = [
            (index, units[rng.randrange(len(units))]) for index in range(len(units))
        ]
        return metric(selected)

    candidate_bootstrap, baseline_bootstrap, distribution = (
        _resample_primary_metrics(n_resamples=n_resamples, draw=draw)
    )
    row_ids = [row["row_hash"] for row in rows]
    return EndpointBootstrapResult(
        evaluator_id="cell_donor_balanced_macro_f1_v1",
        observed_effect=observed_effect,
        candidate_primary_metric=candidate_primary,
        baseline_primary_metric=baseline_primary,
        paired_difference_bootstrap=distribution,
        candidate_primary_bootstrap=candidate_bootstrap,
        baseline_primary_bootstrap=baseline_bootstrap,
        n_units=len(units),
        independent_unit_counts={"donor": len(units)},
        unit_set_sha256=_hash_ids("donor_hashes", units),
        n_rows=len(rows),
        row_set_sha256=_hash_ids("row_hashes", row_ids),
        strata=roster,
        seed=seed,
        n_resamples=n_resamples,
    )


def _bulk_spearman(
    path: Path, n_resamples: int, seed: int, parameters: Mapping[str, object]
) -> EndpointBootstrapResult:
    _strict_parameters(parameters, frozenset())
    rows = _read_table(
        path,
        ("row_hash", "unit_hash", "observed", "candidate", "baseline"),
    )
    units = [_hashed(row, "unit_hash") for row in rows]
    if len(set(units)) != len(units):
        raise EndpointPowerError("bulk endpoint requires exactly one row per donor")
    if len(units) < 3:
        raise EndpointPowerError("bulk rank endpoint requires at least three donors")
    observed = [_finite(row, "observed") for row in rows]
    candidate = [_finite(row, "candidate") for row in rows]
    baseline = [_finite(row, "baseline") for row in rows]

    def metric(indices: Sequence[int]) -> tuple[float, float]:
        truth = [observed[index] for index in indices]
        return (
            spearman_correlation(truth, [candidate[index] for index in indices]),
            spearman_correlation(truth, [baseline[index] for index in indices]),
        )

    candidate_primary, baseline_primary = metric(tuple(range(len(rows))))
    observed_effect = candidate_primary - baseline_primary
    rng = Random(seed)

    def draw() -> tuple[float, float]:
        return metric([rng.randrange(len(rows)) for _ in rows])

    candidate_bootstrap, baseline_bootstrap, distribution = (
        _resample_primary_metrics(n_resamples=n_resamples, draw=draw)
    )
    row_ids = [row["row_hash"] for row in rows]
    return EndpointBootstrapResult(
        evaluator_id="bulk_paired_spearman_gain_v1",
        observed_effect=observed_effect,
        candidate_primary_metric=candidate_primary,
        baseline_primary_metric=baseline_primary,
        paired_difference_bootstrap=distribution,
        candidate_primary_bootstrap=candidate_bootstrap,
        baseline_primary_bootstrap=baseline_bootstrap,
        n_units=len(units),
        independent_unit_counts={"donor": len(units)},
        unit_set_sha256=_hash_ids("donor_hashes", units),
        n_rows=len(rows),
        row_set_sha256=_hash_ids("row_hashes", row_ids),
        strata=(),
        seed=seed,
        n_resamples=n_resamples,
    )


def _variant_fisher_z(
    path: Path, n_resamples: int, seed: int, parameters: Mapping[str, object]
) -> EndpointBootstrapResult:
    raw = _strict_parameters(parameters, frozenset({"strata"}))
    strata = _string_roster(raw["strata"], "strata")
    rows = _read_table(
        path,
        (
            "row_hash",
            "unit_hash",
            "block_hash",
            "stratum",
            "observed",
            "candidate",
            "baseline",
        ),
    )
    blocks = sorted({_hashed(row, "block_hash") for row in rows})
    for row in rows:
        _hashed(row, "unit_hash")
        if row["stratum"] not in strata:
            raise EndpointPowerError("variant row is outside the frozen strata roster")
        for field in ("observed", "candidate", "baseline"):
            _finite(row, field)
    if len(blocks) < 2:
        raise EndpointPowerError("variant endpoint requires at least two LD blocks")
    by_block: dict[str, list[dict[str, str]]] = {block: [] for block in blocks}
    for row in rows:
        by_block[row["block_hash"]].append(row)

    def metric(selected_blocks: Sequence[str]) -> tuple[float, float]:
        selected = [row for block in selected_blocks for row in by_block[block]]
        candidate_correlations: list[float] = []
        baseline_correlations: list[float] = []
        for stratum in strata:
            subset = [row for row in selected if row["stratum"] == stratum]
            if len(subset) < 2:
                raise EndpointPowerError("resample lacks a frozen variant stratum")
            truth = [_finite(row, "observed") for row in subset]
            candidate_correlations.append(
                spearman_correlation(truth, [_finite(row, "candidate") for row in subset])
            )
            baseline_correlations.append(
                spearman_correlation(truth, [_finite(row, "baseline") for row in subset])
            )
        return (
            fisher_z_mean(candidate_correlations),
            fisher_z_mean(baseline_correlations),
        )

    candidate_primary, baseline_primary = metric(blocks)
    observed_effect = candidate_primary - baseline_primary
    rng = Random(seed)

    def draw() -> tuple[float, float]:
        return metric([blocks[rng.randrange(len(blocks))] for _ in blocks])

    candidate_bootstrap, baseline_bootstrap, distribution = (
        _resample_primary_metrics(n_resamples=n_resamples, draw=draw)
    )
    row_ids = [row["row_hash"] for row in rows]
    return EndpointBootstrapResult(
        evaluator_id="variant_ld_block_fisher_z_spearman_gain_v1",
        observed_effect=observed_effect,
        candidate_primary_metric=candidate_primary,
        baseline_primary_metric=baseline_primary,
        paired_difference_bootstrap=distribution,
        candidate_primary_bootstrap=candidate_bootstrap,
        baseline_primary_bootstrap=baseline_bootstrap,
        n_units=len(blocks),
        independent_unit_counts={"ld_block": len(blocks)},
        unit_set_sha256=_hash_ids("ld_block_hashes", blocks),
        n_rows=len(rows),
        row_set_sha256=_hash_ids("row_hashes", row_ids),
        strata=strata,
        seed=seed,
        n_resamples=n_resamples,
    )


def _mpra_locus_spearman(
    path: Path, n_resamples: int, seed: int, parameters: Mapping[str, object]
) -> EndpointBootstrapResult:
    """Evaluate locus-held MPRA allelic effects without implying eQTL validity."""

    _strict_parameters(parameters, frozenset())
    rows = _read_table(
        path,
        (
            "row_hash",
            "unit_hash",
            "block_hash",
            "observed",
            "candidate",
            "baseline",
        ),
    )
    blocks = sorted({_hashed(row, "block_hash") for row in rows})
    for row in rows:
        _hashed(row, "unit_hash")
        for field in ("observed", "candidate", "baseline"):
            _finite(row, field)
    if len(blocks) < 3:
        raise EndpointPowerError(
            "MPRA proxy endpoint requires at least three merged loci"
        )
    by_block: dict[str, list[dict[str, str]]] = {block: [] for block in blocks}
    for row in rows:
        by_block[row["block_hash"]].append(row)

    def metric(selected_blocks: Sequence[str]) -> tuple[float, float]:
        selected = [row for block in selected_blocks for row in by_block[block]]
        if len(selected) < 3:
            raise EndpointPowerError(
                "MPRA proxy endpoint requires at least three selected observations"
            )
        truth = [_finite(row, "observed") for row in selected]
        return (
            spearman_correlation(
                truth, [_finite(row, "candidate") for row in selected]
            ),
            spearman_correlation(
                truth, [_finite(row, "baseline") for row in selected]
            ),
        )

    candidate_primary, baseline_primary = metric(blocks)
    observed_effect = candidate_primary - baseline_primary
    rng = Random(seed)

    def draw() -> tuple[float, float]:
        return metric([blocks[rng.randrange(len(blocks))] for _ in blocks])

    candidate_bootstrap, baseline_bootstrap, distribution = (
        _resample_primary_metrics(n_resamples=n_resamples, draw=draw)
    )
    row_ids = [row["row_hash"] for row in rows]
    return EndpointBootstrapResult(
        evaluator_id="locus_heldout_allelic_spearman_v1",
        observed_effect=observed_effect,
        candidate_primary_metric=candidate_primary,
        baseline_primary_metric=baseline_primary,
        paired_difference_bootstrap=distribution,
        candidate_primary_bootstrap=candidate_bootstrap,
        baseline_primary_bootstrap=baseline_bootstrap,
        n_units=len(blocks),
        independent_unit_counts={"merged_locus": len(blocks)},
        unit_set_sha256=_hash_ids("merged_locus_hashes", blocks),
        n_rows=len(rows),
        row_set_sha256=_hash_ids("row_hashes", row_ids),
        strata=(),
        seed=seed,
        n_resamples=n_resamples,
    )


def _graph_auprc(
    path: Path, n_resamples: int, seed: int, parameters: Mapping[str, object]
) -> EndpointBootstrapResult:
    _strict_parameters(parameters, frozenset())
    rows = _read_table(
        path,
        (
            "row_hash",
            "unit_hash",
            "block_hash",
            "observed_binary",
            "candidate",
            "baseline",
        ),
    )
    blocks = sorted({_hashed(row, "block_hash") for row in rows})
    for row in rows:
        _hashed(row, "unit_hash")
        if row["observed_binary"] not in {"0", "1"}:
            raise EndpointPowerError("observed_binary must contain only 0 or 1")
        _finite(row, "candidate")
        _finite(row, "baseline")
    if len(blocks) < 2:
        raise EndpointPowerError("graph endpoint requires at least two held blocks")
    by_block: dict[str, list[dict[str, str]]] = {block: [] for block in blocks}
    for row in rows:
        by_block[row["block_hash"]].append(row)

    def metric(selected_blocks: Sequence[str]) -> tuple[float, float]:
        selected = [row for block in selected_blocks for row in by_block[block]]
        truth = [int(row["observed_binary"]) for row in selected]
        return (
            average_precision(truth, [_finite(row, "candidate") for row in selected]),
            average_precision(truth, [_finite(row, "baseline") for row in selected]),
        )

    candidate_primary, baseline_primary = metric(blocks)
    observed_effect = candidate_primary - baseline_primary
    rng = Random(seed)

    def draw() -> tuple[float, float]:
        return metric([blocks[rng.randrange(len(blocks))] for _ in blocks])

    candidate_bootstrap, baseline_bootstrap, distribution = (
        _resample_primary_metrics(n_resamples=n_resamples, draw=draw)
    )
    row_ids = [row["row_hash"] for row in rows]
    return EndpointBootstrapResult(
        evaluator_id="graph_ld_block_auprc_gain_v1",
        observed_effect=observed_effect,
        candidate_primary_metric=candidate_primary,
        baseline_primary_metric=baseline_primary,
        paired_difference_bootstrap=distribution,
        candidate_primary_bootstrap=candidate_bootstrap,
        baseline_primary_bootstrap=baseline_bootstrap,
        n_units=len(blocks),
        independent_unit_counts={"held_block": len(blocks)},
        unit_set_sha256=_hash_ids("held_block_hashes", blocks),
        n_rows=len(rows),
        row_set_sha256=_hash_ids("row_hashes", row_ids),
        strata=(),
        seed=seed,
        n_resamples=n_resamples,
    )


def _rna_atac_deviance(
    path: Path, n_resamples: int, seed: int, parameters: Mapping[str, object]
) -> EndpointBootstrapResult:
    raw = _strict_parameters(parameters, frozenset({"strata"}))
    strata = _string_roster(raw["strata"], "strata")
    rows = _read_table(
        path,
        (
            "row_hash",
            "donor_hash",
            "block_hash",
            "stratum",
            "observed",
            "candidate",
            "baseline",
        ),
    )
    donors = sorted({_hashed(row, "donor_hash") for row in rows})
    blocks = sorted({_hashed(row, "block_hash") for row in rows})
    for row in rows:
        if row["stratum"] not in strata:
            raise EndpointPowerError("RNA-ATAC row is outside the frozen strata roster")
        for field in ("observed", "candidate", "baseline"):
            if _finite(row, field) < 0.0:
                raise EndpointPowerError("RNA-ATAC profiles must be nonnegative")
    if len(donors) < 2 or len(blocks) < 2:
        raise EndpointPowerError(
            "RNA-ATAC endpoint requires at least two donors and two genomic blocks"
        )
    by_profile: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (row["donor_hash"], row["block_hash"], row["stratum"])
        by_profile.setdefault(key, []).append(row)

    missing_strata = [
        stratum
        for stratum in strata
        if not any(key[2] == stratum for key in by_profile)
    ]
    if missing_strata:
        raise EndpointPowerError(
            "RNA-ATAC table lacks frozen strata: " + ", ".join(missing_strata)
        )

    def metric(
        selected_donors: Sequence[str], selected_blocks: Sequence[str]
    ) -> tuple[float, float]:
        reductions_by_stratum: dict[str, list[float]] = {
            stratum: [] for stratum in strata
        }
        # Each observed donor x lineage x genomic-block profile is one macro
        # unit. Sparse topology stays sparse; absent profiles are never
        # fabricated as biological zeros. Duplicate draws retain bootstrap
        # multiplicity because the loops traverse the sampled unit vectors.
        for donor in selected_donors:
            for block in selected_blocks:
                for stratum in strata:
                    subset = by_profile.get((donor, block, stratum), ())
                    if not subset:
                        continue
                    observed = [_finite(row, "observed") for row in subset]
                    candidate = [_finite(row, "candidate") for row in subset]
                    baseline = [_finite(row, "baseline") for row in subset]
                    candidate_deviance = multinomial_deviance(observed, candidate)
                    baseline_deviance = multinomial_deviance(observed, baseline)
                    if not isfinite(baseline_deviance) or baseline_deviance <= 0.0:
                        raise EndpointPowerError(
                            "baseline profile deviance must be finite and positive"
                        )
                    reductions_by_stratum[stratum].append(
                        (baseline_deviance - candidate_deviance) / baseline_deviance
                    )
        lineage_scores: list[float] = []
        for stratum in strata:
            reductions = reductions_by_stratum[stratum]
            if not reductions:
                raise EndpointPowerError("resample lacks a frozen RNA-ATAC stratum")
            lineage_scores.append(sum(reductions) / len(reductions))
        return sum(lineage_scores) / len(lineage_scores), 0.0

    candidate_primary, baseline_primary = metric(donors, blocks)
    observed_effect = candidate_primary - baseline_primary
    rng = Random(seed)

    def draw() -> tuple[float, float]:
        sampled_donors = [donors[rng.randrange(len(donors))] for _ in donors]
        sampled_blocks = [blocks[rng.randrange(len(blocks))] for _ in blocks]
        return metric(sampled_donors, sampled_blocks)

    candidate_bootstrap, baseline_bootstrap, distribution = (
        _resample_primary_metrics(n_resamples=n_resamples, draw=draw)
    )
    row_ids = [row["row_hash"] for row in rows]
    unit_identity = {"donor_hashes": donors, "genomic_block_hashes": blocks}
    return EndpointBootstrapResult(
        evaluator_id="rna_atac_two_way_deviance_reduction_v1",
        observed_effect=observed_effect,
        candidate_primary_metric=candidate_primary,
        baseline_primary_metric=baseline_primary,
        paired_difference_bootstrap=distribution,
        candidate_primary_bootstrap=candidate_bootstrap,
        baseline_primary_bootstrap=baseline_bootstrap,
        n_units=min(len(donors), len(blocks)),
        independent_unit_counts={"donor": len(donors), "genomic_block": len(blocks)},
        unit_set_sha256=canonical_sha256(unit_identity),
        n_rows=len(rows),
        row_set_sha256=_hash_ids("row_hashes", row_ids),
        strata=strata,
        seed=seed,
        n_resamples=n_resamples,
    )


_EVALUATORS: Mapping[
    str,
    Callable[[Path, int, int, Mapping[str, object]], EndpointBootstrapResult],
] = {
    "bulk_paired_spearman_gain_v1": _bulk_spearman,
    "cell_donor_balanced_macro_f1_v1": _cell_macro_f1,
    "graph_ld_block_auprc_gain_v1": _graph_auprc,
    "locus_heldout_allelic_spearman_v1": _mpra_locus_spearman,
    "rna_atac_two_way_deviance_reduction_v1": _rna_atac_deviance,
    "variant_ld_block_fisher_z_spearman_gain_v1": _variant_fisher_z,
}

_CACHE_MAX_ENTRIES = 64
_CACHE_LOCK = RLock()
_ENDPOINT_BOOTSTRAP_CACHE: OrderedDict[
    tuple[str, str, str, int, int], EndpointBootstrapResult
] = OrderedDict()


def _copy_result(result: EndpointBootstrapResult) -> EndpointBootstrapResult:
    return replace(
        result,
        independent_unit_counts=dict(result.independent_unit_counts),
    )


def _stat_identity(path: Path) -> tuple[int, int, int, int, int, int, int]:
    status = path.stat()
    return (
        status.st_dev,
        status.st_ino,
        status.st_mode,
        status.st_size,
        status.st_mtime_ns,
        status.st_ctime_ns,
        status.st_nlink,
    )


def _stable_table_snapshot(
    table_path: str | Path,
) -> tuple[Path, str, tuple[int, int, int, int, int, int, int]]:
    candidate = Path(table_path)
    if ".." in candidate.parts:
        raise EndpointPowerError("endpoint table path may not contain '..'")
    try:
        lexical = reject_symlink_components(candidate, label="endpoint table")
        if not lexical.is_file():
            raise EndpointPowerError(
                f"endpoint table is not a regular file: {lexical}"
            )
        resolved = lexical.resolve(strict=True)
        before = _stat_identity(lexical)
        digest = sha256_file(lexical)
        after = _stat_identity(lexical)
        checked = reject_symlink_components(lexical, label="endpoint table")
        if (
            not checked.is_file()
            or checked.resolve(strict=True) != resolved
            or before != after
        ):
            raise EndpointPowerError(
                "endpoint table changed while its content identity was computed"
            )
    except EndpointPowerError:
        raise
    except (ArtifactError, OSError, ValueError) as error:
        raise EndpointPowerError(f"cannot bind endpoint table: {error}") from error
    return resolved, digest, after


def _clear_endpoint_cache_for_testing() -> None:
    with _CACHE_LOCK:
        _ENDPOINT_BOOTSTRAP_CACHE.clear()


def registered_endpoint_evaluators() -> tuple[str, ...]:
    return tuple(sorted(_EVALUATORS))


def recompute_endpoint_bootstrap(
    *,
    evaluator_id: str,
    table_path: str | Path,
    n_resamples: int,
    seed: int,
    parameters: Mapping[str, object],
) -> EndpointBootstrapResult:
    """Recompute one registered paired endpoint from strict development rows."""

    if evaluator_id not in _EVALUATORS:
        raise EndpointPowerError(f"unregistered endpoint evaluator: {evaluator_id}")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise EndpointPowerError("seed must be an integer")
    if (
        not isinstance(n_resamples, int)
        or isinstance(n_resamples, bool)
        or n_resamples < 1
    ):
        raise EndpointPowerError("n_resamples must be a positive integer")
    try:
        serialized_parameters = canonical_json(parameters)
        copied_parameters = json.loads(serialized_parameters)
    except (HashingError, json.JSONDecodeError) as error:
        raise EndpointPowerError(
            f"evaluator parameters are not canonical JSON: {error}"
        ) from error
    if not isinstance(copied_parameters, dict):
        raise EndpointPowerError("evaluator parameters must be an object")

    path, table_sha256, table_identity = _stable_table_snapshot(table_path)
    parameter_sha256 = sha256_bytes(serialized_parameters.encode("utf-8"))
    cache_key = (
        table_sha256,
        evaluator_id,
        parameter_sha256,
        seed,
        n_resamples,
    )
    with _CACHE_LOCK:
        cached = _ENDPOINT_BOOTSTRAP_CACHE.get(cache_key)
        if cached is not None:
            _ENDPOINT_BOOTSTRAP_CACHE.move_to_end(cache_key)
            cached_copy = _copy_result(cached)
        else:
            cached_copy = None
    if cached_copy is not None:
        checked_path, checked_sha256, checked_identity = _stable_table_snapshot(path)
        if (
            checked_path != path
            or checked_sha256 != table_sha256
            or checked_identity != table_identity
        ):
            raise EndpointPowerError(
                "endpoint table changed while a cached result was verified"
            )
        return cached_copy

    result = _EVALUATORS[evaluator_id](
        path, n_resamples, seed, copied_parameters
    )
    checked_path, checked_sha256, checked_identity = _stable_table_snapshot(path)
    if (
        checked_path != path
        or checked_sha256 != table_sha256
        or checked_identity != table_identity
    ):
        raise EndpointPowerError(
            "endpoint table changed while the registered evaluator was running"
        )
    stored = _copy_result(result)
    with _CACHE_LOCK:
        _ENDPOINT_BOOTSTRAP_CACHE[cache_key] = stored
        _ENDPOINT_BOOTSTRAP_CACHE.move_to_end(cache_key)
        while len(_ENDPOINT_BOOTSTRAP_CACHE) > _CACHE_MAX_ENTRIES:
            _ENDPOINT_BOOTSTRAP_CACHE.popitem(last=False)
    return _copy_result(stored)


__all__ = [
    "EndpointBootstrapResult",
    "EndpointPowerError",
    "recompute_endpoint_bootstrap",
    "registered_endpoint_evaluators",
]
