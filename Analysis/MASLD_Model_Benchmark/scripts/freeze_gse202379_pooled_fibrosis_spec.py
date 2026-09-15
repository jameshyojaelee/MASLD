#!/usr/bin/env python3
"""Freeze the prespecified cohort-stratified fibrosis spec before any fitting.

This is a new test on a question that has already returned a negative, so
without a frozen spec it is post hoc by construction. Endpoints, folds,
multiplicity family, decision rule and power are fixed here, before GSE202379
molecular data is opened.

Provenance of the graded set. The 40-donor SAF grading below was derived here,
from the primary GEO series matrix, by this script. It was not handed down and
corrected: the parse builds a per-sample dictionary by splitting every cell on
its own ``key:`` prefix, because the ``!Sample_characteristics_ch1`` rows are
positionally ragged and keying a row by its first cell silently mixes fields.
A parser that makes that mistake loses two real gradings.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from math import atanh, sqrt, tanh
from pathlib import Path
from random import Random
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from masld_bench.evaluators.detectable_effect import (  # noqa: E402
    BH_SINGLE_TRUE_ALPHA,
    NOMINAL_ALPHA,
    TARGET_POWER,
    _quantile,
    _standardised_ranks,
)
from masld_bench.evaluators.metrics import spearman_correlation  # noqa: E402


class SpecError(RuntimeError):
    """Raised when a source or structural expectation does not hold."""


SAF_TRIPLE = re.compile(r"^S(\d+)A(\d+)F(\d+)$")
SERIES_MATRIX = (
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
    "SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/"
    "GSE202379_series_matrix.txt.gz"
)
NULL_DRAWS = 200_000
REPLICATES = 4_000
SEED = 20260825
RHO_GRID = tuple(round(0.05 * step, 2) for step in range(0, 19))

BOUND = {
    "phenotype_endpoints": (
        "executions/gse296875-phenotype-endpoints-20260825",
        "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    ),
    "campaign_spec": (
        "executions/gse296875-phenotype-campaign-spec-21107179",
        "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    ),
    "null_reference_note_v3": (
        "executions/ranking-metric-null-reference-note-v3-21109019",
        "3f371ba28847f28c5c17295d50dc6ca65f51c4bde6039ee31e2bcf764603d40f",
    ),
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def parse_saf(path: Path) -> tuple[dict[str, int], dict[str, str], dict[str, str]]:
    """Per-donor SAF fibrosis grade from the primary series matrix.

    Two traps are handled explicitly. The characteristics rows are positionally
    ragged, so every cell is split on its own key. And the ``saf score`` field is
    never absent: for nineteen samples its value is the literal string
    ``end stage`` or ``healthy control`` where a triple belongs, so testing for
    a missing field finds nothing missing and treats status words as data.
    """

    accessions: list[str] = []
    rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                accessions = [
                    cell.strip().strip('"') for cell in line.split("\t")[1:]
                ]
            elif line.startswith("!Sample_characteristics_ch1"):
                rows.append(
                    [
                        cell.strip().strip('"')
                        for cell in line.rstrip("\n").split("\t")[1:]
                    ]
                )
    if len(accessions) != 59:
        raise SpecError(f"expected 59 GSMs, found {len(accessions)}")
    records: list[dict[str, str]] = [dict() for _ in accessions]
    for row in rows:
        if len(row) != len(accessions):
            raise SpecError("characteristics row width differs from the GSM axis")
        for index, cell in enumerate(row):
            if ":" not in cell:
                continue
            key, value = cell.split(":", 1)
            records[index][key.strip().lower()] = value.strip()

    graded: dict[str, int] = {}
    status: dict[str, str] = {}
    ungraded_value: dict[str, str] = {}
    for record in records:
        donor = record.get("patient id")
        if not donor:
            raise SpecError("a sample carries no patient id")
        status[donor] = record.get("disease status", "")
        raw = record.get("saf score")
        if raw is None:
            raise SpecError("the saf score field is absent; the parse assumption changed")
        match = SAF_TRIPLE.match(raw.replace(" ", ""))
        if match:
            stage = int(match.group(3))
            if donor in graded and graded[donor] != stage:
                raise SpecError(f"conflicting SAF grades within donor {donor}")
            graded[donor] = stage
        else:
            ungraded_value[donor] = raw
    return graded, status, ungraded_value


def fisher_meta(values: list[tuple[float, int]]) -> float:
    """Inverse-variance fixed-effect combination of Fisher-z Spearman."""

    numerator = 0.0
    denominator = 0.0
    for correlation, n in values:
        clipped = max(min(correlation, 0.999999), -0.999999)
        weight = max(n - 3, 1)
        numerator += weight * atanh(clipped)
        denominator += weight
    return numerator / denominator


def power_curve(
    cohorts: list[tuple[str, list[float]]],
    *,
    replicates: int,
    null_draws: int,
    seed: int,
) -> tuple[list[dict[str, object]], dict[str, float]]:
    """Detectable-effect curve for the meta-analysed rank association."""

    references = {name: _standardised_ranks(values) for name, values in cohorts}
    rng = Random(seed)
    null: list[float] = []
    for _ in range(null_draws):
        combined = fisher_meta(
            [
                (
                    spearman_correlation(
                        [rng.gauss(0.0, 1.0) for _ in values], values
                    ),
                    len(values),
                )
                for name, values in cohorts
            ]
        )
        null.append(abs(combined))
    null.sort()
    critical = {
        "nominal": _quantile(null, 1.0 - NOMINAL_ALPHA),
        "bh_single_true": _quantile(null, 1.0 - BH_SINGLE_TRUE_ALPHA),
    }
    points: list[dict[str, object]] = []
    for latent in RHO_GRID:
        residual = sqrt(max(0.0, 1.0 - latent * latent))
        hits = {"nominal": 0, "bh_single_true": 0}
        observed: list[float] = []
        for _ in range(replicates):
            per_cohort = []
            for name, values in cohorts:
                anchor = references[name]
                predictor = [
                    latent * a + residual * rng.gauss(0.0, 1.0) for a in anchor
                ]
                per_cohort.append(
                    (spearman_correlation(predictor, values), len(values))
                )
            combined = abs(fisher_meta(per_cohort))
            observed.append(combined)
            for rule, value in critical.items():
                if combined >= value:
                    hits[rule] += 1
        points.append(
            {
                "planted_rho": latent,
                "mean_meta_z": sum(observed) / len(observed),
                "mean_meta_rho": tanh(sum(observed) / len(observed)),
                "power_nominal": hits["nominal"] / replicates,
                "power_bh_single_true": hits["bh_single_true"] / replicates,
            }
        )
    return points, critical


def minimum_detectable(points: list[dict[str, object]], rule: str) -> float | None:
    field = f"power_{rule}"
    previous = None
    for point in points:
        if float(point[field]) >= TARGET_POWER:
            if previous is None:
                return float(point["planted_rho"])
            span = float(point[field]) - float(previous[field])
            weight = 0.0 if span == 0 else (TARGET_POWER - float(previous[field])) / span
            return float(previous["planted_rho"]) + weight * (
                float(point["planted_rho"]) - float(previous["planted_rho"])
            )
        previous = point
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    root = arguments.root

    bound: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in BOUND.items():
        observed = digest(root / relative / "ARTIFACTS.json")
        if observed != expected:
            raise SpecError(f"bound artifact changed: {relative}")
        bound[name] = {"path": relative, "artifacts_sha256": observed}

    series = Path(SERIES_MATRIX)
    graded, status, ungraded_value = parse_saf(series)
    donors = sorted(set(status))
    ungraded = sorted(set(donors) - set(graded))
    if len(donors) != 47 or len(graded) != 40 or len(ungraded) != 7:
        raise SpecError(
            f"SAF census differs: {len(donors)} donors, {len(graded)} graded"
        )
    distribution = {f"F{stage}": 0 for stage in range(5)}
    for stage in graded.values():
        distribution[f"F{stage}"] += 1
    if distribution != {"F0": 3, "F1": 9, "F2": 12, "F3": 12, "F4": 4}:
        raise SpecError(f"SAF distribution differs: {distribution}")

    endpoint_rows = list(
        csv.DictReader(
            (
                root
                / "executions/gse296875-phenotype-endpoints-20260825/endpoint_lock"
                / "donor_endpoints.tsv"
            ).open(encoding="utf-8", newline=""),
            delimiter="\t",
        )
    )
    binary = [
        1.0 if row["fibrosis_any"] == "true" else 0.0
        for row in endpoint_rows
        if row["fibrosis_observed"] == "true"
    ]
    ordinal = [float(graded[donor]) for donor in sorted(graded)]

    solo_202379, crit_202379 = power_curve(
        [("gse202379", ordinal)],
        replicates=REPLICATES,
        null_draws=NULL_DRAWS,
        seed=SEED,
    )
    solo_296875, crit_296875 = power_curve(
        [("gse296875", binary)],
        replicates=REPLICATES,
        null_draws=NULL_DRAWS,
        seed=SEED + 1,
    )
    meta, crit_meta = power_curve(
        [("gse202379", ordinal), ("gse296875", binary)],
        replicates=REPLICATES,
        null_draws=NULL_DRAWS,
        seed=SEED + 2,
    )

    mde = {
        "gse202379_alone": {
            rule: minimum_detectable(solo_202379, rule)
            for rule in ("nominal", "bh_single_true")
        },
        "gse296875_alone": {
            rule: minimum_detectable(solo_296875, rule)
            for rule in ("nominal", "bh_single_true")
        },
        "cohort_stratified_meta": {
            rule: minimum_detectable(meta, rule)
            for rule in ("nominal", "bh_single_true")
        },
    }

    control_arm = {
        donor: graded.get(donor, "ungraded")
        for donor in donors
        if status[donor] == "Healthy control"
    }

    spec = {
        "schema_version": "masld-bench-gse202379-pooled-fibrosis-spec-v1",
        "frozen_before_any_gse202379_molecular_data_was_opened": True,
        "why_a_frozen_spec": (
            "This is a new test on a question that has already returned a "
            "negative. Without a spec frozen in advance it is post hoc by "
            "construction."
        ),
        "graded_set": {
            "source": "primary GEO series matrix saf score field",
            "series_matrix": SERIES_MATRIX,
            "series_matrix_sha256": digest(series),
            "derived_by": "this script, independently, not adopted from a lane",
            "donors": len(donors),
            "gsms": 59,
            "graded_donors": len(graded),
            "distribution": distribution,
            "ungraded_donors": len(ungraded),
            "ungraded": {donor: ungraded_value.get(donor, "") for donor in ungraded},
            "within_donor_conflicts": 0,
        },
        "parsing_rules_that_are_load_bearing": {
            "ragged_characteristics_rows": (
                "!Sample_characteristics_ch1 rows are positionally ragged: one "
                "row carries gender for some columns and liver lobe for others. "
                "Every cell is split on its own key: prefix. Keying a row by its "
                "first cell mixes fields and loses two real gradings."
            ),
            "saf_score_is_never_absent": (
                "The field is present on all 59 samples. For nineteen of them "
                "the value is the literal string end stage or healthy control "
                "where a triple belongs. A missingness test finds nothing "
                "missing and treats status words as data. Ungraded is defined "
                "as failing to match ^S\\d+A\\d+F\\d+$, never as an absent field."
            ),
        },
        "donor_fstage_documented_tsv_is_not_a_stage_source": {
            "verified_against_the_primary_source": True,
            "agree": 38,
            "overridden": {"P98": {"source": "F1", "table": "F0"}},
            "invented": {
                "end_stage_to_F4": ["PCL16", "PCL17", "PCL18", "PCL103", "PCL104"],
                "healthy_control_to_F0": ["PHL1", "PHL2"],
            },
            "dropped": {"P70": "source F3, absent from the table"},
            "table_total": 46,
            "prohibition": (
                "Not used as a stage source in this spec. All seven invented "
                "rows carry series-matrix provenance for values that file does "
                "not contain, and the invented cells sit at F0 and F4, the "
                "extremes that drive an ordinal correlation."
            ),
            "keyed_on_runs_not_donors": True,
        },
        "endpoints": {
            "gse202379": {
                "endpoint": "source_native_SAF_fibrosis_stage_F0_to_F4",
                "scale": "ordinal, five levels",
                "n": len(ordinal),
                "metric": "cross_fitted_spearman",
                "retires_the_zero_positive_fold_problem": True,
            },
            "gse296875": {
                "endpoint": "source_described_any_fibrosis",
                "scale": "binary, two levels",
                "n": len(binary),
                "positives": int(sum(binary)),
                "negatives": int(len(binary) - sum(binary)),
                "metric": "cross_fitted_spearman_on_the_native_binary_scale",
                "kept_as_its_own_stratum": True,
            },
            "common_scale_rationale": (
                "Both cohorts are scored by Spearman on their own native scale, "
                "so no cross-cohort recode is needed to combine them. The binary "
                "endpoint is a two-level ordinal; it is lower resolution, not a "
                "different quantity."
            ),
        },
        "harmonisation": {
            "silent_recode_prohibited": True,
            "required_columns_on_every_row": [
                "native_scale_value",
                "native_scale_name",
                "f_stage_source",
                "cross_cohort_comparability",
            ],
            "f_stage_source_values": ["saf_graded", "imputed_from_disease_status"],
            "primary_set": "saf_graded only, 40 donors",
            "sensitivity_set": (
                "47 donors including the imputed 7, labelled, never pooled "
                "silently into the primary"
            ),
            "forty_six_is_never_used": True,
            "precedent": (
                "GSE213621's harmonised stage silently pools F0F1 and F3F4 with "
                "no level 4, so identically named levels meant different "
                "clinical things across cohorts. Native scale travels with every "
                "value here for that reason."
            ),
        },
        "combination_rule": {
            "method": "cohort_stratified_fixed_effect_meta_on_fisher_z_spearman",
            "naive_pooling_prohibited": True,
            "why": (
                "GSE202379 is severity-skewed biopsy and transplant tissue; "
                "GSE296875 is a mild cardiometabolic donor cohort aged 13 to 75. "
                "Cohort is partly collinear with severity, so a pooled estimate "
                "would be confounded by construction."
            ),
            "leave_one_cohort_out": True,
            "weights": "n minus three per cohort",
        },
        "unit_of_inference": "donor",
        "donors_not_samples": (
            "47 donors across 59 GSMs. Sample-level rows are never an "
            "inferential unit."
        ),
        "stated_weaknesses": {
            "control_arm": {
                "healthy_control_donors": control_arm,
                "graded_controls": sum(
                    1 for value in control_arm.values() if value != "ungraded"
                ),
                "statement": (
                    "GSE202379 has four healthy control donors and only two "
                    "carry a SAF grade. It cannot carry a control-versus-disease "
                    "contrast. It can carry graded severity, which GSE296875 "
                    "cannot. They are complementary and neither alone is "
                    "sufficient."
                ),
            },
            "f0_is_thin": {
                "n": distribution["F0"],
                "donors": sorted(d for d in graded if graded[d] == 0),
                "statement": (
                    "F0 is three donors, two NAFLD and one healthy control. This "
                    "is stated here rather than discovered in a fold. P30 and "
                    "P98 are the only near-normal donors carrying documented "
                    "gradings, which is why a parse that dropped them would have "
                    "mattered rather than merely being untidy."
                ),
            },
            "no_end_stage_donor_is_graded": (
                "All five end-stage donors are ungraded, so the graded set has "
                "no explant tissue and F4 rests on four NASH-with-cirrhosis "
                "donors."
            ),
        },
        "power": {
            "derived_from_the_frozen_design_not_from_any_result": True,
            "decision_rules": {
                "nominal": {"alpha": NOMINAL_ALPHA, "critical_absolute_meta_z": crit_meta["nominal"]},
                "bh_single_true": {
                    "alpha": BH_SINGLE_TRUE_ALPHA,
                    "critical_absolute_meta_z": crit_meta["bh_single_true"],
                },
            },
            "minimum_detectable_rho_at_80_percent_power": mde,
            "curves": {
                "gse202379_alone": "power_gse202379_alone.tsv",
                "gse296875_alone": "power_gse296875_alone.tsv",
                "cohort_stratified_meta": "power_cohort_stratified_meta.tsv",
            },
            "combined_donors": len(ordinal) + len(binary),
            "statement": (
                "Adding GSE202379 improves the detectable effect but does not "
                "reach the regime this estimand needs. The improvement and its "
                "insufficiency are both stated here, before any pooled result "
                "exists."
            ),
        },
        "two_readings_that_the_design_cannot_separate": (
            "The GSE296875 fibrosis average precisions sat below the continuous "
            "random-scorer null mean. One reading is a real null. The other is "
            "noise, since that null's 95th percentile is 0.612 and an observed "
            "0.42 sits inside one standard deviation of it. Both survive the "
            "data, and the honest statement is that the design could not "
            "distinguish them. This sentence is frozen before the pooled result "
            "exists and must read identically whichever way that result lands."
        ),
        "realised_null_correction": (
            "An earlier report of that campaign described the observed average "
            "precisions as below the random-scorer null mean of 0.458 and "
            "therefore below chance. That overstated the case in the direction "
            "of the conclusion. Simulation of the frozen pipeline at a planted "
            "effect of exactly zero gives a mean realised average precision of "
            "0.415 to 0.417, so the continuous random-scorer null is the wrong "
            "centre for a cross-fitted pipeline whose predictions are noisier "
            "and more tied. The correct statement is that the observed 0.315 to "
            "0.421 sits almost exactly on this pipeline's own null."
        ),
        "firewall": {
            "phenotype_is_evaluator_only_on_every_model_path": True,
            "transforms_fit_inside_donor_training_folds_only": True,
            "predictions_frozen_and_hash_committed_before_any_outcome_join": True,
            "hotspot_programs_project_only_after_selection_is_frozen": True,
        },
        "claims": {
            "allowed": "cross-sectional histology-associated development across two donor cohorts",
            "forbidden": [
                "MASLD diagnosis",
                "MASH diagnosis",
                "NAS",
                "standard fibrosis stage",
                "longitudinal progression",
                "external validation",
                "champion confirmation",
            ],
        },
        "bound_artifacts": bound,
        "draws": {"null_draws": NULL_DRAWS, "replicates": REPLICATES, "seed": SEED},
        "status": "pass_pooled_fibrosis_spec",
    }

    arguments.output.mkdir(parents=True, exist_ok=False)
    for name, rows in (
        ("power_gse202379_alone.tsv", solo_202379),
        ("power_gse296875_alone.tsv", solo_296875),
        ("power_cohort_stratified_meta.tsv", meta),
    ):
        with (arguments.output / name).open("x", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    with (arguments.output / "graded_donors.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=("donor_id", "disease_status", "saf_fibrosis_stage", "f_stage_source"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for donor in sorted(donors):
            writer.writerow(
                {
                    "donor_id": donor,
                    "disease_status": status[donor],
                    "saf_fibrosis_stage": graded.get(donor, ""),
                    "f_stage_source": "saf_graded" if donor in graded else "ungraded",
                }
            )
    (arguments.output / "pooled_fibrosis_spec.json").write_text(
        json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": spec["status"], "mde": mde}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
