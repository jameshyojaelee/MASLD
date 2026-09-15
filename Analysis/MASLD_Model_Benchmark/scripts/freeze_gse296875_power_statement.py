#!/usr/bin/env python3
"""Freeze the prespecified power statement for the GSE296875 endpoints.

Everything is derived from the frozen design: the outcome vectors and their real
denominators (steatosis 38 of 39, any-fibrosis 37 of 39 at 15 positives and 22
negatives), the donor-grouped folds, the pooled out-of-fold metric, and the
Benjamini-Hochberg family of twelve. **No observed model performance is read.** A
power calculation taken off a realised effect is circular; this one plants an
effect of stated size and measures how often the design detects it.

Reference effect sizes are conventional external benchmarks, chosen before
looking at anything: rho of 0.3 and 0.5 for a rank association, and a true AUROC
of 0.70 and 0.80 for discrimination. They are not derived from this cohort.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from masld_bench.evaluators.detectable_effect import (  # noqa: E402
    BH_SINGLE_TRUE_ALPHA,
    CONFIRMATORY_FAMILY_SIZE,
    NOMINAL_ALPHA,
    TARGET_POWER,
    auprc_critical_values,
    auprc_power_curve,
    minimum_detectable_effect,
    power_verdict,
    required_sample_size,
    required_sample_size_binary,
    spearman_critical_values,
    spearman_power_curve,
)


class PowerStatementError(RuntimeError):
    """Raised when a bound input does not verify."""


NULL_DRAWS = 40_000
REPLICATES = 4_000
SWEEP_REPLICATES = 1_500
SWEEP_NULL_DRAWS = 12_000
SEED = 20260825

RHO_GRID = tuple(round(0.05 * step, 2) for step in range(0, 19))
SEPARATION_GRID = tuple(round(0.15 * step, 3) for step in range(0, 19))
SAMPLE_SIZES = (39, 50, 75, 100, 150, 200, 300, 400, 600, 800)

REFERENCE_RHO = (0.30, 0.50)
REFERENCE_AUROC = (0.70, 0.80)

BOUND = {
    "phenotype_endpoints": (
        "executions/gse296875-phenotype-endpoints-20260825",
        "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    ),
    "donor_folds": (
        "executions/gse296875-donor-folds-20260825",
        "e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416",
    ),
    "campaign_spec": (
        "executions/gse296875-phenotype-campaign-spec-21107179",
        "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    ),
    "null_reference_note_v3": (
        "executions/ranking-metric-null-reference-note-v3-21109019",
        "3f371ba28847f28c5c17295d50dc6ca65f51c4bde6039ee31e2bcf764603d40f",
    ),
    "decision_thresholds": (
        "executions/gse296875-decision-thresholds-v2-21114548",
        "d1343044e312dd88c1026f4f6af1b6f284d20f4f4b8e155c6fb3e45964d956d8",
    ),
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    # The header is taken from the first row, so a heterogeneous list would
    # silently drop or misalign columns.  This project has twice been bitten by
    # positionally ragged tables; fail loudly here instead of writing one.
    fieldnames = tuple(rows[0])
    for index, row in enumerate(rows):
        if tuple(row) != fieldnames:
            raise PowerStatementError(
                f"{path.name}: row {index} has key set {sorted(row)}, "
                f"header is {sorted(fieldnames)}. Write one table per estimand "
                "rather than a ragged union."
            )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def separation_for_auroc(auroc: float) -> float:
    return sqrt(2.0) * NormalDist().inv_cdf(auroc)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--realised", type=Path, nargs="*", default=[])
    arguments = parser.parse_args()
    root = arguments.root

    bound: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in BOUND.items():
        observed = digest(root / relative / "ARTIFACTS.json")
        if observed != expected:
            raise PowerStatementError(f"bound artifact changed: {relative}")
        bound[name] = {"path": relative, "artifacts_sha256": observed}

    endpoint_rows = read_tsv(
        root
        / "executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv"
    )
    steatosis = [
        float(row["steatosis_numeric"])
        for row in endpoint_rows
        if row["steatosis_observed"] == "true"
    ]
    fibrosis = [
        1 if row["fibrosis_any"] == "true" else 0
        for row in endpoint_rows
        if row["fibrosis_observed"] == "true"
    ]
    positives = sum(fibrosis)
    negatives = len(fibrosis) - positives
    if len(steatosis) != 38 or len(fibrosis) != 37 or positives != 15:
        raise PowerStatementError("frozen denominators differ")

    adults = {row["donor_id"] for row in endpoint_rows if row["is_adult"] == "true"}
    adult_steatosis = [
        float(row["steatosis_numeric"])
        for row in endpoint_rows
        if row["steatosis_observed"] == "true" and row["donor_id"] in adults
    ]
    adult_fibrosis = [
        1 if row["fibrosis_any"] == "true" else 0
        for row in endpoint_rows
        if row["fibrosis_observed"] == "true" and row["donor_id"] in adults
    ]

    # Thresholds are read from the frozen derived output file, never recomputed
    # here and never supplied by hand. The Benjamini-Hochberg level is the
    # 99.583rd percentile, which needs far more than the draw count that is
    # adequate for a nominal p95, so resolving it is a separate frozen job.
    derived = json.loads(
        (
            root
            / BOUND["decision_thresholds"][0]
            / "thresholds"
            / "decision_thresholds.json"
        ).read_text()
    )
    if derived["observed_model_performance_read"] is not False:
        raise PowerStatementError("threshold artifact claims to read a result")
    for key, record in derived["thresholds"].items():
        for field in ("label_vector_ordering", "quantile_convention", "seed"):
            if not record.get(field):
                raise PowerStatementError(
                    f"threshold {key} does not record {field}; draw count and "
                    "seed alone are not a provenance record"
                )
    rank_critical = {
        rule: derived["thresholds"]["steatosis_all_donors"][rule]["frozen_value"]
        for rule in ("nominal", "bh_single_true")
    }
    auprc_critical = {
        rule: derived["thresholds"]["fibrosis_all_donors"][rule]["frozen_value"]
        for rule in ("nominal", "bh_single_true")
    }

    rank_ceiling = spearman_power_curve(
        steatosis, RHO_GRID, replicates=REPLICATES, null_draws=NULL_DRAWS, seed=SEED
    )
    auprc_ceiling = auprc_power_curve(
        positives,
        negatives,
        SEPARATION_GRID,
        replicates=REPLICATES,
        null_draws=NULL_DRAWS,
        seed=SEED,
    )
    adult_rank_ceiling = spearman_power_curve(
        adult_steatosis,
        RHO_GRID,
        replicates=REPLICATES,
        null_draws=NULL_DRAWS,
        seed=SEED,
    )
    adult_auprc_ceiling = auprc_power_curve(
        sum(adult_fibrosis),
        len(adult_fibrosis) - sum(adult_fibrosis),
        SEPARATION_GRID,
        replicates=REPLICATES,
        null_draws=NULL_DRAWS,
        seed=SEED,
    )

    sweeps: list[dict[str, object]] = []
    for rho in REFERENCE_RHO:
        for row in required_sample_size(
            steatosis,
            SAMPLE_SIZES,
            rho,
            replicates=SWEEP_REPLICATES,
            null_draws=SWEEP_NULL_DRAWS,
            seed=SEED,
        ):
            sweeps.append({"endpoint": "steatosis", "reference": f"rho_{rho}", **row})
    for auroc in REFERENCE_AUROC:
        separation = separation_for_auroc(auroc)
        for row in required_sample_size_binary(
            positives / len(fibrosis),
            SAMPLE_SIZES,
            separation,
            replicates=SWEEP_REPLICATES,
            null_draws=SWEEP_NULL_DRAWS,
            seed=SEED,
        ):
            sweeps.append({"endpoint": "fibrosis", "reference": f"auroc_{auroc}", **row})

    def smallest_n(endpoint: str, reference: str, rule: str) -> int | None:
        field = "power_nominal" if rule == "nominal" else "power_bh_single_true"
        candidates = [
            int(row["n"])
            for row in sweeps
            if row["endpoint"] == endpoint
            and row["reference"] == reference
            and float(row[field]) >= TARGET_POWER
        ]
        return min(candidates) if candidates else None

    required = {
        f"{endpoint}|{reference}": {
            "nominal": smallest_n(endpoint, reference, "nominal"),
            "bh_single_true": smallest_n(endpoint, reference, "bh_single_true"),
            "cohort_has": 38 if endpoint == "steatosis" else 37,
        }
        for endpoint, references in (
            ("steatosis", [f"rho_{r}" for r in REFERENCE_RHO]),
            ("fibrosis", [f"auroc_{a}" for a in REFERENCE_AUROC]),
        )
        for reference in references
    }

    realised: list[dict[str, object]] = []
    realised_bound: list[dict[str, str]] = []
    realised_scopes: set[str] = set()
    for folder in arguments.realised:
        audit = json.loads((folder / "power" / "audit.json").read_text())
        if audit["observed_model_performance_read"]:
            raise PowerStatementError("a realised-power input read an observed result")
        realised_scopes.add(audit["lineage_scope"])
        realised.extend(read_tsv(folder / "power" / "realised_power.tsv"))
        realised_bound.append(
            {
                "path": str(folder.relative_to(root)),
                "artifacts_sha256": digest(folder / "ARTIFACTS.json"),
            }
        )

    # The prose published beside this curve names the lineage and argues from its
    # depth that the curve is an upper bound across scopes. Derive the scope from
    # the inputs rather than restating it, so a sweep run on a different lineage
    # fails closed instead of being silently mislabelled.
    if realised and len(realised_scopes) != 1:
        raise PowerStatementError(
            "realised-power inputs must share one lineage scope; got "
            f"{sorted(realised_scopes)}"
        )
    realised_scope = realised_scopes.pop() if realised_scopes else None
    if realised and realised_scope != "hepatocyte":
        raise PowerStatementError(
            "the published upper-bound argument is written for hepatocyte, the "
            f"deepest lineage; got {realised_scope!r}. Rewrite the rationale "
            "before publishing a different scope."
        )

    smallest_values = [
        value
        for record in required.values()
        for value in (record["bh_single_true"], record["nominal"])
        if value is not None
    ]
    unreached = [
        key
        for key, record in required.items()
        if record["bh_single_true"] is None
    ]
    verdict = power_verdict(smallest_values, cohort_n=39)

    statement = {
        "schema_version": "masld-bench-gse296875-power-statement-v1",
        "dataset_id": "gse296875",
        "unit_of_inference": "donor",
        "derived_from": "frozen_design_only",
        "observed_model_performance_read": False,
        "circularity_note": (
            "A power calculation computed off a realised effect is circular. "
            "Every number here plants an effect of stated size and measures the "
            "detection rate under the frozen decision rule. No candidate metric, "
            "prediction table, or score artifact is opened by any component."
        ),
        "design": {
            "donors": 39,
            "steatosis_observed": len(steatosis),
            "fibrosis_observed": len(fibrosis),
            "fibrosis_positive": positives,
            "fibrosis_negative": negatives,
            "adult_only_steatosis_observed": len(adult_steatosis),
            "adult_only_fibrosis_observed": len(adult_fibrosis),
            "adult_only_fibrosis_positive": sum(adult_fibrosis),
            "outer_folds": 5,
            "pooled_out_of_fold_metric_only": True,
            "fold_3_has_zero_fibrosis_positives": True,
            "per_fold_metric_averaging_is_a_different_estimand_and_invalid_here": True,
            "leave_one_well_out_undefined_for": ["well3", "well5", "well7"],
            "donor_is_nested_within_well": True,
            "bootstrap_replicates": 10_000,
        },
        "decision_rules": {
            "nominal": {
                "alpha": NOMINAL_ALPHA,
                "role": "upper bound on power; ignores multiplicity",
                "steatosis_critical_absolute_rho": rank_critical["nominal"],
                "fibrosis_critical_auprc": auprc_critical["nominal"],
            },
            "bh_single_true": {
                "alpha": BH_SINGLE_TRUE_ALPHA,
                "family_size": CONFIRMATORY_FAMILY_SIZE,
                "role": (
                    "lower bound on power. Benjamini-Hochberg rejects the "
                    "smallest p only if it clears alpha/m, so a family of twelve "
                    "with exactly one true effect reduces to Bonferroni. If more "
                    "scopes were truly non-null BH would be less stringent, so "
                    "the two rules bracket the operating characteristic instead "
                    "of pretending to one number."
                ),
                "steatosis_critical_absolute_rho": rank_critical["bh_single_true"],
                "fibrosis_critical_auprc": auprc_critical["bh_single_true"],
            },
        },
        "information_ceiling": {
            "meaning": (
                "Power when an oracle predictor of the stated strength is handed "
                "to the evaluator. No model can beat this at this n, because the "
                "limit is the sampling noise of the metric rather than the "
                "quality of the learner. A design that cannot detect an effect "
                "here cannot detect it at all."
            ),
            "steatosis_curve": "ceiling_steatosis.tsv",
            "fibrosis_curve": "ceiling_fibrosis.tsv",
            "adult_only_steatosis_curve": "ceiling_steatosis_adult_only.tsv",
            "adult_only_fibrosis_curve": "ceiling_fibrosis_adult_only.tsv",
            "minimum_detectable_effect_at_80_percent_power": {
                "steatosis_nominal": minimum_detectable_effect(
                    rank_ceiling, rule="nominal"
                ),
                "steatosis_bh_single_true": minimum_detectable_effect(
                    rank_ceiling, rule="bh_single_true"
                ),
                "fibrosis_nominal": minimum_detectable_effect(
                    auprc_ceiling, rule="nominal"
                ),
                "fibrosis_bh_single_true": minimum_detectable_effect(
                    auprc_ceiling, rule="bh_single_true"
                ),
                "adult_only_steatosis_bh_single_true": minimum_detectable_effect(
                    adult_rank_ceiling, rule="bh_single_true"
                ),
                "adult_only_fibrosis_bh_single_true": minimum_detectable_effect(
                    adult_auprc_ceiling, rule="bh_single_true"
                ),
            },
        },
        "realised_power": {
            "meaning": (
                "Power after the predictor has to be learned from 39 donors by "
                "the frozen pipeline, on the real donor-by-lineage expression "
                "matrix so that feature covariance and the n-to-p ratio are "
                "real. Always at or below the ceiling; the gap is the price of "
                "estimation."
            ),
            "curve": "realised_power.tsv" if realised else None,
            "inputs": realised_bound,
            "scope": realised_scope,
            "scope_is_the_favourable_case_not_the_typical_one": (
                "Hepatocyte is the deepest lineage: all 39 donor units clear the "
                "frozen twenty-nucleus threshold, so no block is imputed and no "
                "power is lost to a missing-lineage indicator. Every other scope "
                "is thinner. Cholangiocyte loses 11 of 39 units to the threshold "
                "and t_cell loses 10, and an imputed block contributes no signal "
                "while still costing a degree of freedom. **This curve is "
                "therefore an upper bound on realised power across scopes, not a "
                "typical scope.** Read any other lineage as strictly worse."
            ),
            "planted_variance_fraction_swept": [0.02, 0.10],
            "sweep_rationale": (
                "How much donor-level expression variance a planted biological "
                "program occupies is the one quantity not fixed by the frozen "
                "design. It is swept rather than chosen so its influence is "
                "visible. If the two curves barely differ, that is itself the "
                "answer to the first question a reader will ask."
            ),
        },
        "required_sample_size": {
            "meaning": (
                "Donors needed for 80 percent power at conventional external "
                "effect-size benchmarks, chosen independently of this cohort."
            ),
            "reference_effects": {
                "steatosis": [f"rho_{r}" for r in REFERENCE_RHO],
                "fibrosis": [f"true_auroc_{a}" for a in REFERENCE_AUROC],
            },
            "curve": {
                "steatosis": "required_sample_size_steatosis.tsv",
                "fibrosis": "required_sample_size_fibrosis.tsv",
            },
            "smallest_n_reaching_80_percent": required,
            "cohort_has_donors": 39,
            "comparison_to_the_program_arbitration_finding": {
                "parallel_result": (
                    "A separate analysis in this project found that single-cell "
                    "data cannot arbitrate abundance versus activity for any of "
                    "117 programs, because that question needs roughly 206 to "
                    "483 donors per cell type and the data have 13 to 37."
                ),
                "required_donors_per_cell_type_there": [206, 483],
                "donors_available_there": [13, 37],
                "why_the_comparison_is_worth_making": (
                    "The two results approach the same cohort from opposite "
                    "directions: one from arbitrating programs, one from "
                    "detecting an endpoint. If the required sample sizes land in "
                    "the same range, the design is underpowered by roughly an "
                    "order of magnitude and the two corroborate each other. If "
                    "the endpoint requirement lands near 39, then the design was "
                    "adequately powered and the negative is a statement about "
                    "the biology rather than about the cohort. The table decides "
                    "which; it is not assumed either way."
                ),
                "verdict": verdict,
                "reference_effects_never_reaching_80_percent_by_800_donors": unreached,
            },
        },
        "reference_effects_are_external": (
            "rho of 0.3 and 0.5, and true AUROC of 0.70 and 0.80, are "
            "conventional benchmarks. None is derived from this cohort's data or "
            "from any model's result."
        ),
        "scope_and_claims": {
            "framing": "cross-sectional histology-associated development",
            "forbidden": [
                "MASLD diagnosis",
                "MASH diagnosis",
                "NAS",
                "standard fibrosis stage",
                "longitudinal progression",
                "external validation",
                "champion confirmation",
            ],
            "phenotype_is_evaluator_only_on_every_model_path": True,
            "transforms_fit_inside_donor_training_folds_only": True,
            "hotspot_programs_project_only_after_selection_is_frozen": True,
            "adult_only_is_a_prespecified_mask_no_donor_silently_dropped": True,
            "donors_under_18": 5,
        },
        "bound_artifacts": bound,
        "draws": {
            "null_draws": NULL_DRAWS,
            "replicates_per_grid_point": REPLICATES,
            "sweep_replicates": SWEEP_REPLICATES,
            "seed": SEED,
        },
        "status": "pass_power_statement",
    }

    arguments.output.mkdir(parents=True, exist_ok=False)
    write_tsv(
        arguments.output / "ceiling_steatosis.tsv",
        [point.to_dict() for point in rank_ceiling],
    )
    write_tsv(
        arguments.output / "ceiling_fibrosis.tsv",
        [point.to_dict() for point in auprc_ceiling],
    )
    write_tsv(
        arguments.output / "ceiling_steatosis_adult_only.tsv",
        [point.to_dict() for point in adult_rank_ceiling],
    )
    write_tsv(
        arguments.output / "ceiling_fibrosis_adult_only.tsv",
        [point.to_dict() for point in adult_auprc_ceiling],
    )
    # Steatosis sweeps carry Spearman-rho fields and fibrosis sweeps carry
    # binary separation/AUROC/AUPRC fields.  They are different estimands and
    # do not share a header.
    write_tsv(
        arguments.output / "required_sample_size_steatosis.tsv",
        [row for row in sweeps if row["endpoint"] == "steatosis"],
    )
    write_tsv(
        arguments.output / "required_sample_size_fibrosis.tsv",
        [row for row in sweeps if row["endpoint"] == "fibrosis"],
    )
    if realised:
        write_tsv(arguments.output / "realised_power.tsv", realised)
    (arguments.output / "power_statement.json").write_text(
        json.dumps(statement, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": statement["status"], "required": required}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
