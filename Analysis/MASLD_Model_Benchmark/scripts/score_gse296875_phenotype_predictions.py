#!/usr/bin/env python3
"""Score frozen GSE296875 phenotype predictions, pooled and once per endpoint.

This process runs after the predictions are frozen and hash-committed.  It is
the only place a prediction is joined to an outcome.  Every metric here is
computed once over a complete registered roster; the enforcement lives in
``masld_bench.gse296875_phenotype_scoring`` and is exercised by its unit tests.

Two references accompany every number, because neither is optional at this
sample size.  For fibrosis the comparison is the continuous random-score
distribution, not the prevalence line, since average precision is upward biased
with fifteen positives in thirty-seven donors.  For steatosis it is the
correlation a random scorer reaches five percent of the time.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from masld_bench.evaluators.stats import benjamini_hochberg  # noqa: E402
from masld_bench.gse296875_phenotype_scoring import (  # noqa: E402
    METRIC_BY_ENDPOINT,
    PRIMARY_LINEAGES_DEFAULT,
    PooledMetricViolation,
    confirmatory_family,
    endpoint_values,
    load_predictions,
    permutation_null,
    pooled_donor_bootstrap,
    pooled_metric,
    random_score_reference,
    read_tsv,
    registered_roster,
)


class ScoringError(RuntimeError):
    """Raised when a frozen input does not verify."""


ARMS = ("molecular", "metadata", "molecular_metadata")
SECONDARY_LINEAGES = ("endothelial_cell", "b_cell")
BOOTSTRAP_REPLICATES = 10_000
SEED = 20260825

BOUND = {
    "phenotype_endpoints": "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    "donor_folds": "e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416",
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def scoped(predictions: dict[str, float], roster) -> dict[str, float]:
    """Restrict a prediction map to a registered roster, never the reverse."""

    return {
        donor: value for donor, value in predictions.items() if donor in roster.donor_ids
    }


def score_one(
    table: Path,
    arm: str,
    scope_name: str,
    endpoint: str,
    roster,
    outcomes: dict[str, float],
    baseline_table: Path | None = None,
    baseline_arm: str | None = None,
) -> dict[str, object]:
    kind = METRIC_BY_ENDPOINT[endpoint]
    predictions = scoped(
        load_predictions(table, arm=arm, lineage_scope=scope_name, endpoint=endpoint),
        roster,
    )
    record: dict[str, object] = {
        "arm": arm,
        "lineage_scope": scope_name,
        "endpoint": endpoint,
        "metric": kind,
        "scope": roster.scope,
        "donors": len(roster.donor_ids),
    }
    try:
        record["estimate"] = pooled_metric(kind, predictions, outcomes, roster)
    except PooledMetricViolation as error:
        record["estimate"] = None
        record["undefined_reason"] = str(error)
        return record

    bootstrap = pooled_donor_bootstrap(
        kind,
        predictions,
        outcomes,
        roster,
        n_resamples=BOOTSTRAP_REPLICATES,
        seed=SEED,
    )
    record["bootstrap"] = bootstrap.to_dict()
    null = permutation_null(
        kind, predictions, outcomes, roster, n_permutations=BOOTSTRAP_REPLICATES, seed=SEED
    )
    record["label_permutation_null"] = null.to_dict()
    record["permutation_p_value"] = null.p_value
    reference = random_score_reference(
        kind, outcomes, roster, n_draws=BOOTSTRAP_REPLICATES, seed=SEED
    )
    record["continuous_random_score_reference"] = reference.to_dict()
    record["exceeds_random_score_95th_percentile"] = bool(
        record["estimate"] > reference.percentile_95
        if kind == "auprc"
        else abs(record["estimate"]) > reference.percentile_95
    )

    if baseline_table is not None and baseline_arm is not None:
        baseline_scope = "all_lineage" if baseline_arm == "metadata" else scope_name
        baseline = scoped(
            load_predictions(
                baseline_table,
                arm=baseline_arm,
                lineage_scope=baseline_scope,
                endpoint=endpoint,
            ),
            roster,
        )
        paired = pooled_donor_bootstrap(
            kind,
            predictions,
            outcomes,
            roster,
            baseline=baseline,
            n_resamples=BOOTSTRAP_REPLICATES,
            seed=SEED,
        )
        record["paired_difference_vs_" + baseline_arm] = paired.to_dict()
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--predictions-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    root = arguments.root
    endpoints_root = root / "executions/gse296875-phenotype-endpoints-20260825"
    folds_root = root / "executions/gse296875-donor-folds-20260825"
    if digest(endpoints_root / "ARTIFACTS.json") != BOUND["phenotype_endpoints"]:
        raise ScoringError("endpoint lock changed")
    if digest(folds_root / "ARTIFACTS.json") != BOUND["donor_folds"]:
        raise ScoringError("donor fold lock changed")
    if digest(arguments.predictions / "ARTIFACTS.json") != arguments.predictions_sha256:
        raise ScoringError("prediction artifact changed since it was frozen")

    endpoint_rows = read_tsv(endpoints_root / "endpoint_lock" / "donor_endpoints.tsv")
    fold_rows = read_tsv(folds_root / "split_lock" / "donor_folds.tsv")
    outcomes = {
        endpoint: endpoint_values(endpoint, endpoint_rows)
        for endpoint in METRIC_BY_ENDPOINT
    }

    table = arguments.predictions / "predictions" / "predictions.tsv"
    adult_table = arguments.predictions / "predictions" / "predictions_adult_refit.tsv"
    well_table = (
        arguments.predictions / "predictions" / "predictions_leave_one_well_out.tsv"
    )

    results: list[dict[str, object]] = []

    # Confirmatory: molecular arm, pooled scope, twelve tests.
    family = confirmatory_family(PRIMARY_LINEAGES_DEFAULT)
    family_records: list[dict[str, object]] = []
    for arm, scope_name, endpoint in family:
        roster = registered_roster(
            endpoint, "pooled_all_donors", endpoint_rows, fold_rows
        )
        record = score_one(
            table,
            arm,
            scope_name,
            endpoint,
            roster,
            outcomes[endpoint],
            baseline_table=table,
            baseline_arm="metadata",
        )
        record["confirmatory"] = True
        family_records.append(record)
    corrected = benjamini_hochberg(
        {
            f"{r['lineage_scope']}|{r['endpoint']}": r["permutation_p_value"]
            for r in family_records
            if r.get("permutation_p_value") is not None
        }
    )
    for record in family_records:
        key = f"{record['lineage_scope']}|{record['endpoint']}"
        record["bh_q_value"] = corrected.get(key)
        record["bh_significant_at_q_0_05"] = (
            corrected.get(key) is not None and corrected[key] < 0.05
        )
    results.extend(family_records)

    # Comparator arms and secondary lineages: reported, never family members.
    for endpoint in METRIC_BY_ENDPOINT:
        roster = registered_roster(
            endpoint, "pooled_all_donors", endpoint_rows, fold_rows
        )
        for arm in ("metadata", "molecular_metadata"):
            scopes = ("all_lineage",) if arm == "metadata" else (
                "all_lineage",
                *PRIMARY_LINEAGES_DEFAULT,
            )
            for scope_name in scopes:
                record = score_one(
                    table, arm, scope_name, endpoint, roster, outcomes[endpoint]
                )
                record["confirmatory"] = False
                record["role"] = "comparator_arm"
                results.append(record)
        for scope_name in SECONDARY_LINEAGES:
            record = score_one(
                table, "molecular", scope_name, endpoint, roster, outcomes[endpoint]
            )
            record["confirmatory"] = False
            record["role"] = "secondary_lineage_not_a_gate_member"
            results.append(record)

    # Sensitivities, outside the family, nominal only.
    sensitivities: list[dict[str, object]] = []
    for endpoint in METRIC_BY_ENDPOINT:
        adult_roster = registered_roster(
            endpoint, "adult_only_scored", endpoint_rows, fold_rows
        )
        for scope_name in ("all_lineage", *PRIMARY_LINEAGES_DEFAULT):
            record = score_one(
                table, "molecular", scope_name, endpoint, adult_roster, outcomes[endpoint]
            )
            record["role"] = "adult_only_scored_from_all_age_model"
            sensitivities.append(record)
            refit_roster = registered_roster(
                endpoint, "adult_only_refit", endpoint_rows, fold_rows
            )
            record = score_one(
                adult_table,
                "molecular",
                scope_name,
                endpoint,
                refit_roster,
                outcomes[endpoint],
            )
            record["role"] = "adult_only_refit_inside_adult_folds"
            sensitivities.append(record)
        for well in sorted({row["well_id"] for row in fold_rows}):
            well_roster = registered_roster(
                endpoint, f"leave_one_well_out_{well}", endpoint_rows, fold_rows
            )
            record = score_one(
                well_table,
                "molecular",
                "all_lineage",
                endpoint,
                well_roster,
                outcomes[endpoint],
            )
            record["role"] = "leave_one_well_out_descriptive_only"
            record["well_is_not_an_independent_biological_stratum"] = True
            sensitivities.append(record)

    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "confirmatory_and_comparator_results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (arguments.output / "sensitivity_results.json").write_text(
        json.dumps(sensitivities, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    significant = [r for r in family_records if r["bh_significant_at_q_0_05"]]
    audit = {
        "schema_version": "masld-bench-gse296875-phenotype-scores-v1",
        "dataset_id": "gse296875",
        "unit_of_inference": "donor",
        "pooled_out_of_fold_metric_computed_once_per_endpoint": True,
        "per_fold_metric_averaging_is_unreachable_in_code": True,
        "confirmatory_family_size": len(family_records),
        "confirmatory_arm": "molecular",
        "multiple_testing": "BH_q_0.05_within_GSE296875_histopathology_secondary_family",
        "confirmatory_significant_at_q_0_05": len(significant),
        "significant_scopes": [
            f"{r['lineage_scope']}|{r['endpoint']}" for r in significant
        ],
        "comparator_and_secondary_records": len(results) - len(family_records),
        "sensitivity_records": len(sensitivities),
        "undefined_sensitivity_records": sum(
            1 for r in sensitivities if r.get("estimate") is None
        ),
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "seed": SEED,
        "prevalence_is_not_the_auprc_null": True,
        "champion_eligible": False,
        "role": "secondary_development_only",
        "status": "pass_pooled_scores",
    }
    (arguments.output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
