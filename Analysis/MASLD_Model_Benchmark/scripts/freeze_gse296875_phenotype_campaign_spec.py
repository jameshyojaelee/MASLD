#!/usr/bin/env python3
"""Freeze the GSE296875 donor-by-lineage phenotype campaign specification.

Everything the campaign is allowed to do is decided here, before any molecular
data is read and before any model is fitted.  The power position, the
Benjamini-Hochberg family, the minimum-cell threshold, and the claim boundary
are written into a read-only output file so that a later result cannot quietly
redefine what was being asked.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from masld_bench.gse296875_phenotype_scoring import (  # noqa: E402
    confirmatory_family,
    endpoint_values,
    random_score_reference,
    read_tsv,
    registered_roster,
)


class CampaignSpecError(RuntimeError):
    """Raised when a bound input or a frozen decision does not verify."""


# Hash-bound inputs.  Every one of these is frozen upstream of this campaign.
BOUND_FIXTURES = {
    "donor_supplement": (
        "executions/gse296875-donor-supplement-20260825-v2",
        "06986191ab76f9d8bfba0dbda109ecfdeff05354f37a34fb9a5d4bf49dd50fae",
    ),
    "barcode_donor_join": (
        "executions/gse296875-barcode-donor-join-20260825",
        "ead9a54099fe2dfb47517497c072feec53fa5ebf9c0392e38ee34a164f8d64e8",
    ),
    "phenotype_endpoints": (
        "executions/gse296875-phenotype-endpoints-20260825",
        "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    ),
    "donor_folds": (
        "executions/gse296875-donor-folds-20260825",
        "e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416",
    ),
    "fragment_membership": (
        "executions/gse296875-fragment-membership-21063829",
        "582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64",
    ),
    "corgi_context_counts_cross_check": (
        "executions/corgi-gse296875-context-counts-21066278",
        "6946c329e74cd1a6cd72fe37faaf805183ff24c729ef9e549da11c49ca02fef8",
    ),
}

# The frozen partition.  These five are bound by hash into both the
# barcode-donor join receipt and the donor-folds receipt.
PRIMARY_LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")

# Permitted by `secondary_lineages_are_not_primary_gate_members = true`.
# Reported, never a check member, never in the confirmatory family.
SECONDARY_LINEAGES = ("endothelial_cell", "b_cell")

# Author label to frozen lineage.  Recorded so the two are reconcilable.
COLLAPSE_MAP = {
    "Cholangiocytes": "cholangiocyte",
    "Mesenchymal": "fibroblast",
    "Hepatocytes": "hepatocyte",
    "Kupffer": "macrophage",
    "NK-T": "t_cell",
    "LSEC": "endothelial_cell",
    "B cells": "b_cell",
}

MINIMUM_CELLS_PER_UNIT = 20
CANDIDATE_THRESHOLDS = (10, 20, 25, 30, 50)
BOOTSTRAP_REPLICATES = 10_000
SEED = 20260825


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_bound_fixtures(root: Path) -> dict[str, dict[str, str]]:
    verified: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in BOUND_FIXTURES.items():
        manifest = root / relative / "ARTIFACTS.json"
        observed = digest(manifest)
        if observed != expected:
            raise CampaignSpecError(
                f"bound fixture changed: {relative} expected {expected} got {observed}"
            )
        verified[name] = {"path": relative, "artifacts_sha256": observed}
    return verified


def lineage_census(membership_root: Path) -> tuple[dict[str, dict[str, int]], list[str]]:
    path = membership_root / "cell_membership.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    counts: dict[str, dict[str, int]] = {}
    donors: set[str] = set()
    observed_labels: set[str] = set()
    for row in rows:
        donors.add(row["donor_id"])
        observed_labels.add(row["source_label"])
        counts.setdefault(row["lineage_id"], {})
        counts[row["lineage_id"]][row["donor_id"]] = (
            counts[row["lineage_id"]].get(row["donor_id"], 0) + 1
        )
    if observed_labels != set(COLLAPSE_MAP):
        raise CampaignSpecError(
            f"author label set differs: {sorted(observed_labels)}"
        )
    if set(counts) != set(PRIMARY_LINEAGES) | set(SECONDARY_LINEAGES):
        raise CampaignSpecError(f"frozen lineage set differs: {sorted(counts)}")
    return counts, sorted(donors)


def census_rows(
    counts: dict[str, dict[str, int]], donors: list[str]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for lineage in PRIMARY_LINEAGES + SECONDARY_LINEAGES:
        per_donor = [counts[lineage].get(donor, 0) for donor in donors]
        rows.append(
            {
                "lineage_id": lineage,
                "analysis_role": (
                    "primary" if lineage in PRIMARY_LINEAGES else "secondary"
                ),
                "in_confirmatory_family": lineage in PRIMARY_LINEAGES,
                "total_nuclei": sum(per_donor),
                "donor_units_with_any_nuclei": sum(1 for v in per_donor if v > 0),
                "minimum_nuclei": min(per_donor),
                "median_nuclei": int(statistics.median(per_donor)),
                "maximum_nuclei": max(per_donor),
                "units_below_threshold": sum(
                    1 for v in per_donor if v < MINIMUM_CELLS_PER_UNIT
                ),
            }
        )
    return rows


def threshold_sensitivity(
    counts: dict[str, dict[str, int]], donors: list[str]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for threshold in CANDIDATE_THRESHOLDS:
        primary_masked = sum(
            1
            for donor in donors
            for lineage in PRIMARY_LINEAGES
            if counts[lineage].get(donor, 0) < threshold
        )
        secondary_masked = sum(
            1
            for donor in donors
            for lineage in SECONDARY_LINEAGES
            if counts[lineage].get(donor, 0) < threshold
        )
        by_lineage = {
            lineage: sum(
                1 for donor in donors if counts[lineage].get(donor, 0) < threshold
            )
            for lineage in PRIMARY_LINEAGES + SECONDARY_LINEAGES
        }
        rows.append(
            {
                "minimum_nuclei_per_unit": threshold,
                "selected": threshold == MINIMUM_CELLS_PER_UNIT,
                "primary_units_total": len(donors) * len(PRIMARY_LINEAGES),
                "primary_units_masked": primary_masked,
                "primary_units_observed": len(donors) * len(PRIMARY_LINEAGES)
                - primary_masked,
                "secondary_units_total": len(donors) * len(SECONDARY_LINEAGES),
                "secondary_units_masked": secondary_masked,
                "masked_by_lineage": json.dumps(by_lineage, sort_keys=True),
            }
        )
    return rows


def power_statement(root: Path) -> dict[str, object]:
    endpoint_rows = read_tsv(
        root
        / "executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv"
    )
    fold_rows = read_tsv(
        root / "executions/gse296875-donor-folds-20260825/split_lock/donor_folds.tsv"
    )
    steatosis = endpoint_values("steatosis", endpoint_rows)
    fibrosis = endpoint_values("fibrosis", endpoint_rows)
    steatosis_roster = registered_roster(
        "steatosis", "pooled_all_donors", endpoint_rows, fold_rows
    )
    fibrosis_roster = registered_roster(
        "fibrosis", "pooled_all_donors", endpoint_rows, fold_rows
    )
    adult_steatosis = registered_roster(
        "steatosis", "adult_only_refit", endpoint_rows, fold_rows
    )
    adult_fibrosis = registered_roster(
        "fibrosis", "adult_only_refit", endpoint_rows, fold_rows
    )

    steatosis_reference = random_score_reference(
        "spearman", steatosis, steatosis_roster, n_draws=BOOTSTRAP_REPLICATES, seed=SEED
    )
    fibrosis_reference = random_score_reference(
        "auprc", fibrosis, fibrosis_roster, n_draws=BOOTSTRAP_REPLICATES, seed=SEED
    )

    values = sorted(steatosis.values())
    positives = int(sum(fibrosis.values()))

    return {
        "position": (
            "Thirty-nine donors, fifteen fibrosis positives and twenty-two "
            "negatives, one cohort, project-exposed, development only."
        ),
        "unit_of_inference": "donor",
        "donors": 39,
        "inferential_n_is_never_the_unit_count": True,
        "steatosis": {
            "observed_donors": len(steatosis),
            "missing_donors": 39 - len(steatosis),
            "median_percent": statistics.median(values),
            "donors_at_zero_percent": sum(1 for v in values if v == 0.0),
            "donors_at_or_below_ten_percent": sum(1 for v in values if v <= 10.0),
            "maximum_percent": max(values),
            "donors_at_maximum": sum(1 for v in values if v == max(values)),
            "adult_only_observed_donors": len(adult_steatosis.donor_ids),
            "random_score_reference": steatosis_reference.to_dict(),
            "detectable_threshold_note": (
                "A continuous random scorer reaches |rho| of "
                f"{steatosis_reference.percentile_95:.3f} five percent of the "
                "time on this roster. Report that threshold beside every "
                "reported correlation."
            ),
            "leverage_note": (
                f"{sum(1 for v in values if v <= 10.0)} of {len(values)} donors "
                f"sit at or below ten percent and "
                f"{sum(1 for v in values if v == 0.0)} sit at exactly zero, so "
                f"the rank statistic is carried by the "
                f"{sum(1 for v in values if v >= 45.0)} donors above forty-five "
                "percent. This leverage must be reported in the results, not "
                "only in the specification."
            ),
        },
        "fibrosis": {
            "observed_donors": len(fibrosis),
            "missing_donors": 39 - len(fibrosis),
            "positives": positives,
            "negatives": len(fibrosis) - positives,
            "prevalence": positives / len(fibrosis),
            "adult_only_observed_donors": len(adult_fibrosis.donor_ids),
            "adult_only_positives": int(
                sum(fibrosis[donor] for donor in adult_fibrosis.donor_ids)
            ),
            "random_score_reference": fibrosis_reference.to_dict(),
            "prevalence_is_not_the_null": True,
            "baseline_note": (
                "Average precision is upward biased at this sample size. A "
                "continuous random scorer averages "
                f"{fibrosis_reference.mean:.3f} on this roster and exceeds "
                f"{fibrosis_reference.percentile_95:.3f} one time in twenty, so "
                f"the prevalence line of {positives / len(fibrosis):.3f} is not "
                "the null. Every reported AUPRC carries this reference beside "
                "it."
            ),
        },
        "can_detect": [
            "a large, lineage-localised association: Spearman |rho| at or above "
            "roughly 0.35 to 0.40 on steatosis, AUPRC at or above roughly 0.65 "
            "to 0.70 on fibrosis",
            "a clean, bounded negative on either endpoint",
        ],
        "cannot_detect": [
            "any modest association",
            "a between-arm difference smaller than roughly 0.15 AUPRC or 0.20 rho",
            "any interaction or covariate-adjusted effect",
            "any per-fold, per-well, or single-donor statement",
            "separation between lineages whose predictions are correlated",
            "whether a non-hepatocyte signal is ambient hepatocyte RNA rather "
            "than lineage-intrinsic biology; 46,286 of 68,398 nuclei are "
            "hepatocyte and these are raw, not decontaminated, counts",
        ],
        "acceptable_outcome": (
            "A well-characterised negative is an acceptable and publishable "
            "outcome of this campaign. The cohort is one hundred percent "
            "development-exposed and cannot confirm a champion, so a bounded "
            "null with a stated detection floor is worth more than an "
            "overstated positive."
        ),
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build_spec(root: Path) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    fixtures = verify_bound_fixtures(root)
    membership_root = root / BOUND_FIXTURES["fragment_membership"][0]
    counts, donors = lineage_census(membership_root)
    census = census_rows(counts, donors)
    sensitivity = threshold_sensitivity(counts, donors)

    selected = next(row for row in sensitivity if row["selected"])
    primary_nuclei = sum(
        row["total_nuclei"] for row in census if row["analysis_role"] == "primary"
    )
    secondary_nuclei = sum(
        row["total_nuclei"] for row in census if row["analysis_role"] == "secondary"
    )
    family = confirmatory_family(PRIMARY_LINEAGES)

    spec = {
        "schema_version": "masld-bench-gse296875-donor-lineage-phenotype-v1",
        "dataset_id": "gse296875",
        "role": "secondary_development_only",
        "champion_eligible": False,
        "external_or_sealed": False,
        "unit_of_inference": "donor",
        "donors": len(donors),
        "frozen_before_any_fitting": True,
        "bound_fixtures": fixtures,
        "lineages": {
            "partition": "frozen_fragment_membership_primary_five",
            "primary": list(PRIMARY_LINEAGES),
            "secondary": list(SECONDARY_LINEAGES),
            "secondary_are_reported_outside_the_confirmatory_family": True,
            "author_label_collapse_map": COLLAPSE_MAP,
            "primary_nuclei": primary_nuclei,
            "secondary_nuclei": secondary_nuclei,
            "total_nuclei": primary_nuclei + secondary_nuclei,
            "coverage_note": (
                f"The frozen primary five cover {primary_nuclei} nuclei. The "
                f"remaining {secondary_nuclei} are exactly LSEC "
                f"{sum(counts['endothelial_cell'].values())} plus B cell "
                f"{sum(counts['b_cell'].values())}. Those two are "
                "reported as labelled secondary scopes and are never gate "
                "members."
            ),
        },
        "fixture": {
            "kind": "new_full_transcriptome_donor_by_lineage_pseudobulk",
            "rationale": (
                "The Corgi context counts were selected for a different task. "
                "Inheriting that feature selection would bake another task's "
                "choice into this analysis, so the fixture is built fresh over "
                "all Gene Expression features and the Corgi artifact is bound "
                "as a declared cross-check instead."
            ),
            "counts": "raw Cell Ranger ARC UMIs, Gene Expression features only",
            "decontaminated": False,
            "minimum_nuclei_per_unit": MINIMUM_CELLS_PER_UNIT,
            "primary_units_total": selected["primary_units_total"],
            "primary_units_observed": selected["primary_units_observed"],
            "primary_units_masked": selected["primary_units_masked"],
            "masked_state": "insufficient_cells",
            "masked_units_carry_observed_cell_count": True,
            "missing_is_never_zero_and_never_an_absent_row": True,
            "threshold_burden_note": (
                "The threshold is not spread evenly. Cholangiocyte and t_cell "
                "carry almost the whole loss, so the two thinnest lineages "
                "absorb the cost of the choice. State this in the results "
                "rather than in a footnote."
            ),
        },
        "arms": {
            "molecular": {
                "features": "fold-fitted principal components of lineage pseudobulk",
                "lineage_observed_indicator": True,
            },
            "metadata": {"features": ["age_in_yr", "reported_sex", "BMI"]},
            "molecular_plus_metadata": {
                "features": "molecular principal components concatenated with metadata"
            },
            "race_is_not_fitted": True,
            "race_policy": (
                "Reported descriptively as an error-audit stratum only. With "
                "five Black donors a stable coefficient is not estimable, and "
                "the frozen TaskSpec forbids fitting one."
            ),
        },
        "preprocessing_fit_inside_training_folds_only": [
            "per-unit counts per million",
            "log1p",
            "detection filter at fifty percent of training units of that lineage",
            "top two thousand training-variance highly variable genes",
            "training mean and standard deviation z-score",
            "principal components, k = min(20, n_train - 1)",
            "missing lineage block imputed at the training-fold lineage mean "
            "with a binary observed indicator, never zero",
        ],
        "estimators": {
            "steatosis": {
                "model": "ridge_regression",
                "alpha_grid": [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0],
                "inner_selection": "leave_one_donor_out_within_training_fold",
                "inner_criterion": "mean_squared_error",
            },
            "fibrosis": {
                "model": "l2_logistic_regression",
                "C_grid": [0.001, 0.01, 0.1, 1.0, 10.0],
                "inner_selection": "leave_one_donor_out_within_training_fold",
                "inner_criterion": "log_loss",
                "calibration": "none; average precision is rank based",
            },
        },
        "metrics": {
            "steatosis": "cross_fitted_spearman",
            "fibrosis": "cross_fitted_auprc",
            "computed_once_on_pooled_out_of_fold_predictions": True,
            "per_fold_averaging_is_a_different_estimand": True,
            "enforcement_module": "src/masld_bench/gse296875_phenotype_scoring.py",
            "enforcement": [
                "the frozen prediction artifact carries no metric",
                "load_predictions projects the table to donor_id and prediction, "
                "so the fold column is never bound and a per-fold split is "
                "inexpressible downstream",
                "registered_roster serves only a closed scope registry that "
                "contains no per-fold scope and cannot be extended at call time",
                "pooled_metric requires exact roster equality, so one fold's "
                "donors raise as a proper subset and a masked donor raises as a "
                "superset",
                "no per-fold metric function is defined, so there is nothing to "
                "call by accident",
                "the fitting process runs on a Python without the benchmark "
                "package installed and cannot import the metric module at all",
            ],
            "fold_3_note": (
                "Outer fold 3 holds four donors and zero fibrosis positives, so "
                "a within-fold AUPRC is undefined there. The pooled requirement "
                "is independent of that and would hold even if every fold were "
                "defined."
            ),
        },
        "uncertainty": {
            "method": "pooled_donor_bootstrap",
            "contract_field": "paired_donor_cluster_bootstrap",
            "contract_satisfied": True,
            "existing_implementation_inapplicable": (
                "evaluators.stats.paired_cluster_bootstrap reduces within each "
                "cluster and then averages the differences, which assumes an "
                "additive metric. Pooled Spearman and pooled average precision "
                "are not additive over donors, so a new implementation with its "
                "own tests is supplied. The contract is met; only that "
                "algorithm does not apply."
            ),
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": SEED,
            "resampling_unit": "donor",
            "degenerate_resamples": "counted and reported, never silently dropped",
        },
        "multiple_testing": {
            "procedure": "BH_q_0.05_within_GSE296875_histopathology_secondary_family",
            "family_size": len(family),
            "family": [
                {"arm": arm, "lineage_scope": scope, "endpoint": endpoint}
                for arm, scope, endpoint in family
            ],
            "comparator_arms_are_not_family_members": True,
            "comparator_note": (
                "Metadata-only and molecular-plus-metadata enter as paired "
                "difference intervals against the molecular arm. They are "
                "comparators, not independent hypotheses. Carrying all three "
                "arms would make thirty-six tests at n=39, which guarantees "
                "nothing survives and reports a power limit as if it were a "
                "result."
            ),
        },
        "sensitivities_outside_the_family": {
            "adult_only_scored": (
                "all-age predictions restricted to the adult mask; asks how the "
                "fitted model behaves on adults"
            ),
            "adult_only_refit": (
                "refit inside the frozen adult-only folds; asks what an "
                "adult-only campaign would have found. A different question "
                "from the above and labelled as such."
            ),
            "leave_one_well_out": (
                "descriptive only. Donor is perfectly nested within well, with "
                "zero donors spanning wells, so a held-out well is the same "
                "donors under a different partition and confounds batch with "
                "donor composition. AUPRC is undefined for well3 and well5, "
                "which have no positives, and for well7, which has no "
                "negatives."
            ),
        },
        "ambient_rna_diagnostic": {
            "registered_before_fitting": True,
            "blocking": False,
            "risk": (
                "The ambient pool in liver single-nucleus RNA is dominated by "
                "hepatocyte transcript. A non-hepatocyte lineage could track "
                "steatosis through contamination alone, which would break "
                "precisely the claim this campaign exists to make."
            ),
            "diagnostic": (
                "For every non-hepatocyte unit, report the UMI fraction in the "
                "top fifty hepatocyte-restricted genes, defined inside training "
                "folds only, and report whether the association survives "
                "adjustment for that fraction."
            ),
            "escalation": (
                "If a non-hepatocyte association appears, the decontaminated "
                "dXbg re-analysis becomes the load-bearing check."
            ),
        },
        "power": power_statement(root),
        "claims": {
            "allowed": (
                "cross-sectional histology-associated representation or "
                "prediction in this donor cohort"
            ),
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
        "hotspot_programs": {
            "used_as_targets": False,
            "used_as_selected_features": False,
            "used_as_tuning_endpoints": False,
            "projection_permitted_after_selection_is_frozen": True,
        },
        "label_firewall": {
            "labels_are_evaluator_only_on_every_model_path": True,
            "predictions_frozen_and_hash_committed_before_outcome_join": True,
        },
        "missingness": "explicit_endpoint_mask_never_absence",
        "missing_sentinel": "NA",
        "deposited_donors_not_analyzed": ["381"],
        "donors_under_18_are_masked_not_dropped": True,
        "environment_note": (
            "PYTHONNOUSERSITE=1 is mandatory. The user site-packages sklearn is "
            "built against NumPy 1.x and hard-fails against the project "
            "environment's NumPy 2.2.6 with AttributeError: _ARRAY_API not "
            "found. Without the flag this costs an hour to diagnose."
        ),
        "python_split_note": (
            "Fitting runs on the scanpy environment, Python 3.10, which has "
            "NumPy and scikit-learn but cannot import the benchmark package "
            "because tomllib is 3.11 and later. Scoring runs on module Python "
            "3.11 with the benchmark package and no NumPy. The metric module is "
            "therefore unimportable from the fitting process, which reinforces "
            "the label firewall by construction rather than by discipline."
        ),
    }
    return spec, census, sensitivity


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    spec, census, sensitivity = build_spec(arguments.root)
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "campaign_spec.json").write_text(
        json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (arguments.output / "power_statement.json").write_text(
        json.dumps(spec["power"], indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_tsv(arguments.output / "donor_lineage_census.tsv", census)
    write_tsv(arguments.output / "min_cell_sensitivity.tsv", sensitivity)
    print(
        json.dumps(
            {
                "status": "campaign_spec_frozen",
                "donors": spec["donors"],
                "primary_lineages": list(PRIMARY_LINEAGES),
                "secondary_lineages": list(SECONDARY_LINEAGES),
                "minimum_nuclei_per_unit": MINIMUM_CELLS_PER_UNIT,
                "primary_units_observed": spec["fixture"]["primary_units_observed"],
                "primary_units_masked": spec["fixture"]["primary_units_masked"],
                "family_size": spec["multiple_testing"]["family_size"],
                "outcomes_read": False,
                "model_training_activated": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
