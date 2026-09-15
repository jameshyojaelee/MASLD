#!/usr/bin/env python3
"""Freeze GSE296875 donor-grouped folds and the prespecified sensitivity splits.

This runs before any fitting and re-uses the donor-grouped outer split already
frozen by the fragment-membership campaign rather than drawing a new one; a
split redrawn after training output files exist would not be a prespecified split.

Three partitions are frozen together.  The donor-grouped outer folds put every
nucleus of a donor in one fold.  Leave-one-well-out is a transport sensitivity
over the technical batch, and because every donor sits in exactly one well,
holding a well out removes a whole donor block; the well effect is aliased with
that block and must never be read as an independent biological replicate.  The
adult-only partition is a mask over all thirty-nine donors, so the five donors
under eighteen stay visible in the output file instead of disappearing from it.

Endpoint availability is counted per fold and per well because a partition with
no fibrosis positive cannot yield a within-partition average precision.  That
is recorded, not repaired: both endpoints are cross-fitted, so the metric
belongs on pooled out-of-fold predictions.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import tomllib
from typing import Any

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


ADULT_AGE_YEARS = 18


class DonorFoldError(ValueError):
    """Raised when a frozen split contradicts the donor roster or the requirements."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _availability(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Count endpoint support without ever imputing an unobserved donor."""

    steatosis = [row for row in rows if row["steatosis_observed"]]
    fibrosis = [row for row in rows if row["fibrosis_observed"]]
    positives = sum(1 for row in fibrosis if row["fibrosis_any"])
    negatives = len(fibrosis) - positives
    return {
        "donors": len(rows),
        "steatosis_observed": len(steatosis),
        "steatosis_distinct_values": len({row["steatosis_numeric"] for row in steatosis}),
        "fibrosis_observed": len(fibrosis),
        "fibrosis_positive": positives,
        "fibrosis_negative": negatives,
        "adult_donors": sum(1 for row in rows if row["is_adult"] is True),
        "donors_under_18": sum(1 for row in rows if row["is_adult"] is False),
        # Average precision needs at least one positive and one negative, and a
        # rank correlation needs at least two distinct scores.
        "within_partition_auprc_defined": positives >= 1 and negatives >= 1,
        "within_partition_spearman_defined": (
            len(steatosis) >= 2 and len({row["steatosis_numeric"] for row in steatosis}) >= 2
        ),
    }


def build_partitions(
    endpoints: list[dict[str, Any]],
    folds: dict[str, int],
    nuclei_by_donor: dict[str, int],
    outer_folds: int,
) -> dict[str, Any]:
    roster = [row["donor_id"] for row in endpoints]
    if len(set(roster)) != len(roster):
        raise DonorFoldError("a donor appears more than once in the endpoint table")
    if sorted(folds) != sorted(roster):
        raise DonorFoldError("the frozen fold assignment does not cover the donor roster exactly")
    if sorted(nuclei_by_donor) != sorted(roster):
        raise DonorFoldError("the donor census does not cover the donor roster exactly")
    observed_folds = sorted(set(folds.values()))
    if observed_folds != list(range(outer_folds)):
        raise DonorFoldError(f"outer folds are not 0..{outer_folds - 1}: {observed_folds}")

    by_donor = {row["donor_id"]: row for row in endpoints}
    for row in endpoints:
        row["outer_fold"] = folds[row["donor_id"]]
        row["nuclei"] = nuclei_by_donor[row["donor_id"]]

    fold_rows = {
        fold: [row for row in endpoints if row["outer_fold"] == fold]
        for fold in observed_folds
    }
    wells = sorted({row["well_id"] for row in endpoints})
    well_rows = {well: [row for row in endpoints if row["well_id"] == well] for well in wells}

    fold_availability = []
    for fold in observed_folds:
        rows = fold_rows[fold]
        train = [row for row in endpoints if row["outer_fold"] != fold]
        fold_availability.append(
            {
                "outer_fold": fold,
                "nuclei": sum(row["nuclei"] for row in rows),
                "wells_represented": sorted({row["well_id"] for row in rows}),
                "test": _availability(rows),
                "train": _availability(train),
            }
        )
    well_availability = []
    for well in wells:
        rows = well_rows[well]
        train = [row for row in endpoints if row["well_id"] != well]
        well_availability.append(
            {
                "held_out_well": well,
                "nuclei": sum(row["nuclei"] for row in rows),
                "folds_represented": sorted({row["outer_fold"] for row in rows}),
                "test": _availability(rows),
                "train": _availability(train),
            }
        )

    adult = [row for row in endpoints if row["is_adult"] is True]
    unknown_age = [row for row in endpoints if row["is_adult"] is None]
    adult_fold_availability = [
        {
            "outer_fold": fold,
            "test": _availability([row for row in fold_rows[fold] if row["is_adult"] is True]),
        }
        for fold in observed_folds
    ]

    fold_auprc_undefined = [
        item["outer_fold"] for item in fold_availability
        if not item["test"]["within_partition_auprc_defined"]
    ]
    well_auprc_undefined = [
        item["held_out_well"] for item in well_availability
        if not item["test"]["within_partition_auprc_defined"]
    ]
    return {
        "endpoints": endpoints,
        "fold_availability": fold_availability,
        "well_availability": well_availability,
        "adult_fold_availability": adult_fold_availability,
        "summary": {
            "outer_folds": outer_folds,
            "donors": len(endpoints),
            "wells": len(wells),
            "nuclei": sum(row["nuclei"] for row in endpoints),
            "fold_sizes": {str(fold): len(fold_rows[fold]) for fold in observed_folds},
            "donors_per_well": {well: len(well_rows[well]) for well in wells},
            "donor_is_nested_within_well": True,
            "well_is_confounded_with_a_donor_block": True,
            "well_is_a_technical_batch_not_a_biological_replicate": True,
            "adult_only_donors": len(adult),
            "donors_under_18": sum(1 for row in endpoints if row["is_adult"] is False),
            "donors_under_18_ids": [
                row["donor_id"] for row in endpoints if row["is_adult"] is False
            ],
            "donors_with_unknown_age": len(unknown_age),
            "adult_only_fold_sizes": {
                str(item["outer_fold"]): item["test"]["donors"]
                for item in adult_fold_availability
            },
            "folds_without_a_within_fold_auprc": fold_auprc_undefined,
            "wells_without_a_within_well_auprc": well_auprc_undefined,
            "pooled_out_of_fold_metric_required": bool(fold_auprc_undefined),
            "per_fold_metric_averaging_is_invalid": bool(fold_auprc_undefined),
            "frozen_before_any_fitting": True,
            "under_18_donors_are_masked_not_dropped": True,
        },
    }


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    path.chmod(0o440)


def _flag(value: Any) -> str:
    if value is None:
        return ""
    return "true" if value else "false"


def run(
    *,
    endpoint_lock: Path,
    join_lock: Path,
    membership: Path,
    taskspec: Path,
    outer_folds: int,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise DonorFoldError("refusing to overwrite donor-fold artifact")
    with taskspec.open("rb") as handle:
        spec = tomllib.load(handle)
    endpoints = json.loads((endpoint_lock / "donor_endpoints.json").read_text(encoding="utf-8"))
    folds = {
        row["donor_id"]: int(row["outer_fold"])
        for row in read_tsv(membership / "donor_folds.tsv")
    }
    census = read_tsv(join_lock / "donor_census.tsv")
    nuclei_by_donor = {row["donor_id"]: int(row["nuclei"]) for row in census}
    contract = json.loads((membership / "contract.json").read_text(encoding="utf-8"))
    if contract["outer_folds"] != outer_folds or contract["split_id"] != "donor_outer":
        raise DonorFoldError("the frozen membership split identity differs")

    # The census carries one row per donor, so a donor that spanned wells would
    # appear twice.  Proving nesting here is what makes leave-one-well-out a
    # whole-donor-block hold-out rather than a within-donor split.
    census_wells: dict[str, list[str]] = {}
    for row in census:
        census_wells.setdefault(row["donor_id"], []).append(row["well_id"])
    spanning = sorted(donor for donor, wells in census_wells.items() if len(set(wells)) > 1)
    if spanning:
        raise DonorFoldError(f"a donor spans more than one well: {spanning}")
    disagreeing = sorted(
        row["donor_id"]
        for row in endpoints
        if census_wells.get(row["donor_id"], [None])[0] != row["well_id"]
    )
    if disagreeing:
        raise DonorFoldError(f"endpoint well disagrees with the frozen join: {disagreeing}")

    partitions = build_partitions(endpoints, folds, nuclei_by_donor, outer_folds)
    summary = partitions["summary"]
    if summary["donors"] != spec["analyzed_donors"]:
        raise DonorFoldError("donor count differs from the TaskSpec")
    if summary["wells"] != spec["technical_batches"]:
        raise DonorFoldError("well count differs from the TaskSpec")
    if summary["donors_under_18"] != spec["donors_under_18"]:
        raise DonorFoldError("donors under 18 differ from the TaskSpec")
    if summary["adult_only_donors"] + summary["donors_under_18"] != summary["donors"]:
        raise DonorFoldError("the adult mask does not partition the roster")

    output.mkdir(mode=0o750, parents=True)
    rows = partitions["endpoints"]
    _write_tsv(
        output / "donor_folds.tsv",
        [
            "donor_id",
            "well_id",
            "outer_fold",
            "nuclei",
            "steatosis_observed",
            "fibrosis_observed",
            "is_adult",
            "in_adult_only_sensitivity",
        ],
        [
            [
                row["donor_id"],
                row["well_id"],
                str(row["outer_fold"]),
                str(row["nuclei"]),
                _flag(row["steatosis_observed"]),
                _flag(row["fibrosis_observed"]),
                _flag(row["is_adult"]),
                _flag(row["is_adult"] is True),
            ]
            for row in rows
        ],
    )
    _write_tsv(
        output / "leave_one_well_out.tsv",
        [
            "held_out_well",
            "test_donors",
            "test_nuclei",
            "test_steatosis_observed",
            "test_fibrosis_observed",
            "test_fibrosis_positive",
            "train_donors",
            "train_fibrosis_positive",
            "within_well_auprc_defined",
        ],
        [
            [
                item["held_out_well"],
                str(item["test"]["donors"]),
                str(item["nuclei"]),
                str(item["test"]["steatosis_observed"]),
                str(item["test"]["fibrosis_observed"]),
                str(item["test"]["fibrosis_positive"]),
                str(item["train"]["donors"]),
                str(item["train"]["fibrosis_positive"]),
                _flag(item["test"]["within_partition_auprc_defined"]),
            ]
            for item in partitions["well_availability"]
        ],
    )
    _write_tsv(
        output / "fold_endpoint_availability.tsv",
        [
            "outer_fold",
            "test_donors",
            "test_nuclei",
            "test_steatosis_observed",
            "test_fibrosis_observed",
            "test_fibrosis_positive",
            "test_fibrosis_negative",
            "test_adult_donors",
            "within_fold_auprc_defined",
            "within_fold_spearman_defined",
        ],
        [
            [
                str(item["outer_fold"]),
                str(item["test"]["donors"]),
                str(item["nuclei"]),
                str(item["test"]["steatosis_observed"]),
                str(item["test"]["fibrosis_observed"]),
                str(item["test"]["fibrosis_positive"]),
                str(item["test"]["fibrosis_negative"]),
                str(item["test"]["adult_donors"]),
                _flag(item["test"]["within_partition_auprc_defined"]),
                _flag(item["test"]["within_partition_spearman_defined"]),
            ]
            for item in partitions["fold_availability"]
        ],
    )
    write_json_exclusive(
        output / "partitions.json",
        {
            "fold_availability": partitions["fold_availability"],
            "well_availability": partitions["well_availability"],
            "adult_fold_availability": partitions["adult_fold_availability"],
        },
    )
    receipt = {
        "dataset_id": "gse296875",
        "unit_of_inference": "donor",
        "split_id": contract["split_id"],
        "split_seed": contract["split_seed"],
        "split_origin": "reused_from_frozen_fragment_membership_campaign_not_redrawn",
        "membership_artifacts_sha256": sha256_file(membership / "ARTIFACTS.json"),
        "membership_donor_folds_sha256": sha256_file(membership / "donor_folds.tsv"),
        "endpoint_lock_artifacts_sha256": sha256_file(endpoint_lock / "ARTIFACTS.json"),
        "join_lock_artifacts_sha256": sha256_file(join_lock / "ARTIFACTS.json"),
        "well_policy": spec["covariates"]["well_policy"],
        "age_scope_policy": spec["covariates"]["age_scope_policy"],
        "race_policy": spec["covariates"]["race_policy"],
        "selection_use": spec["inference"]["selection_use"],
        "uncertainty": spec["inference"]["uncertainty"],
        "bootstrap_replicates": spec["inference"]["bootstrap_replicates"],
        "multiple_testing": spec["inference"]["multiple_testing"],
        **summary,
    }
    write_json_exclusive(output / "split_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse296875_donor_grouped_split_lock",
            "dataset_id": "gse296875",
            "unit_of_inference": "donor",
            "split_id": contract["split_id"],
            "split_seed": contract["split_seed"],
            "outer_folds": outer_folds,
            "donors": summary["donors"],
            "wells": summary["wells"],
            "frozen_before_any_fitting": True,
            "pooled_out_of_fold_metric_required": summary["pooled_out_of_fold_metric_required"],
            "under_18_donors_are_masked_not_dropped": True,
            "champion_eligible": False,
            "external_or_sealed": False,
        },
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint-lock", required=True, type=Path)
    parser.add_argument("--join-lock", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--taskspec", required=True, type=Path)
    parser.add_argument("--outer-folds", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run(
        endpoint_lock=arguments.endpoint_lock,
        join_lock=arguments.join_lock,
        membership=arguments.membership,
        taskspec=arguments.taskspec,
        outer_folds=arguments.outer_folds,
        output=arguments.output,
    )
    print(json.dumps({"output": arguments.output.as_posix(), "receipt": receipt}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
