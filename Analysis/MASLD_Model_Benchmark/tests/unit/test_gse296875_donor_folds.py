#!/usr/bin/env python3
"""Checks for the GSE296875 donor-grouped folds and sensitivity partitions."""

from __future__ import annotations

from pathlib import Path
import unittest

from scripts import freeze_gse296875_donor_folds as splits


def _donor(
    donor: str,
    well: str,
    *,
    steatosis: float | None = 5.0,
    fibrosis: bool | None = False,
    adult: bool | None = True,
) -> dict:
    return {
        "donor_id": donor,
        "well_id": well,
        "steatosis_observed": steatosis is not None,
        "steatosis_numeric": steatosis,
        "fibrosis_observed": fibrosis is not None,
        "fibrosis_any": fibrosis,
        "is_adult": adult,
    }


class GSE296875DonorFoldTests(unittest.TestCase):
    def test_a_fold_with_no_positive_is_flagged_not_repaired(self) -> None:
        endpoints = [
            _donor("1", "well1", fibrosis=False),
            _donor("2", "well1", fibrosis=False),
            _donor("3", "well2", fibrosis=True),
            _donor("4", "well2", fibrosis=False),
        ]
        folds = {"1": 0, "2": 0, "3": 1, "4": 1}
        nuclei = {donor: 100 for donor in folds}
        result = splits.build_partitions(endpoints, folds, nuclei, 2)
        summary = result["summary"]
        self.assertEqual(summary["folds_without_a_within_fold_auprc"], [0])
        self.assertTrue(summary["pooled_out_of_fold_metric_required"])
        self.assertTrue(summary["per_fold_metric_averaging_is_invalid"])
        self.assertFalse(result["fold_availability"][0]["test"]["within_partition_auprc_defined"])
        self.assertTrue(result["fold_availability"][1]["test"]["within_partition_auprc_defined"])

    def test_under_18_donors_are_masked_and_still_present(self) -> None:
        endpoints = [
            _donor("1", "well1", adult=True),
            _donor("2", "well1", adult=False),
            _donor("3", "well2", adult=True),
            _donor("4", "well2", adult=True),
        ]
        folds = {"1": 0, "2": 0, "3": 1, "4": 1}
        result = splits.build_partitions(endpoints, folds, {d: 10 for d in folds}, 2)
        summary = result["summary"]
        self.assertEqual(summary["donors"], 4)
        self.assertEqual(summary["adult_only_donors"], 3)
        self.assertEqual(summary["donors_under_18"], 1)
        self.assertEqual(summary["donors_under_18_ids"], ["2"])
        self.assertTrue(summary["under_18_donors_are_masked_not_dropped"])
        self.assertEqual(len(result["endpoints"]), 4)
        self.assertEqual(summary["adult_only_fold_sizes"], {"0": 1, "1": 2})

    def test_missing_endpoints_never_count_toward_a_partition(self) -> None:
        endpoints = [
            _donor("1", "well1", steatosis=None, fibrosis=None),
            _donor("2", "well1", steatosis=5.0, fibrosis=True),
        ]
        folds = {"1": 0, "2": 0}
        result = splits.build_partitions(endpoints, folds, {"1": 10, "2": 10}, 1)
        test = result["fold_availability"][0]["test"]
        self.assertEqual(test["donors"], 2)
        self.assertEqual(test["steatosis_observed"], 1)
        self.assertEqual(test["fibrosis_observed"], 1)
        self.assertEqual(test["fibrosis_positive"], 1)
        self.assertEqual(test["fibrosis_negative"], 0)
        self.assertFalse(test["within_partition_auprc_defined"])

    def test_leave_one_well_out_removes_a_whole_donor_block(self) -> None:
        endpoints = [
            _donor("1", "well1"),
            _donor("2", "well1"),
            _donor("3", "well2", fibrosis=True),
        ]
        folds = {"1": 0, "2": 1, "3": 0}
        result = splits.build_partitions(endpoints, folds, {d: 10 for d in folds}, 2)
        self.assertTrue(result["summary"]["donor_is_nested_within_well"])
        self.assertTrue(result["summary"]["well_is_confounded_with_a_donor_block"])
        by_well = {item["held_out_well"]: item for item in result["well_availability"]}
        self.assertEqual(by_well["well1"]["test"]["donors"], 2)
        self.assertEqual(by_well["well1"]["train"]["donors"], 1)
        self.assertEqual(by_well["well2"]["test"]["fibrosis_positive"], 1)
        self.assertEqual(by_well["well2"]["train"]["fibrosis_positive"], 0)

    def test_a_repeated_donor_row_fails_closed(self) -> None:
        endpoints = [_donor("1", "well1"), _donor("1", "well2")]
        with self.assertRaisesRegex(splits.DonorFoldError, "more than once"):
            splits.build_partitions(endpoints, {"1": 0}, {"1": 10}, 1)

    def test_a_fold_assignment_that_misses_a_donor_fails_closed(self) -> None:
        endpoints = [_donor("1", "well1"), _donor("2", "well1")]
        with self.assertRaisesRegex(splits.DonorFoldError, "cover the donor roster"):
            splits.build_partitions(endpoints, {"1": 0}, {"1": 10, "2": 10}, 1)

    def test_a_single_valued_steatosis_fold_has_no_within_fold_spearman(self) -> None:
        endpoints = [_donor("1", "well1", steatosis=0.0), _donor("2", "well1", steatosis=0.0)]
        result = splits.build_partitions(endpoints, {"1": 0, "2": 0}, {"1": 10, "2": 10}, 1)
        self.assertFalse(
            result["fold_availability"][0]["test"]["within_partition_spearman_defined"]
        )


if __name__ == "__main__":
    unittest.main()
