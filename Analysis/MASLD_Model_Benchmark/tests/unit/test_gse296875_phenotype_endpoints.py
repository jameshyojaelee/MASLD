#!/usr/bin/env python3
"""Checks for the GSE296875 endpoint rows and typed missingness masks."""

from __future__ import annotations

from pathlib import Path
import tomllib
import unittest

from scripts import freeze_gse296875_phenotype_endpoints as endpoints


ROOT = Path(__file__).parents[2]
TASKSPEC = ROOT / "config/evaluation/gse296875_histopathology.toml"


def _record(donor: str, **overrides) -> dict[str, str | None]:
    record = {
        "sample_id": donor,
        "well": "well1",
        "reported_sex": "F",
        "age_in_yr": "42",
        "reported_race": "W",
        "height_m": "1.65",
        "weight_kg": "70.0",
        "BMI": "25.7",
        "cause_of_death": "Anoxia",
        "genetic_similarity": "EUR",
        "fibrosis_categorical": "None",
        "fibrosis_status": "No significant fibrosis",
        "steatosis_categorical": "None",
        "steatosis_numeric": "0",
        "steatosis_status": "0",
    }
    record.update(overrides)
    return record


class GSE296875PhenotypeEndpointTests(unittest.TestCase):
    def setUp(self) -> None:
        with TASKSPEC.open("rb") as handle:
            self.spec = tomllib.load(handle)

    def test_the_na_sentinel_becomes_a_typed_mask_never_a_zero(self) -> None:
        records = [
            _record("331"),
            _record(
                "151",
                steatosis_numeric="NA",
                steatosis_categorical="NA",
                steatosis_status="NA",
                fibrosis_categorical="NA",
                fibrosis_status="NA",
            ),
        ]
        rows, summary = endpoints.build_endpoints(records, ["331", "151"], self.spec)
        missing = rows[1]
        self.assertFalse(missing["steatosis_observed"])
        self.assertIsNone(missing["steatosis_numeric"])
        self.assertIsNone(missing["steatosis_categorical"])
        self.assertFalse(missing["fibrosis_observed"])
        self.assertIsNone(missing["fibrosis_any"])
        self.assertNotEqual(missing["steatosis_numeric"], 0)
        self.assertNotEqual(missing["fibrosis_any"], False)
        self.assertEqual(summary["steatosis_observed_donors"], 1)
        self.assertEqual(summary["steatosis_missing_donors"], 1)
        self.assertEqual(summary["fibrosis_missing_donor_ids"], ["151"])

    def test_a_true_zero_steatosis_stays_observed_and_distinct_from_missing(self) -> None:
        records = [_record("331", steatosis_numeric="0"), _record("342", steatosis_numeric="NA")]
        rows, _ = endpoints.build_endpoints(records, ["331", "342"], self.spec)
        self.assertTrue(rows[0]["steatosis_observed"])
        self.assertEqual(rows[0]["steatosis_numeric"], 0.0)
        self.assertFalse(rows[1]["steatosis_observed"])
        self.assertIsNone(rows[1]["steatosis_numeric"])

    def test_fibrosis_uses_only_the_frozen_level_lists(self) -> None:
        records = [
            _record("331", fibrosis_categorical="None"),
            _record("342", fibrosis_categorical="mild"),
            _record("346", fibrosis_categorical="yes"),
        ]
        rows, summary = endpoints.build_endpoints(records, ["331", "342", "346"], self.spec)
        self.assertEqual([row["fibrosis_any"] for row in rows], [False, True, True])
        self.assertEqual(summary["fibrosis_positive_donors"], 2)
        self.assertEqual(summary["fibrosis_negative_donors"], 1)

    def test_an_unlisted_fibrosis_level_fails_closed_rather_than_becoming_negative(self) -> None:
        records = [_record("331", fibrosis_categorical="moderate")]
        with self.assertRaisesRegex(endpoints.PhenotypeEndpointError, "neither frozen list"):
            endpoints.build_endpoints(records, ["331"], self.spec)

    def test_no_ordinal_stage_is_inferred_from_a_prose_description(self) -> None:
        records = [_record("331", fibrosis_status="Mild portal and periportal fibrosis")]
        rows, _ = endpoints.build_endpoints(records, ["331"], self.spec)
        self.assertNotIn("fibrosis_stage", rows[0])
        self.assertIs(self.spec["fibrosis"]["ordinal_stage_inference"], False)
        self.assertEqual(self.spec["steatosis"]["recoding"], "none")

    def test_a_donor_deposited_but_not_analyzed_is_reported_and_excluded(self) -> None:
        records = [_record("331"), _record("381", well="well3")]
        rows, summary = endpoints.build_endpoints(records, ["331"], self.spec)
        self.assertEqual([row["donor_id"] for row in rows], ["331"])
        self.assertEqual(summary["deposited_donors_not_analyzed"], ["381"])

    def test_an_analyzed_donor_absent_from_the_supplement_fails_closed(self) -> None:
        with self.assertRaisesRegex(endpoints.PhenotypeEndpointError, "absent from the supplement"):
            endpoints.build_endpoints([_record("331")], ["331", "999"], self.spec)

    def test_adulthood_is_a_mask_not_a_silent_exclusion(self) -> None:
        records = [
            _record("331", age_in_yr="42"),
            _record("457", age_in_yr="13"),
            _record("342", age_in_yr="NA"),
        ]
        rows, summary = endpoints.build_endpoints(records, ["331", "457", "342"], self.spec)
        self.assertEqual([row["is_adult"] for row in rows], [True, False, None])
        self.assertEqual(summary["donors_under_18"], 1)
        self.assertEqual(summary["donors_under_18_ids"], ["457"])
        self.assertEqual(summary["adult_donors"], 1)
        self.assertEqual(summary["covariate_observed_donors"]["age_in_yr"], 2)

    def test_covariate_masks_are_counted_per_field(self) -> None:
        records = [_record("331"), _record("342", BMI="NA", reported_race="NA")]
        _, summary = endpoints.build_endpoints(records, ["331", "342"], self.spec)
        self.assertEqual(summary["covariate_observed_donors"]["BMI"], 1)
        self.assertEqual(summary["covariate_observed_donors"]["reported_race"], 1)
        self.assertEqual(summary["covariate_observed_donors"]["reported_sex"], 2)
        self.assertEqual(summary["reported_race_counts"], {"W": 1})


if __name__ == "__main__":
    unittest.main()
