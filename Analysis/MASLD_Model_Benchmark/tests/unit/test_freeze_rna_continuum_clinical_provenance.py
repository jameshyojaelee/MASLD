from __future__ import annotations

import unittest

from scripts.freeze_rna_continuum_clinical_provenance import (
    CROSS_COHORT_STAGE_COMPARABLE,
    FIBROSIS_SCALE_NATIVE,
    GSE213621_STAGE_RECODE,
    SEX_SOURCE,
    ProvenanceError,
    build_reproduction,
    derive_gse162694,
    derive_gse213621,
    parse_optional_int,
    summarize,
)


class ParseTests(unittest.TestCase):
    def test_recorded_na_tokens_are_missing_not_zero(self) -> None:
        for token in ("", "NA", "N/A", "None", "--"):
            self.assertIsNone(parse_optional_int(token))

    def test_stage_zero_is_a_value_not_a_missing_marker(self) -> None:
        self.assertEqual(parse_optional_int("0"), 0)

    def test_integral_float_is_accepted_and_fractional_is_not_rounded(self) -> None:
        self.assertEqual(parse_optional_int("3.0"), 3)
        self.assertIsNone(parse_optional_int("3.5"))


class GSE162694Tests(unittest.TestCase):
    def test_normal_liver_histology_is_missing_not_stage_zero(self) -> None:
        derived = derive_gse162694({"fibrosis_stage": "normal liver histology",
                                    "nas_score": "2"})
        self.assertIsNone(derived["fibrosis_stage"])
        self.assertEqual(derived["fibrosis_native_value"], "normal liver histology")

    def test_recorded_stage_zero_is_retained(self) -> None:
        derived = derive_gse162694({"fibrosis_stage": "0", "nas_score": "4"})
        self.assertEqual(derived["fibrosis_stage"], 0)

    def test_fibrosis_and_nas_missingness_are_independent(self) -> None:
        # The two complete sets are not nested: a sample may carry NAS without a stage.
        derived = derive_gse162694({"fibrosis_stage": "normal liver histology",
                                    "nas_score": "1"})
        self.assertIsNone(derived["fibrosis_stage"])
        self.assertEqual(derived["nas_score"], 1)


class GSE213621ScaleTests(unittest.TestCase):
    def test_collapsed_scale_has_no_level_four(self) -> None:
        self.assertNotIn(4, set(GSE213621_STAGE_RECODE.values()))

    def test_level_one_pools_f0_with_f1_and_level_three_pools_f3_with_f4(self) -> None:
        self.assertEqual(GSE213621_STAGE_RECODE["F0F1"], 1)
        self.assertEqual(GSE213621_STAGE_RECODE["F3F4"], 3)

    def test_nas_is_unavailable_cohort_wide(self) -> None:
        derived = derive_gse213621({"fibrotic_stage": "F2"})
        self.assertIsNone(derived["nas_score"])

    def test_unrecognised_native_value_becomes_missing_not_zero(self) -> None:
        self.assertIsNone(derive_gse213621({"fibrotic_stage": "F9"})["fibrosis_stage"])

    def test_collapsed_cohort_is_flagged_not_stage_comparable(self) -> None:
        self.assertFalse(CROSS_COHORT_STAGE_COMPARABLE["GSE213621"])
        self.assertEqual(FIBROSIS_SCALE_NATIVE["GSE213621"], "gse213621_collapsed_0_3")

    def test_kleiner_cohorts_are_comparable_to_each_other(self) -> None:
        self.assertTrue(CROSS_COHORT_STAGE_COMPARABLE["GSE130970"])
        self.assertTrue(CROSS_COHORT_STAGE_COMPARABLE["GSE135251"])


class SexProvenanceTests(unittest.TestCase):
    def test_both_continuum_evaluation_cohorts_carry_inferred_sex(self) -> None:
        self.assertEqual(SEX_SOURCE["GSE162694"], "annotated")
        self.assertEqual(SEX_SOURCE["GSE213621"], "inferred_kmeans")
        self.assertEqual(SEX_SOURCE["GSE135251"], "inferred_kmeans")


class ReproductionTests(unittest.TestCase):
    def _manifest(self) -> list[dict[str, str]]:
        return [
            {"sample_id": "SRR1", "analysis_unit_id": "GSE130970::SRR1",
             "dataset": "GSE130970", "fibrosis_stage": "2", "nas_score": "5"},
            {"sample_id": "SRR2", "analysis_unit_id": "GSE213621::SRR2",
             "dataset": "GSE213621", "fibrosis_stage": "3", "nas_score": "NA"},
        ]

    def _sources(self) -> dict[str, list[dict[str, str]]]:
        return {
            "GSE130970": [
                {"Run": "SRR1", "fibrosis_stage": "2", "nafld_activity_score": "5"},
                {"Run": "SRR9", "fibrosis_stage": "0", "nafld_activity_score": "1"},
            ],
            "GSE213621": [{"Run": "SRR2", "fibrotic_stage": "F3F4"}],
        }

    def test_matching_values_reproduce(self) -> None:
        rows, orphans = build_reproduction(self._manifest(), self._sources())
        self.assertTrue(all(row["fibrosis_agrees"] for row in rows))
        self.assertTrue(all(row["nas_agrees"] for row in rows))
        self.assertEqual([o["sra_run_accession"] for o in orphans], ["SRR9"])

    def test_every_row_names_the_sample_unit_and_discloses_the_absent_donor_key(self) -> None:
        rows, _ = build_reproduction(self._manifest(), self._sources())
        for row in rows:
            self.assertEqual(row["unit_of_observation"], "sequencing_sample")
            self.assertIsNone(row["n_donors"])
            self.assertEqual(row["n_donors_status"], "untestable_no_donor_key_on_disk")

    def test_a_changed_source_value_is_caught_not_silently_accepted(self) -> None:
        sources = self._sources()
        sources["GSE130970"][0]["fibrosis_stage"] = "4"
        rows, _ = build_reproduction(self._manifest(), sources)
        self.assertFalse(rows[0]["fibrosis_agrees"])
        self.assertEqual(summarize(rows)["totals"]["n_fibrosis_disagreements"], 1)

    def test_manifest_sample_missing_from_source_is_an_error(self) -> None:
        sources = self._sources()
        sources["GSE130970"] = []
        with self.assertRaises(ProvenanceError):
            build_reproduction(self._manifest(), sources)

    def test_duplicate_join_key_is_rejected(self) -> None:
        sources = self._sources()
        sources["GSE130970"].append({"Run": "SRR1", "fibrosis_stage": "2",
                                     "nafld_activity_score": "5"})
        with self.assertRaises(ProvenanceError):
            build_reproduction(self._manifest(), sources)

    def test_summary_counts_fibrosis_and_nas_completeness_separately(self) -> None:
        rows, _ = build_reproduction(self._manifest(), self._sources())
        totals = summarize(rows)["totals"]
        self.assertEqual(totals["n_samples"], 2)
        self.assertEqual(totals["n_fibrosis_complete"], 2)
        self.assertEqual(totals["n_nas_complete"], 1)
        self.assertEqual(totals["n_both_complete"], 1)


if __name__ == "__main__":
    unittest.main()
