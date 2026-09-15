from __future__ import annotations

import unittest

from scripts.adjudicate_gse135251_source_roster_admission import (
    GROUP_ORDER,
    NASH_GROUPS,
    AdmissionAdjudicationError,
    cross_tabulate,
    nas_threshold_separability,
    read_rowname_offset_tsv,
)

# The deposited NAS distribution, as re-derived from the frozen copy.
DEPOSITED = {
    "control": {0: 8},
    "NAFL": {1: 6, 2: 16, 3: 9, 4: 8, 5: 2},
    "NASH_F0-F1": {3: 5, 4: 5, 5: 10, 6: 8},
    "NASH_F2": {3: 3, 4: 11, 5: 16, 6: 10, 7: 5, 8: 2},
    "NASH_F3": {3: 3, 4: 9, 5: 11, 6: 12, 7: 7, 8: 2},
    "NASH_F4": {1: 1, 3: 2, 4: 1, 5: 3, 6: 2, 7: 1, 8: 2},
}


def expand(group: str) -> list[int]:
    return [value for value, count in DEPOSITED[group].items() for _ in range(count)]


class SeparabilityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.nafl = expand("NAFL")
        self.nash = [value for group in NASH_GROUPS for value in expand(group)]

    def test_the_deposited_census_is_131_against_41(self) -> None:
        self.assertEqual(len(self.nash), 131)
        self.assertEqual(len(self.nafl), 41)

    def test_no_nas_threshold_reproduces_the_deposited_call(self) -> None:
        """This is the whole adjudication: a derived rule would score exactly 1.0."""

        result = nas_threshold_separability(nafl=self.nafl, nash=self.nash)
        self.assertFalse(result["call_is_a_nas_sum_threshold"])
        self.assertLess(result["best_accuracy"], 1.0)
        self.assertGreater(result["best_misclassified"], 0)

    def test_the_classes_overlap_on_the_nas_sum(self) -> None:
        result = nas_threshold_separability(nafl=self.nafl, nash=self.nash)
        self.assertEqual(result["overlapping_nas_values"], [1, 3, 4, 5])
        self.assertGreater(result["nafl_inside_overlap"], 0)
        self.assertGreater(result["nash_inside_overlap"], 0)

    def test_a_genuinely_threshold_derived_call_is_detected(self) -> None:
        """If the deposit had been a threshold, the guard must catch it."""

        nafl = [0, 1, 2, 2, 3]
        nash = [4, 5, 6, 7, 8]
        result = nas_threshold_separability(nafl=nafl, nash=nash)
        self.assertTrue(result["call_is_a_nas_sum_threshold"])
        self.assertEqual(result["best_accuracy"], 1.0)
        self.assertEqual(result["overlapping_nas_values"], [])

    def test_the_best_rule_is_reported_with_its_error_counts(self) -> None:
        result = nas_threshold_separability(nafl=self.nafl, nash=self.nash)
        row = next(
            item
            for item in result["grid"]
            if item["nas_threshold"] == result["best_threshold"]
        )
        self.assertEqual(
            row["misclassified"], row["false_positive"] + row["false_negative"]
        )
        self.assertEqual(result["evaluated_participants"], 172)

    def test_an_empty_contrast_is_rejected(self) -> None:
        with self.assertRaises(AdmissionAdjudicationError):
            nas_threshold_separability(nafl=[], nash=[])


class CrossTabTests(unittest.TestCase):
    def test_cross_tabulation_covers_every_group_in_order(self) -> None:
        rows = [
            {"group_in_paper": group, "Stage": "early"}
            for group in GROUP_ORDER
        ]
        emitted, table = cross_tabulate(rows, field="Stage")
        self.assertEqual([row["group_in_paper"] for row in emitted], list(GROUP_ORDER))
        self.assertEqual(set(table), set(GROUP_ORDER))

    def test_a_field_determined_by_stage_is_visible_as_such(self) -> None:
        """source_name tracks fibrosis stage; the record must not call it a population field."""

        rows = (
            [{"group_in_paper": "NAFL", "source_name": "NAFLD_early_liver biopsy"}] * 41
            + [{"group_in_paper": "NASH_F3", "source_name": "NAFLD_moderate_liver biopsy"}] * 44
        )
        _, table = cross_tabulate(rows, field="source_name")
        self.assertEqual(len(table["NAFL"]), 1)
        self.assertEqual(len(table["NASH_F3"]), 1)
        self.assertNotEqual(set(table["NAFL"]), set(table["NASH_F3"]))


class OffsetReaderTests(unittest.TestCase):
    def test_the_rowname_offset_is_required(self) -> None:
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "metadata.tsv"
            path.write_text("a\tb\nx\t1\t2\n", encoding="utf-8")
            header, rows = read_rowname_offset_tsv(path)
            self.assertEqual(header, ["a", "b"])
            self.assertEqual(rows[0], {"a": "1", "b": "2"})
            path.write_text("a\tb\n1\t2\n", encoding="utf-8")
            with self.assertRaises(AdmissionAdjudicationError):
                read_rowname_offset_tsv(path)


if __name__ == "__main__":
    unittest.main()
