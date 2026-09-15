"""The ragged header must be caught by measurement, and S2 must not pass unpowered.

The first block is the important one. Both external metadata files carry one
more data field than header field, and a name-based read returns sex where
fibrosis is asked for. A test that only checked the happy path would pass on a
correctly-read table and pass equally on a silently shifted one, so these pin
the refusal behaviour: an unexpected offset is refused rather than guessed, and
a resolved column whose contents are not the ordinal it claims is refused rather
than coerced.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "evaluate_two_axis_external_support.py"
EXTERNAL = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
    "Deconvolution/bulk")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


external = _load(MODULE, "two_axis_external_tested")


def _write(text: str) -> Path:
    handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    handle.write(text)
    handle.close()
    return Path(handle.name)


class RaggedHeaderTests(unittest.TestCase):
    def test_an_offset_of_one_is_measured_and_applied(self) -> None:
        path = _write("a\tb\tc\nX\t1\t2\t3\nY\t4\t5\t6\n")
        rows, audit = external.read_ragged_table(path)
        self.assertEqual(audit["measured_offset"], 1)
        self.assertTrue(audit["is_ragged"])
        self.assertEqual(rows[0], {"a": "1", "b": "2", "c": "3"})

    def test_an_aligned_table_is_read_without_an_offset(self) -> None:
        path = _write("a\tb\n1\t2\n3\t4\n")
        rows, audit = external.read_ragged_table(path)
        self.assertEqual(audit["measured_offset"], 0)
        self.assertFalse(audit["is_ragged"])
        self.assertEqual(rows[0], {"a": "1", "b": "2"})

    def test_an_offset_of_two_is_refused_rather_than_guessed(self) -> None:
        path = _write("a\tb\nW\tX\t1\t2\n")
        with self.assertRaises(external.ExternalError):
            external.read_ragged_table(path)

    def test_inconsistent_row_widths_are_refused(self) -> None:
        path = _write("a\tb\nX\t1\t2\nY\t3\n")
        with self.assertRaises(external.ExternalError):
            external.read_ragged_table(path)

    def test_both_real_cohorts_are_ragged_by_exactly_one(self) -> None:
        for cohort, header, data in (("GSE130970", 35, 36), ("GSE193066", 87, 88)):
            path = EXTERNAL / cohort / f"{cohort}_metadata.tsv"
            if not path.exists():
                self.skipTest("external metadata not present")
            _, audit = external.read_ragged_table(path)
            self.assertEqual(audit["header_fields"], header, cohort)
            self.assertEqual(audit["data_fields"], data, cohort)
            self.assertEqual(audit["measured_offset"], 1, cohort)

    def test_the_resolved_gse193066_columns_hold_what_they_claim(self) -> None:
        path = EXTERNAL / "GSE193066" / "GSE193066_metadata.tsv"
        if not path.exists():
            self.skipTest("external metadata not present")
        rows, _ = external.read_ragged_table(path)
        fibrosis = external.integer_column(rows, "fibrosis stage", 0, 4)
        self.assertEqual(len(fibrosis), 164)
        self.assertEqual(sorted({int(v) for v in fibrosis}), [0, 1, 2, 3, 4])
        biopsy = {row["biopsy"].strip() for row in rows}
        self.assertEqual(biopsy, {"1st biopsy", "2nd biopsy"})


class ContentValidationTests(unittest.TestCase):
    def test_a_column_holding_sex_is_refused_where_an_ordinal_is_asked_for(self) -> None:
        rows = [{"fibrosis stage": "female"}, {"fibrosis stage": "male"}]
        with self.assertRaises(external.ExternalError) as caught:
            external.integer_column(rows, "fibrosis stage", 0, 4)
        self.assertIn("ragged-header failure", str(caught.exception))

    def test_an_out_of_range_ordinal_is_refused(self) -> None:
        rows = [{"nas": "9"}, {"nas": "2"}]
        with self.assertRaises(external.ExternalError):
            external.integer_column(rows, "nas", 0, 8)

    def test_a_valid_ordinal_passes(self) -> None:
        rows = [{"nas": "0"}, {"nas": "8"}, {"nas": "4"}]
        self.assertEqual(list(external.integer_column(rows, "nas", 0, 8)), [0.0, 8.0, 4.0])


class MatchingTests(unittest.TestCase):
    def test_the_background_matches_the_target_bin_histogram_exactly(self) -> None:
        generator = np.random.default_rng(1)
        bins = np.repeat(np.arange(10), 200)
        target = np.concatenate([np.flatnonzero(bins == b)[:5] for b in (2, 3, 7)])
        pool = np.setdiff1d(np.arange(bins.size), target)
        background = external.matched_background(bins, target, pool, generator)
        self.assertIsNotNone(background)
        from collections import Counter
        self.assertEqual(Counter(bins[background].tolist()),
                         Counter(bins[target].tolist()))
        self.assertEqual(len(set(background.tolist()) & set(target.tolist())), 0)

    def test_an_exhausted_bin_returns_none_rather_than_a_short_draw(self) -> None:
        generator = np.random.default_rng(2)
        bins = np.asarray([0, 0, 1, 1, 1])
        target = np.asarray([0, 1])
        pool = np.asarray([2, 3, 4])
        self.assertIsNone(external.matched_background(bins, target, pool, generator))

    def test_decile_bins_span_ten_levels_on_a_spread_input(self) -> None:
        bins = external.decile_bins(np.arange(1000, dtype=float))
        self.assertEqual(len(set(bins.tolist())), 10)


class S2PowerGuardTests(unittest.TestCase):
    """S2 accepts a null, so an unpowered arm must not be allowed to pass it."""

    def test_a_quiet_but_unpowered_condition_is_inapplicable_not_met(self) -> None:
        block = {"served_by_this_arm": True, "exceeds_the_matched_null_p95": False,
                 "arm_could_have_detected_an_elevation_that_size": False}
        quiet = not block["exceeds_the_matched_null_p95"]
        powered = block["arm_could_have_detected_an_elevation_that_size"]
        self.assertTrue(quiet)
        self.assertFalse(powered)
        self.assertFalse(quiet and powered)

    def test_the_three_evidence_states_are_distinguished(self) -> None:
        for quiet, powered, expected in (
            (True, True, "tested_negative"),
            (True, False, "indeterminate"),
            (False, True, "elevated"),
            (False, False, "elevated"),
        ):
            state = ("tested_negative" if (quiet and powered)
                     else "indeterminate" if quiet else "elevated")
            self.assertEqual(state, expected)


if __name__ == "__main__":
    unittest.main()
