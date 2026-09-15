"""The aspect-separability prespecification must be checkable, not merely present.

A criterion derived after seeing the answer is not a criterion, and a
prespecification that only *could* have been written first proves nothing.
These tests pin the three criteria, their thresholds and the marginal label
distributions that were inspected before the freeze, so that a later edit to
the criteria is a test failure rather than a silent rewrite of the gate.

They also pin the exclusions that decide what the three aspects are:
``lobular_necrosis`` is deposited in the endpoint table but is not a NASH-CRN
NAS component, and fibrosis is a fourth axis that is never folded into the
three.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "freeze_aspect_separability_prespecification.py"
ENDPOINTS = (
    ROOT
    / "executions"
    / "model-data-061-21079623"
    / "activation"
    / "participant_endpoints.tsv"
)
SPEC = importlib.util.spec_from_file_location("aspect_prespec_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
prespec = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prespec)

THREE_ASPECTS = ("steatosis", "ballooning", "lobular_inflammation")


def _endpoint_rows() -> list[dict[str, str]]:
    lines = ENDPOINTS.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


class PrespecificationShapeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.payload = prespec.build()

    def test_the_three_criteria_are_present_and_named(self) -> None:
        self.assertEqual(
            sorted(self.payload["criteria"]),
            [
                "c1_labels_are_not_redundant",
                "c2_effective_dimensionality",
                "c3_aspects_rank_genes_differently",
            ],
        )

    def test_thresholds_are_pinned(self) -> None:
        criteria = self.payload["criteria"]
        self.assertEqual(
            criteria["c1_labels_are_not_redundant"]["threshold"],
            "every pair |rho| < 0.80",
        )
        self.assertEqual(
            criteria["c2_effective_dimensionality"]["threshold"],
            ">= 2.0 of a possible 3.0",
        )
        self.assertEqual(
            criteria["c3_aspects_rank_genes_differently"]["threshold"],
            "at least one pair with |rho| < 0.80",
        )

    def test_every_threshold_is_declared_a_judgment_call(self) -> None:
        for name, criterion in self.payload["criteria"].items():
            with self.subTest(criterion=name):
                self.assertIs(criterion["threshold_is_a_judgment_call"], True)

    def test_the_gene_ordering_criterion_is_the_decisive_one(self) -> None:
        criteria = self.payload["criteria"]
        self.assertIs(criteria["c3_aspects_rank_genes_differently"]["decisive"], True)
        for name in ("c1_labels_are_not_redundant", "c2_effective_dimensionality"):
            with self.subTest(criterion=name):
                self.assertNotIn("decisive", criteria[name])

    def test_a_stop_is_recorded_as_a_result(self) -> None:
        rule = self.payload["decision_rule"]
        self.assertEqual(rule["go"], "c1 AND c2 AND c3 all met")
        self.assertIs(rule["a_stop_is_a_result_not_a_failure"], True)
        self.assertIn("stop_labels_collapse", rule)
        self.assertIn("stop_same_gene_ordering", rule)

    def test_necrosis_and_fibrosis_exclusions_are_recorded(self) -> None:
        controls = " ".join(self.payload["controls_that_apply"])
        self.assertIn("lobular_necrosis", controls)
        self.assertIn("not a NASH-CRN NAS component", controls)
        self.assertIn("fibrosis is analysed as a fourth axis", controls)
        self.assertIn("measured permutation nulls", controls)

    def test_no_gene_expression_was_opened_before_the_freeze(self) -> None:
        inspected = self.payload["what_was_inspected_before_freezing"]
        self.assertIs(inspected["no_gene_expression_was_opened"], True)
        self.assertIs(inspected["no_cross_tabulation_or_correlation_was_computed"], True)
        self.assertIs(inspected["marginal_label_distributions_only"], True)

    def test_the_unit_of_inference_is_the_participant(self) -> None:
        self.assertEqual(self.payload["cohort"]["unit_of_inference"], "participant")
        self.assertEqual(self.payload["cohort"]["n_participants"], 99)


class RecordedMarginalsMatchTheEndpointTableTests(unittest.TestCase):
    """The recorded marginals must be the deposited ones, not remembered ones."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = _endpoint_rows()
        cls.payload = prespec.build()

    def test_a_row_count_is_not_a_unit_count(self) -> None:
        self.assertEqual(len(self.rows), 99)
        self.assertEqual(len({row["participant_id"] for row in self.rows}), 99)

    def test_recorded_marginals_reproduce(self) -> None:
        inspected = self.payload["what_was_inspected_before_freezing"]
        for column in THREE_ASPECTS + ("fibrosis",):
            observed: dict[str, int] = {}
            for row in self.rows:
                observed[row[column]] = observed.get(row[column], 0) + 1
            with self.subTest(column=column):
                self.assertEqual(observed, inspected[column])

    def test_the_all_zero_tied_block_reproduces(self) -> None:
        zeros = sum(1 for row in self.rows if row["nash_crn_component_sum"] == "0")
        self.assertEqual(zeros, 24)
        self.assertEqual(
            self.payload["what_was_inspected_before_freezing"][
                "nash_crn_component_sum_zero_participants"
            ],
            24,
        )

    def test_the_component_sum_is_the_three_aspects_and_excludes_necrosis(self) -> None:
        three_matches = 0
        with_necrosis_matches = 0
        for row in self.rows:
            total = sum(int(row[column]) for column in THREE_ASPECTS)
            deposited = int(row["nash_crn_component_sum"])
            three_matches += total == deposited
            with_necrosis_matches += (
                total + int(row["lobular_necrosis"]) == deposited
            )
        self.assertEqual(three_matches, 99)
        self.assertLess(with_necrosis_matches, 99)

    def test_the_endpoint_table_the_prespec_names_is_the_one_read(self) -> None:
        self.assertEqual(
            ROOT / self.payload["cohort"]["endpoint_table"], ENDPOINTS
        )


class PrespecificationWriteGuardTests(unittest.TestCase):
    def test_it_refuses_to_overwrite_an_existing_prespecification(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "nested" / "prespec.json"
            argv = sys.argv
            try:
                sys.argv = ["prespec", "--output", str(target)]
                self.assertEqual(prespec.main(), 0)
                loaded = json.loads(target.read_text(encoding="utf-8"))
                self.assertEqual(
                    loaded["prespec_id"], "gse267145_aspect_separability_v1"
                )
                with self.assertRaises(prespec.PrespecError):
                    prespec.main()
            finally:
                sys.argv = argv


if __name__ == "__main__":
    unittest.main()
