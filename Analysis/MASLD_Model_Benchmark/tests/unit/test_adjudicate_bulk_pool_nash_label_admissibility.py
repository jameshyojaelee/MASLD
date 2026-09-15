from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.adjudicate_bulk_pool_nash_label_admissibility import (
    EXPECTED_CENSUS,
    PROHIBITED_NEGATIVE,
    PROHIBITED_POSITIVE,
    LabelAdmissibilityError,
    adjudicate_cohort,
    classify_label_values,
    read_rowname_offset_tsv,
)


def verdict_for(cohort: str) -> dict:
    arms = classify_label_values(cohort=cohort, values=EXPECTED_CENSUS[cohort])
    return adjudicate_cohort(cohort=cohort, arms=arms)


class Gse162694Tests(unittest.TestCase):
    """NASH positives exist, but the only negative arm is a healthy control."""

    def test_it_does_not_pass(self) -> None:
        result = verdict_for("GSE162694")
        self.assertEqual(result["verdict"], "does_not_pass")
        self.assertFalse(result["admissible_as_a_nash_training_source"])

    def test_it_has_an_explicit_nash_positive(self) -> None:
        result = verdict_for("GSE162694")
        self.assertTrue(result["has_explicit_nash_positive"])
        self.assertEqual(sum(result["nash_values"].values()), 112)

    def test_it_fails_only_on_the_negative_clause(self) -> None:
        result = verdict_for("GSE162694")
        clauses = [item["clause"] for item in result["failed_clauses"]]
        self.assertEqual(clauses, ["[source_training] negative_required"])
        self.assertEqual(
            result["failed_clauses"][0]["prohibited_substitution"], "healthy_control"
        )

    def test_its_only_negative_is_the_control_arm(self) -> None:
        result = verdict_for("GSE162694")
        self.assertTrue(result["only_negative_is_a_healthy_control_arm"])
        self.assertEqual(result["healthy_control_values"], {"Control": 31})
        self.assertEqual(result["explicit_no_nash_values"], {})


class Gse240729Tests(unittest.TestCase):
    """Fibrosis stage only, which the contract names as a prohibited positive."""

    def test_it_does_not_pass(self) -> None:
        result = verdict_for("GSE240729")
        self.assertEqual(result["verdict"], "does_not_pass")

    def test_it_fails_the_positive_clause_on_fibrosis_stage_only(self) -> None:
        result = verdict_for("GSE240729")
        substitutions = [
            item["prohibited_substitution"] for item in result["failed_clauses"]
        ]
        self.assertIn("fibrosis_stage_only", substitutions)
        self.assertTrue(result["label_is_fibrosis_stage_only"])
        self.assertFalse(result["has_explicit_nash_positive"])

    def test_it_has_no_negative_arm_at_all(self) -> None:
        result = verdict_for("GSE240729")
        self.assertEqual(result["healthy_control_values"], {})
        self.assertEqual(result["explicit_no_nash_values"], {})
        self.assertEqual(result["clauses_failed"], 2)


class Gse213621Tests(unittest.TestCase):
    """Fibrosis stage only AND a healthy control arm: both clauses fail."""

    def test_it_does_not_pass(self) -> None:
        result = verdict_for("GSE213621")
        self.assertEqual(result["verdict"], "does_not_pass")

    def test_it_fails_both_clauses(self) -> None:
        result = verdict_for("GSE213621")
        clauses = sorted(item["clause"] for item in result["failed_clauses"])
        self.assertEqual(
            clauses,
            ["[source_training] negative_required", "[source_training] positive_required"],
        )

    def test_the_pooled_stage_values_are_not_read_as_nash(self) -> None:
        result = verdict_for("GSE213621")
        self.assertEqual(
            sorted(result["fibrosis_stage_or_other_values"]), ["F0F1", "F2", "F3F4"]
        )
        self.assertEqual(result["healthy_control_values"], {"Control": 68})


class PassingShapeTests(unittest.TestCase):
    def test_a_cohort_with_nash_and_nafl_arms_passes(self) -> None:
        """GSE135251's admitted shape, as the positive control for this logic."""

        arms = classify_label_values(
            cohort="GSE135251",
            values={"NASH_F2": 47, "NASH_F3": 44, "NAFL": 41, "control": 8},
        )
        result = adjudicate_cohort(cohort="GSE135251", arms=arms)
        self.assertEqual(result["verdict"], "passes")
        self.assertTrue(result["admissible_as_a_nash_training_source"])
        self.assertEqual(result["explicit_no_nash_values"], {"NAFL": 41})

    def test_a_control_only_negative_never_passes(self) -> None:
        arms = classify_label_values(
            cohort="X", values={"NASH_F2": 47, "control": 8}
        )
        self.assertEqual(adjudicate_cohort(cohort="X", arms=arms)["verdict"], "does_not_pass")


class ContractLiteralTests(unittest.TestCase):
    def test_every_cited_substitution_is_a_contract_literal(self) -> None:
        allowed = set(PROHIBITED_POSITIVE) | set(PROHIBITED_NEGATIVE)
        for cohort in ("GSE162694", "GSE240729", "GSE213621"):
            for failure in verdict_for(cohort)["failed_clauses"]:
                self.assertIn(failure["prohibited_substitution"], allowed)

    def test_an_invented_substitution_is_rejected(self) -> None:
        arms = classify_label_values(cohort="X", values={"NASH_F2": 5, "NAFL": 5})
        arms["has_admissible_negative"] = False
        arms["only_negative_is_a_healthy_control_arm"] = False
        result = adjudicate_cohort(cohort="X", arms=arms)
        self.assertEqual(
            result["failed_clauses"][0]["prohibited_substitution"], "missing"
        )


class OffsetReaderTests(unittest.TestCase):
    def test_a_table_without_the_offset_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "m.tsv"
            path.write_text("a\tb\nx\t1\t2\n", encoding="utf-8")
            header, rows = read_rowname_offset_tsv(path)
            self.assertEqual(rows[0], {"a": "1", "b": "2"})
            path.write_text("a\tb\n1\t2\n", encoding="utf-8")
            with self.assertRaises(LabelAdmissibilityError):
                read_rowname_offset_tsv(path)


if __name__ == "__main__":
    unittest.main()
