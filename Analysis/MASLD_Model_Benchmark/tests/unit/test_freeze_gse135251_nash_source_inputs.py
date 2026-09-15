from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.freeze_gse135251_nash_source_inputs import (
    ARMS,
    EXPECTED_CENSUS,
    NASH_GROUPS,
    SourceFreezeError,
    assert_declared_census,
    assign_outer_folds,
    build_crosswalk,
    read_rowname_offset_tsv,
)

HEADER = ["sample_id", "disease", "Fibrosis_stage", "group_in_paper", "nas_score"]


def write_offset_table(path: Path, rows: list[list[str]], *, offset: bool = True) -> None:
    lines = ["\t".join(HEADER)]
    for index, row in enumerate(rows):
        payload = row if not offset else [f"SRR{index:07d}"] + row
        lines.append("\t".join(payload))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class RownameOffsetReaderTests(unittest.TestCase):
    def test_offset_rows_are_read_without_shifting_columns(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "metadata.tsv"
            write_offset_table(
                path, [["SRR0000000", "NAFLD", "2", "NASH_F2", "5"]]
            )
            header, rows = read_rowname_offset_tsv(path)
            self.assertEqual(header, HEADER)
            self.assertEqual(rows[0]["Fibrosis_stage"], "2")
            self.assertEqual(rows[0]["group_in_paper"], "NASH_F2")
            self.assertEqual(rows[0]["sample_id"], "SRR0000000")

    def test_a_table_without_the_offset_fails_closed(self) -> None:
        """A naive N-field table must raise, not silently parse one column over."""

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "metadata.tsv"
            write_offset_table(
                path, [["SRR0000000", "NAFLD", "2", "NASH_F2", "5"]], offset=False
            )
            with self.assertRaises(SourceFreezeError):
                read_rowname_offset_tsv(path)


class CensusTests(unittest.TestCase):
    def _rows(self) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        plan = [
            ("control", "Control", ["0"] * 6 + ["1", "2"]),
            ("NAFL", "NAFLD", ["0"] * 25 + ["1"] * 16),
            ("NASH_F0-F1", "NAFLD", ["0"] * 4 + ["1"] * 24),
            ("NASH_F2", "NAFLD", ["2"] * 47),
            ("NASH_F3", "NAFLD", ["3"] * 44),
            ("NASH_F4", "NAFLD", ["4"] * 12),
        ]
        nas_pool = [
            value
            for score, count in EXPECTED_CENSUS["nas_score"].items()
            for value in [score] * count
        ]
        for group, disease, stages in plan:
            for stage in stages:
                rows.append(
                    {
                        "group_in_paper": group,
                        "disease": disease,
                        "Fibrosis_stage": stage,
                        "nas_score": nas_pool[len(rows)],
                    }
                )
        return rows

    def test_declared_census_reproduces(self) -> None:
        audit = assert_declared_census(self._rows())
        self.assertEqual(audit["group_in_paper"]["NASH_F2"], 47)
        self.assertTrue(audit["group_in_paper_agrees_with_fibrosis_stage"])
        self.assertFalse(audit["age_field_present"])
        self.assertFalse(audit["sex_field_present"])

    def test_a_shifted_group_to_stage_pairing_is_rejected(self) -> None:
        rows = self._rows()
        for row in rows:
            if row["group_in_paper"] == "NASH_F3":
                row["Fibrosis_stage"] = "2"
                break
        with self.assertRaises(SourceFreezeError):
            assert_declared_census(rows)

    def test_a_wrong_census_is_rejected(self) -> None:
        rows = self._rows()
        rows[0]["group_in_paper"] = "NAFL"
        with self.assertRaises(SourceFreezeError):
            assert_declared_census(rows)


class CrosswalkTests(unittest.TestCase):
    ID_TO_NAME = {
        "ENSG00000000001": "AAA",
        "ENSG00000000002": "BBB",
        "ENSG00000000003": "BBB",
        "ENSG00000000004": "",
        "ENSG00000000005": "CCC",
    }
    NAME_TO_IDS = {
        "AAA": ["ENSG00000000001"],
        "BBB": ["ENSG00000000002", "ENSG00000000003"],
        "CCC": ["ENSG00000000005"],
    }

    def _resolve(self, keys: list[str]):
        return build_crosswalk(
            row_keys=keys, id_to_name=self.ID_TO_NAME, name_to_ids=self.NAME_TO_IDS
        )

    def test_one_to_one_symbol_and_direct_ensembl_resolve(self) -> None:
        resolved, states = self._resolve(["AAA", "ENSG00000000004"])
        self.assertEqual(resolved["AAA"], "ENSG00000000001")
        self.assertEqual(resolved["ENSG00000000004"], "ENSG00000000004")
        by_key = {record["row_key"]: record["mapping_state"] for record in states}
        self.assertEqual(by_key["AAA"], "symbol_one_to_one")
        self.assertEqual(by_key["ENSG00000000004"], "ensembl_direct")

    def test_a_symbol_naming_two_genes_is_ambiguous(self) -> None:
        resolved, states = self._resolve(["BBB"])
        self.assertNotIn("BBB", resolved)
        self.assertEqual(states[0]["mapping_state"], "ambiguous_symbol")

    def test_two_row_keys_landing_on_one_gene_both_go_unresolved(self) -> None:
        """A collision must drop both rows, never pick one over the other."""

        resolved, states = self._resolve(["CCC", "ENSG00000000005"])
        self.assertEqual(resolved, {})
        self.assertEqual(
            {record["mapping_state"] for record in states},
            {"join_unresolved_collision"},
        )

    def test_absent_and_versioned_identifiers_are_recorded_not_dropped(self) -> None:
        _, states = self._resolve(["ENSG00000009999", "ENSG00000000001.4", "ZZZ"])
        by_key = {record["row_key"]: record["mapping_state"] for record in states}
        self.assertEqual(by_key["ENSG00000009999"], "absent_from_gencode_v49")
        self.assertEqual(by_key["ENSG00000000001.4"], "versioned_identifier")
        self.assertEqual(by_key["ZZZ"], "unmapped_symbol")


class FoldTests(unittest.TestCase):
    def test_every_training_partition_keeps_both_classes(self) -> None:
        labels = np.asarray([1] * 131 + [0] * 41, dtype=np.int64)
        assignment = assign_outer_folds(labels, folds=5, seed=1701)
        self.assertEqual(len(assignment), len(labels))
        for fold in range(5):
            self.assertEqual(
                set(int(v) for v in labels[assignment != fold]), {0, 1}
            )

    def test_assignment_is_deterministic_under_its_seed(self) -> None:
        labels = np.asarray([1] * 20 + [0] * 10, dtype=np.int64)
        first = assign_outer_folds(labels, folds=5, seed=1701)
        second = assign_outer_folds(labels, folds=5, seed=1701)
        self.assertTrue(np.array_equal(first, second))


class ArmTests(unittest.TestCase):
    def test_the_gate_eligible_arm_excludes_the_healthy_control_group(self) -> None:
        arm = ARMS["nash_vs_nafl"]
        self.assertTrue(arm["gate_eligible_arm"])
        self.assertEqual(arm["negative_groups"], ["NAFL"])
        self.assertEqual(arm["excluded_groups"], ["control"])
        self.assertEqual(arm["expected_class_counts"], {"no_nash": 41, "nash": 131})

    def test_the_control_bearing_arm_is_recorded_but_not_gate_eligible(self) -> None:
        arm = ARMS["nash_vs_nafl_plus_control"]
        self.assertFalse(arm["gate_eligible_arm"])
        self.assertEqual(arm["expected_class_counts"], {"no_nash": 49, "nash": 131})

    def test_both_arms_use_the_same_positive_definition(self) -> None:
        self.assertEqual(
            ARMS["nash_vs_nafl"]["positive_groups"],
            ARMS["nash_vs_nafl_plus_control"]["positive_groups"],
        )
        self.assertEqual(tuple(ARMS["nash_vs_nafl"]["positive_groups"]), NASH_GROUPS)


if __name__ == "__main__":
    unittest.main()
