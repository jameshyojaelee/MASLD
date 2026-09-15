from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.build_gse49541_participant_masks import (
    FORBIDDEN_COLUMNS,
    IDENTITY_ALLOWLIST,
    STRUCTURALLY_MISSING,
    ParticipantMaskError,
    build_masks,
    project_identity_columns,
)

FULL_COLUMNS = [
    "cohort_family_id",
    "series",
    "participant_id",
    "sample_accession",
    "source_sample_key",
    "timepoint",
    "repeat_topology",
    "fibrosis_stage_group",
    "exact_fibrosis_stage",
    "nash_status",
    "age",
    "sex",
    "intervention",
    "primary_transfer_evaluable",
    "paired_expression_stress_evaluable",
    "paired_nash_transition_evaluable",
]


def write_table(path: Path, records: int = 72, series: str = "GSE49541") -> list[str]:
    accessions = [f"GSM{789110 + index}" for index in range(records)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FULL_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, accession in enumerate(accessions):
            writer.writerow(
                {
                    "cohort_family_id": "gse31803_gse49541_fibrosis_array",
                    "series": series,
                    "participant_id": f"P{index:03d}",
                    "sample_accession": accession,
                    "source_sample_key": f"key{index}",
                    "timepoint": "baseline",
                    "repeat_topology": "single_record",
                    # Outcome columns carry a poison value: any leak shows up.
                    "fibrosis_stage_group": "advanced_f3_f4" if index % 2 else "mild_f0_f1",
                    "exact_fibrosis_stage": "F3",
                    "nash_status": "NASH",
                    "age": "55",
                    "sex": "female",
                    "intervention": "none",
                    "primary_transfer_evaluable": "true",
                    "paired_expression_stress_evaluable": "false",
                    "paired_nash_transition_evaluable": "false",
                }
            )
    return accessions


class ProjectionTests(unittest.TestCase):
    def test_allowlist_and_forbidden_columns_are_disjoint(self) -> None:
        self.assertEqual(set(IDENTITY_ALLOWLIST) & set(FORBIDDEN_COLUMNS), set())

    def test_projection_drops_every_outcome_column(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "participants.tsv"
            write_table(path)
            rows = project_identity_columns(path)
            self.assertEqual(len(rows), 72)
            for row in rows:
                self.assertEqual(sorted(row), sorted(IDENTITY_ALLOWLIST))
                self.assertNotIn("fibrosis_stage_group", row)

    def test_projection_is_sorted_by_accession(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "participants.tsv"
            accessions = write_table(path)
            rows = project_identity_columns(path)
            self.assertEqual(
                [row["sample_accession"] for row in rows], sorted(accessions)
            )

    def test_wrong_record_count_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "participants.tsv"
            write_table(path, records=71)
            with self.assertRaisesRegex(ParticipantMaskError, "not 72 records"):
                project_identity_columns(path)

    def test_foreign_series_rows_are_not_counted(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "participants.tsv"
            write_table(path, series="GSE83452")
            with self.assertRaisesRegex(ParticipantMaskError, "not 72 records"):
                project_identity_columns(path)


class MaskTests(unittest.TestCase):
    def test_structurally_missing_fields_are_null_and_unobserved(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "participants.tsv"
            accessions = write_table(path)
            masks = build_masks(participants_path=path, summarized_accessions=accessions)
            self.assertEqual(masks["participants"], 72)
            self.assertIs(masks["structurally_missing_imputed"], False)
            self.assertIs(masks["outcome_columns_read"], False)
            self.assertIs(masks["labels_read"], False)
            for record in masks["records"]:  # type: ignore[union-attr]
                for field in STRUCTURALLY_MISSING:
                    self.assertIsNone(record[field])
                    self.assertIs(record[f"{field}_observed"], False)
                self.assertIs(record["prediction_eligible"], True)
                self.assertNotIn("fibrosis_stage_group", record)

    def test_summarized_arrays_must_match_the_participant_table(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "participants.tsv"
            accessions = write_table(path)
            with self.assertRaisesRegex(ParticipantMaskError, "do not match"):
                build_masks(
                    participants_path=path,
                    summarized_accessions=accessions[:-1] + ["GSM000000"],
                )


if __name__ == "__main__":
    unittest.main()
