from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest
import unittest.mock

from scripts.build_gse83452_participant_masks import (
    FORBIDDEN_COLUMNS,
    IDENTITY_ALLOWLIST,
    STRUCTURALLY_WITHHELD,
    ParticipantMaskError,
    build_masks,
    project_identity_columns,
)

COLUMNS = [
    "cohort_family_id", "series", "participant_id", "sample_accession",
    "source_sample_number", "paired_baseline_source_number", "paired_baseline_accession",
    "timepoint", "repeat_topology", "nash_status", "age", "sex", "intervention",
    "scan_date", "primary_transfer_evaluable", "paired_expression_stress_evaluable",
    "paired_nash_transition_evaluable",
]


def write_table(path: Path, *, records: int = 6, participants: int = 4) -> list[str]:
    accessions = [f"GSM{2203254 + i}" for i in range(records)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, accession in enumerate(accessions):
            writer.writerow({
                "cohort_family_id": "antwerp_inserm_shared",
                "series": "GSE83452",
                "participant_id": f"P{index % participants:03d}",
                "sample_accession": accession,
                "source_sample_number": str(index),
                "paired_baseline_source_number": "",
                "paired_baseline_accession": "",
                "timepoint": "baseline" if index < participants else "followup",
                "repeat_topology": "paired" if index >= participants else "single",
                # Outcome and covariates are poisoned: any leak is visible.
                "nash_status": "NASH", "age": "55", "sex": "female",
                "intervention": "treated", "scan_date": "2015-01-01",
                "primary_transfer_evaluable": "true",
                "paired_expression_stress_evaluable": "true",
                "paired_nash_transition_evaluable": "true",
            })
    return accessions


class ProjectionTests(unittest.TestCase):
    def test_allowlist_and_forbidden_columns_are_disjoint(self) -> None:
        self.assertEqual(set(IDENTITY_ALLOWLIST) & set(FORBIDDEN_COLUMNS), set())

    def test_projection_drops_outcome_and_covariates(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            path = base / "records.tsv"
            write_table(path)
            with unittest.mock.patch(
                "scripts.build_gse83452_participant_masks.EXPECTED_RECORDS", 6
            ):
                rows = project_identity_columns(path)
            for row in rows:
                self.assertEqual(sorted(row), sorted(IDENTITY_ALLOWLIST))
                for column in FORBIDDEN_COLUMNS:
                    self.assertNotIn(column, row)


class MaskTests(unittest.TestCase):
    def _build(self, base: Path, *, records: int = 6, participants: int = 4):
        path = base / "records.tsv"
        accessions = write_table(path, records=records, participants=participants)
        with unittest.mock.patch(
            "scripts.build_gse83452_participant_masks.EXPECTED_RECORDS", records
        ), unittest.mock.patch(
            "scripts.build_gse83452_participant_masks.EXPECTED_PARTICIPANTS", participants
        ):
            return build_masks(records_path=path, summarized_accessions=accessions)

    def test_repeats_share_one_outer_group(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            masks = self._build(Path(value))
            self.assertEqual(masks["records"], 6)
            self.assertEqual(masks["participants"], 4)
            self.assertEqual(masks["distinct_outer_groups"], 4)
            self.assertEqual(masks["participants_with_repeat_records"], 2)
            self.assertIs(masks["one_array_per_participant"], False)
            groups = {}
            for record in masks["mask_records"]:
                groups.setdefault(record["participant_id"], set()).add(record["outer_group_key"])
            for participant, keys in groups.items():
                self.assertEqual(len(keys), 1, participant)

    def test_outcome_and_covariates_are_withheld_not_imputed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            masks = self._build(Path(value))
            self.assertIs(masks["withheld_fields_imputed"], False)
            self.assertIs(masks["outcome_columns_read"], False)
            self.assertIs(masks["labels_read"], False)
            for record in masks["mask_records"]:
                for field in STRUCTURALLY_WITHHELD:
                    self.assertIsNone(record[field])
                    self.assertIs(record[f"{field}_observed"], False)
                self.assertNotIn("NASH", str(record.values()))

    def test_cohort_family_firewall_is_recorded(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            masks = self._build(Path(value))
            self.assertEqual(masks["accession_aliases"], ["GSE106737", "GSE83452"])
            self.assertIs(masks["independent_evaluation_between_accession_aliases"], False)
            self.assertIn("78 fingerprint", masks["family_firewall_reason"])

    def test_participant_count_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            path = base / "records.tsv"
            accessions = write_table(path, records=6, participants=4)
            with unittest.mock.patch(
                "scripts.build_gse83452_participant_masks.EXPECTED_RECORDS", 6
            ):
                with self.assertRaisesRegex(ParticipantMaskError, "not 171"):
                    build_masks(records_path=path, summarized_accessions=accessions)

    def test_summarized_arrays_must_match_the_record_table(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            path = base / "records.tsv"
            accessions = write_table(path, records=6, participants=4)
            with unittest.mock.patch(
                "scripts.build_gse83452_participant_masks.EXPECTED_RECORDS", 6
            ), unittest.mock.patch(
                "scripts.build_gse83452_participant_masks.EXPECTED_PARTICIPANTS", 4
            ):
                with self.assertRaisesRegex(ParticipantMaskError, "do not match"):
                    build_masks(
                        records_path=path,
                        summarized_accessions=accessions[:-1] + ["GSM9999999"],
                    )


if __name__ == "__main__":
    unittest.main()
