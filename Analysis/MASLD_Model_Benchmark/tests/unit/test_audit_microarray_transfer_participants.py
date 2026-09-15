from __future__ import annotations

import json
from pathlib import Path
import unittest

from scripts.audit_microarray_transfer_participants import (
    MicroarrayParticipantAuditError,
    build_gse49541,
    build_gse83452,
    parse_characteristics,
)


def source_row(
    series: str, accession: str, title: str, platform: str, characteristics: list[str]
) -> dict[str, str]:
    return {
        "series": series,
        "accession": accession,
        "title": title,
        "source_name": "liver",
        "organism": "Homo sapiens",
        "platform_id": platform,
        "library_strategy": "",
        "description_json": "[]",
        "characteristics_json": json.dumps(characteristics),
        "data_processing_json": "[]",
        "relations_json": "[]",
        "supplementary_files_json": "[]",
    }


class MicroarrayTransferParticipantAuditTests(unittest.TestCase):
    def test_characteristics_reject_duplicate_key(self) -> None:
        with self.assertRaises(MicroarrayParticipantAuditError):
            parse_characteristics(json.dumps(["age: 20", "age: 21"]))

    def test_gse49541_rejects_exact_stage_invention(self) -> None:
        rows = [
            source_row(
                "GSE49541",
                f"GSM{i}",
                f"NAFLD liver biopsy tissue {i}",
                "GPL570",
                [
                    "Stage: mild (fibrosis stage 0-1)" if i < 40 else "Stage: advanced (fibrosis stage 3-4)",
                    "tissue: liver",
                ],
            )
            for i in range(72)
        ]
        values = build_gse49541(rows)
        self.assertEqual({row["exact_fibrosis_stage"] for row in values}, {"structurally_missing"})

    def test_gse83452_rejects_missing_pair_anchor(self) -> None:
        row = source_row(
            "GSE83452",
            "GSM1",
            "liver biopsy 2 (1)",
            "GPL16686",
            [
                "sample name: patient 2 (1)",
                "liver status: no NASH",
                "type of intervention: Diet",
                "time: follow-up",
                "age: 40",
                "gender: female",
                "scan date: 01/01/2014",
                "tissue: liver biopsy",
            ],
        )
        with self.assertRaisesRegex(MicroarrayParticipantAuditError, "pointer"):
            build_gse83452([row])

    def test_frozen_records_separate_expression_stress_from_nash_transition(self) -> None:
        source = (
            Path(__file__).resolve().parents[2]
            / "executions/model-data-062-21079818/samples/GSE83452.tsv"
        )
        if not source.exists():
            self.skipTest("frozen source audit is unavailable")
        from scripts.audit_microarray_transfer_participants import read_samples

        values = build_gse83452(read_samples(source, "GSE83452"))
        self.assertEqual(
            sum(bool(row["paired_expression_stress_evaluable"]) for row in values),
            60,
        )
        self.assertEqual(
            sum(bool(row["paired_nash_transition_evaluable"]) for row in values),
            54,
        )
        self.assertEqual(
            {row["nash_status"] for row in values},
            {"nash", "no_nash", "undefined"},
        )
        self.assertTrue(
            all(
                row["participant_id"].startswith("gse83452::")
                for row in values
            )
        )
        self.assertEqual(
            {row["cohort_family_id"] for row in values},
            {"antwerp_inserm_shared"},
        )


if __name__ == "__main__":
    unittest.main()
