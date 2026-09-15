#!/usr/bin/env python3
"""Synthetic fail-closed tests for external-cohort readiness topology."""

from __future__ import annotations

import copy
import unittest

from scripts.freeze_external_cohort_development_readiness import (
    ExternalCohortReadinessError,
    GSE105127_FOLD_FIELDS,
    GSE268273_PLAN_FIELDS,
    validate_gse105127_fold_rows,
    validate_gse105127_plan_rows,
    validate_gse268273_topology,
)


class ExternalCohortReadinessTests(unittest.TestCase):
    @staticmethod
    def gse105127_rows():
        rna = []
        rrbs = []
        for participant in range(19):
            bundle = str(participant % 8)
            for zone_index, zone in enumerate(("CV", "IZ", "PP")):
                index = participant * 3 + zone_index
                common = {
                    "bundle_id": bundle,
                    "row_id": f"row_{index:02d}",
                    "participant_group_id": f"person_{participant:02d}",
                    "zone": zone,
                    "pairing_topology": "adjacent_section",
                }
                rna.append(
                    {
                        **common,
                        "run_accession": f"SRR{index:07d}",
                        "library_layout": "SINGLE",
                        "read_length": "76",
                        "read_count": "1000",
                        "fastq_bytes": "2000",
                        "fastq_md5": f"{index:032x}",
                    }
                )
                rrbs.append(dict(common))
        return rna, rrbs

    def test_gse105127_accepts_adjacent_section_participant_topology(self) -> None:
        rna, rrbs = self.gse105127_rows()
        observed = validate_gse105127_plan_rows(
            tuple(rna[0]), rna, tuple(rrbs[0]), rrbs
        )
        self.assertEqual(observed["participants"], 19)
        self.assertEqual(observed["pairing_topology"], "adjacent_section")

    def test_gse105127_rejects_false_same_section_pairing(self) -> None:
        rna, rrbs = self.gse105127_rows()
        rna[0]["pairing_topology"] = "same_section"
        rrbs[0]["pairing_topology"] = "same_section"
        with self.assertRaises(ExternalCohortReadinessError):
            validate_gse105127_plan_rows(
                tuple(rna[0]), rna, tuple(rrbs[0]), rrbs
            )

    def test_gse105127_rejects_participant_fold_split(self) -> None:
        rows = []
        for participant in range(19):
            fold = participant % 5
            for zone in ("CV", "IZ", "PP"):
                rows.append(
                    {
                        "row_id": f"{participant}_{zone}",
                        "participant_group_id": f"person_{participant:02d}",
                        "zone": zone,
                        "participant_outer_fold": str(fold),
                        "fold_assignment_uses_outcomes": "False",
                        "available_to_model_as_feature": "False",
                    }
                )
        self.assertEqual(
            sorted(validate_gse105127_fold_rows(GSE105127_FOLD_FIELDS, rows).values()),
            [3, 4, 4, 4, 4],
        )
        rows[1]["participant_outer_fold"] = "4"
        with self.assertRaises(ExternalCohortReadinessError):
            validate_gse105127_fold_rows(GSE105127_FOLD_FIELDS, rows)

    @staticmethod
    def gse268273_rows():
        rows = []
        run_counts = [4] * 36 + [8] * 49 + [12] * 24
        for index, runs in enumerate(run_counts):
            layout = "paired_end" if index < 36 else "single_end"
            rows.append(
                {
                    "bundle_id": str(index % 8),
                    "row_id": f"g268_{index:03d}",
                    "experiment_accession": f"SRX{index:08d}",
                    "biosample_accession": f"SAMN{index:08d}",
                    "technical_runs": str(runs),
                    "fastq_files": str(runs * (2 if layout == "paired_end" else 1)),
                    "fastq_bytes": "1",
                    "effective_library_layout": layout,
                    "aggregation_rule": (
                        "matewise_comma_ordered_technical_runs_then_one_participant_RSEM_library"
                        if layout == "paired_end"
                        else "comma_ordered_single_end_technical_runs_then_one_participant_RSEM_library"
                    ),
                }
            )
        residual = 512881099727 - len(rows)
        rows[-1]["fastq_bytes"] = str(1 + residual)
        return rows

    def test_gse268273_accepts_participant_collapsed_run_topology(self) -> None:
        rows = self.gse268273_rows()
        observed = validate_gse268273_topology(GSE268273_PLAN_FIELDS, rows)
        self.assertEqual(observed["technical_runs"], 824)
        self.assertEqual(observed["fastq_files"], 968)

    def test_gse268273_rejects_run_as_participant_or_mate_loss(self) -> None:
        rows = self.gse268273_rows()
        broken = copy.deepcopy(rows)
        broken[0]["fastq_files"] = broken[0]["technical_runs"]
        with self.assertRaises(ExternalCohortReadinessError):
            validate_gse268273_topology(GSE268273_PLAN_FIELDS, broken)


if __name__ == "__main__":
    unittest.main()
