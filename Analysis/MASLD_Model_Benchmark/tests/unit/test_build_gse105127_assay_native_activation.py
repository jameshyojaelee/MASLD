from __future__ import annotations

import csv
import io
from pathlib import Path
import tomllib
import unittest

from scripts.build_gse105127_assay_native_activation import (
    EXPECTED_RNA_FASTQ_BYTES,
    EXPECTED_RNA_READS,
    GSE105127ActivationError,
    canonical_cpg_interval,
    collapse_cpg_rows,
    liftover_state,
    methylated_read_count,
    parse_ena_runinfo,
    participant_group_id,
    row_id,
)


class GSE105127AssayNativeActivationTests(unittest.TestCase):
    def test_plus_and_shifted_minus_rows_collapse_once(self) -> None:
        observed = collapse_cpg_rows(
            [
                {"start0": 10488, "end0": 10489, "percent": 100.0, "coverage": 9, "strand": "+"},
                {"start0": 10489, "end0": 10490, "percent": 86.36, "coverage": 22, "strand": "-"},
            ]
        )
        self.assertEqual((observed["start0"], observed["end0"]), (10488, 10490))
        self.assertEqual(observed["methylated_reads"], 28)
        self.assertEqual(observed["coverage"], 31)
        self.assertEqual(observed["strand_state"], "both_strands")
        self.assertFalse(observed["auxiliary_columns_used"])

    def test_one_strand_site_is_retained_with_coverage(self) -> None:
        observed = collapse_cpg_rows(
            [{"start0": 20, "end0": 21, "percent": 25.0, "coverage": 4, "strand": "+"}]
        )
        self.assertEqual(observed["methylated_reads"], 1)
        self.assertEqual(observed["coverage"], 4)
        self.assertEqual(observed["strand_state"], "+_only")

    def test_negative_strand_underflow_and_duplicate_strand_fail_closed(self) -> None:
        with self.assertRaises(GSE105127ActivationError):
            canonical_cpg_interval(0, 1, "-")
        with self.assertRaises(GSE105127ActivationError):
            collapse_cpg_rows(
                [
                    {"start0": 10, "end0": 11, "percent": 50.0, "coverage": 2, "strand": "+"},
                    {"start0": 10, "end0": 11, "percent": 50.0, "coverage": 2, "strand": "+"},
                ]
            )

    def test_percentage_must_recover_an_integer_count(self) -> None:
        self.assertEqual(methylated_read_count(86.36, 22), 19)
        with self.assertRaises(GSE105127ActivationError):
            methylated_read_count(33.33, 2)

    def test_liftover_dispositions_are_fail_closed(self) -> None:
        base = {"source_is_cpg": True, "source_contig_recognized": True}
        self.assertEqual(
            liftover_state(
                source_is_cpg=True,
                source_contig_recognized=False,
                mapped_intervals=[],
            ),
            "source_contig_unrecognized",
        )
        self.assertEqual(
            liftover_state(
                source_is_cpg=False,
                source_contig_recognized=True,
                mapped_intervals=[],
            ),
            "source_not_cpg",
        )
        self.assertEqual(liftover_state(**base, mapped_intervals=[]), "unmapped_chain")
        self.assertEqual(
            liftover_state(**base, mapped_intervals=[("chr1", 1, 3), ("chr2", 2, 4)]),
            "mapped_nonunique",
        )
        self.assertEqual(
            liftover_state(**base, mapped_intervals=[("chr1", 1, 4)]),
            "mapped_length_changed",
        )
        self.assertEqual(
            liftover_state(**base, mapped_intervals=[("chrUn", 1, 3)]),
            "mapped_nonprimary",
        )
        self.assertEqual(
            liftover_state(
                **base,
                mapped_intervals=[("chr1", 1, 3)],
                target_primary=True,
            ),
            "target_not_cpg",
        )
        self.assertEqual(
            liftover_state(
                **base,
                mapped_intervals=[("chr1", 1, 3)],
                target_primary=True,
                target_is_cpg=True,
                roundtrip_exact=True,
            ),
            "mapped_unique_cpg",
        )
        self.assertEqual(
            liftover_state(
                **base,
                mapped_intervals=[("chr1", 1, 3)],
                target_primary=True,
                target_is_cpg=True,
                roundtrip_exact=False,
            ),
            "roundtrip_failed",
        )

    def test_row_ids_group_zone_rows_without_exposing_outcomes(self) -> None:
        self.assertEqual(row_id("6610", "CV"), row_id("6610", "CV"))
        self.assertNotEqual(row_id("6610", "CV"), row_id("6610", "IZ"))
        self.assertNotIn("6610", row_id("6610", "CV"))
        self.assertEqual(participant_group_id("6610"), participant_group_id("6610"))
        self.assertNotIn("6610", participant_group_id("6610"))

    def test_ena_parser_rejects_non_single_end_or_wrong_read_length(self) -> None:
        fields = [
            "run_accession", "study_accession", "sample_accession", "experiment_accession",
            "scientific_name", "library_strategy", "library_source", "library_selection",
            "library_layout", "instrument_platform", "instrument_model", "read_count",
            "base_count", "fastq_ftp", "fastq_bytes", "fastq_md5",
        ]
        rows = []
        for index in range(114):
            rna = index < 57
            rows.append(
                {
                    "run_accession": f"SRR{index}",
                    "study_accession": "PRJNA414905",
                    "sample_accession": f"SAMN{index}",
                    "experiment_accession": f"SRX{index}",
                    "scientific_name": "Homo sapiens",
                    "library_strategy": "RNA-Seq" if rna else "Bisulfite-Seq",
                    "library_source": "TRANSCRIPTOMIC" if rna else "GENOMIC",
                    "library_selection": "cDNA" if rna else "Reduced Representation",
                    "library_layout": "SINGLE",
                    "instrument_platform": "ILLUMINA",
                    "instrument_model": "Illumina HiSeq 2500",
                    "read_count": "1",
                    "base_count": "76" if rna else "100",
                    "fastq_ftp": f"ftp.sra.ebi.ac.uk/x/SRR{index}.fastq.gz",
                    "fastq_bytes": "1",
                    "fastq_md5": "0" * 32,
                }
            )
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        with self.assertRaisesRegex(GSE105127ActivationError, "RNA FASTQ byte/read/base census"):
            parse_ena_runinfo(buffer.getvalue().encode())

    def test_ena_parser_accepts_the_frozen_assay_native_census(self) -> None:
        fields = [
            "run_accession", "study_accession", "sample_accession", "experiment_accession",
            "scientific_name", "library_strategy", "library_source", "library_selection",
            "library_layout", "instrument_platform", "instrument_model", "read_count",
            "base_count", "fastq_ftp", "fastq_bytes", "fastq_md5",
        ]
        rows = []
        for index in range(114):
            rna = index < 57
            reads = (
                1
                if rna and index < 56
                else EXPECTED_RNA_READS - 56
                if rna
                else 100
            )
            size = (
                1
                if rna and index < 56
                else EXPECTED_RNA_FASTQ_BYTES - 56
                if rna
                else 100
            )
            rows.append(
                {
                    "run_accession": f"SRR{index}",
                    "study_accession": "PRJNA414905",
                    "sample_accession": f"SAMN{index}",
                    "experiment_accession": f"SRX{index}",
                    "scientific_name": "Homo sapiens",
                    "library_strategy": "RNA-Seq" if rna else "Bisulfite-Seq",
                    "library_source": "TRANSCRIPTOMIC" if rna else "GENOMIC",
                    "library_selection": "cDNA" if rna else "Reduced Representation",
                    "library_layout": "SINGLE",
                    "instrument_platform": "ILLUMINA",
                    "instrument_model": "Illumina HiSeq 2500",
                    "read_count": str(reads),
                    "base_count": str(reads * (76 if rna else 100)),
                    "fastq_ftp": f"ftp.sra.ebi.ac.uk/x/SRR{index}.fastq.gz",
                    "fastq_bytes": str(size),
                    "fastq_md5": "0" * 32,
                }
            )
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        observed = parse_ena_runinfo(buffer.getvalue().encode())
        self.assertEqual(len(observed), 114)
        self.assertEqual(sum(row["library_strategy"] == "RNA-Seq" for row in observed), 57)

    def test_taskspec_preserves_adjacent_section_and_assay_native_target(self) -> None:
        task = tomllib.loads(
            Path("config/evaluation/gse105127_adjacent_section_rna_rrbs_task.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(task["unit_of_inference"], "participant")
        self.assertEqual(task["required_pairing_levels"], ["adjacent_section"])
        self.assertEqual(task["target_modalities"], ["rrbs_methylated_reads", "rrbs_coverage"])
        joined = " ".join(task["admission_gates"]).lower()
        self.assertIn("missing", task["missingness_policy"].lower())
        self.assertIn("never contribute counts", joined)
        self.assertIn("cannot support an rna-only claim", joined)
        self.assertNotIn("same_section", task["required_pairing_levels"])


if __name__ == "__main__":
    unittest.main()
