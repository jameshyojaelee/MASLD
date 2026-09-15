from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
import random
import tempfile
import unittest

import numpy as np

from scripts.build_gse268273_rsem_reference import (
    EXPECTED_PRIMARY_SCAFFOLD_SOURCE_FEATURES,
    gtf_gene_contract,
)
from scripts.download_gse268273_raw_bundle import (
    GSE268273DownloadError,
    validate_member,
)
from scripts.plan_gse268273_raw_campaign import (
    GENCODE_V36_GTF_BYTES,
    GENCODE_V36_GTF_MD5,
    GENCODE_V36_GTF_URL,
    GSE268273RawPlanError,
    assign_participant_bundles,
    expand_run_manifest,
)
from scripts.quantify_gse268273_rsem_bundle import (
    aggregate_v49_expected_counts,
    arrange_participant_fastqs,
    parse_rsem_gene_results,
)


def run_row(
    row_id: str,
    run: str,
    experiment: str,
    biosample: str,
    sizes: tuple[int, ...],
) -> dict[str, str]:
    names = (
        (f"{run}.fastq.gz",)
        if len(sizes) == 1
        else (f"{run}_1.fastq.gz", f"{run}_2.fastq.gz")
    )
    paths = tuple(
        f"ftp.sra.ebi.ac.uk/vol1/fastq/{run}/{name}" for name in names
    )
    return {
        "row_id": row_id,
        "run_accession": run,
        "experiment_accession": experiment,
        "biosample_accession": biosample,
        "library_layout": "SINGLE",
        "fastq_ftp": ";".join(paths),
        "fastq_bytes": ";".join(str(value) for value in sizes),
        "fastq_md5": ";".join(f"{index + 1:032x}" for index in range(len(sizes))),
        "replicate_role": "technical_run_partition",
    }


def ena_read_row(row: dict[str, str], read_count: int = 10) -> dict[str, str]:
    files = len(row["fastq_ftp"].split(";"))
    return {
        "run_accession": row["run_accession"],
        "experiment_accession": row["experiment_accession"],
        "library_layout": row["library_layout"],
        "read_count": str(read_count),
        "base_count": str(read_count * (101 if files == 1 else 202)),
        "fastq_ftp": row["fastq_ftp"],
        "fastq_bytes": row["fastq_bytes"],
        "fastq_md5": row["fastq_md5"],
    }


class GSE268273RawCampaignTests(unittest.TestCase):
    def test_source_annotation_is_the_v36_primary_assembly(self) -> None:
        self.assertTrue(
            GENCODE_V36_GTF_URL.endswith(
                "gencode.v36.primary_assembly.annotation.gtf.gz"
            )
        )
        self.assertEqual(GENCODE_V36_GTF_BYTES, 44_517_553)
        self.assertEqual(GENCODE_V36_GTF_MD5, "5038acd7158686c69a1f88a5e26a5934")
        self.assertEqual(len(EXPECTED_PRIMARY_SCAFFOLD_SOURCE_FEATURES), 7)

    def test_declared_single_two_files_are_admitted_as_paired_mates(self) -> None:
        source = run_row("p1", "SRR1", "SRX1", "SAMN1", (10, 20))
        observed = expand_run_manifest(
            [source],
            [ena_read_row(source)],
            production_contract=False,
        )
        self.assertEqual(len(observed), 2)
        self.assertEqual({row["file_part_count"] for row in observed}, {2})
        self.assertEqual(
            {row["read_role"] for row in observed},
            {"paired_end_mate_1", "paired_end_mate_2"},
        )
        self.assertEqual(
            {row["effective_library_layout"] for row in observed}, {"paired_end"}
        )

    def test_participant_lpt_bundles_are_deterministic_and_indivisible(self) -> None:
        rows = [
            run_row("p1", "SRR1", "SRX1", "SAMN1", (100,)),
            run_row("p2", "SRR2", "SRX2", "SAMN2", (70, 20)),
            run_row("p3", "SRR3", "SRX3", "SAMN3", (60,)),
            run_row("p4", "SRR4", "SRX4", "SAMN4", (50,)),
        ]
        expanded = expand_run_manifest(
            rows, [ena_read_row(row) for row in rows], production_contract=False
        )
        first, _, _ = assign_participant_bundles(expanded, 2)
        random.Random(17).shuffle(expanded)
        second, _, _ = assign_participant_bundles(expanded, 2)
        self.assertEqual(first, second)
        by_participant = {}
        for row in first:
            by_participant.setdefault(row["row_id"], set()).add(row["bundle_id"])
        self.assertTrue(all(len(value) == 1 for value in by_participant.values()))

    def test_raw_plan_rejects_outcome_columns(self) -> None:
        row = run_row("p1", "SRR1", "SRX1", "SAMN1", (10,))
        ena = ena_read_row(row)
        row["fibrosis"] = "F3"
        with self.assertRaises(GSE268273RawPlanError):
            expand_run_manifest([row], [ena], production_contract=False)

    def test_paired_mates_are_not_concatenated_as_single_end(self) -> None:
        rows = [
            run_row("p1", "SRR1", "SRX1", "SAMN1", (10, 20)),
            run_row("p1", "SRR2", "SRX1", "SAMN1", (30, 40)),
        ]
        expanded = expand_run_manifest(
            rows, [ena_read_row(row) for row in rows], production_contract=False
        )
        layout, upstream, downstream = arrange_participant_fastqs(expanded)
        self.assertEqual(layout, "paired_end")
        self.assertEqual(
            [row["run_accession"] for row in upstream], ["SRR1", "SRR2"]
        )
        self.assertEqual(
            [row["read_role"] for row in upstream],
            ["paired_end_mate_1", "paired_end_mate_1"],
        )
        self.assertEqual(
            [row["read_role"] for row in downstream],
            ["paired_end_mate_2", "paired_end_mate_2"],
        )

    def test_unresolved_declared_single_topology_fails_closed(self) -> None:
        source = run_row("p1", "SRR1", "SRX1", "SAMN1", (10, 20))
        ena = ena_read_row(source)
        ena["base_count"] = "1010"
        with self.assertRaises(GSE268273RawPlanError):
            expand_run_manifest([source], [ena], production_contract=False)

    def test_fastq_admission_checks_md5_size_and_gzip(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"@r1\nACGT\n+\n!!!!\n")
            md5 = hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()
            validate_member(path, path.stat().st_size, md5)
            with self.assertRaises(GSE268273DownloadError):
                validate_member(path, path.stat().st_size + 1, md5)

    def test_fastq_admission_checks_first_record_mate_role(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mate2.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"@SRR1.1 instrument/2\nACGT\n+\n!!!!\n")
            md5 = hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()
            validate_member(
                path, path.stat().st_size, md5, "paired_end_mate_2"
            )
            with self.assertRaises(GSE268273DownloadError):
                validate_member(
                    path, path.stat().st_size, md5, "paired_end_mate_1"
                )

    def test_reference_gtf_and_v49_raw_sum_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            gtf = Path(temporary) / "genes.gtf.gz"
            with gzip.open(gtf, "wt", encoding="utf-8") as handle:
                handle.write(
                    'chr1\tGENCODE\tgene\t1\t4\t.\t+\t.\tgene_id "ENSG1.1";\n'
                )
                handle.write(
                    'chr1\tGENCODE\tgene\t7\t9\t.\t+\t.\tgene_id "ENSG2.2";\n'
                )
            exact, stable, contigs = gtf_gene_contract(gtf)
        self.assertEqual(exact, {"ENSG1.1", "ENSG2.2"})
        self.assertEqual(stable, {"ENSG1", "ENSG2"})
        self.assertEqual(set(contigs), {"chr1"})
        candidates = [
            {
                "source_feature_count": "1",
                "source_feature_ids": "ENSG1.1",
                "raw_count_aggregation": "identity",
            },
            {
                "source_feature_count": "2",
                "source_feature_ids": "ENSG2.1;ENSG2.2",
                "raw_count_aggregation": "deterministic_sum_before_normalization",
            },
        ]
        observed = aggregate_v49_expected_counts(
            {"ENSG1.1": 2.5, "ENSG2.1": 3.25, "ENSG2.2": 4.75}, candidates
        )
        np.testing.assert_array_equal(observed, np.asarray([2.5, 8.0]))

    def test_rsem_result_parser_preserves_fractional_expected_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sample.genes.results"
            path.write_text(
                "gene_id\ttranscript_id(s)\tlength\teffective_length\t"
                "expected_count\tTPM\tFPKM\n"
                "ENSG1.1\tENST1.1\t100\t80\t1.25\t2\t3\n",
                encoding="utf-8",
            )
            observed = parse_rsem_gene_results(path)
        self.assertEqual(observed, {"ENSG1.1": 1.25})


if __name__ == "__main__":
    unittest.main()
