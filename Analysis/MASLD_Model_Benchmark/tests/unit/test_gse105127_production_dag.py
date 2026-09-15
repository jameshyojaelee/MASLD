from __future__ import annotations

import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import sqlite3
import tempfile
import unittest

from scripts.build_gse105127_cpg_crosswalk import (
    PENDING_CHAIN,
    PENDING_TARGET_CHECK,
    classify_forward,
    classify_reverse,
    export_crosswalk,
)
from scripts.collapse_gse105127_rrbs_bundle import (
    GSE105127CollapseError,
    collapse_one,
)
from scripts.download_gse105127_rna_bundle import (
    GSE105127DownloadError,
    audit_fastq,
)
from scripts.plan_gse105127_production_dag import assign_participant_bundles
from scripts.quantify_gse105127_rna_bundle import (
    GSE105127QuantificationError,
    audit_gene_results,
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


class GSE105127ProductionDAGTests(unittest.TestCase):
    def test_bundle_assignment_never_splits_three_participant_zones(self) -> None:
        rows = [
            {
                "participant_group_id": f"opaque_{participant:02d}",
                "combined_source_bytes": str((participant + 1) * (zone + 1)),
            }
            for participant in range(19)
            for zone in range(3)
        ]
        assignment, summary = assign_participant_bundles(rows, 8)
        self.assertEqual(len(assignment), 19)
        self.assertEqual(sum(item["participant_zone_rows"] for item in summary), 57)
        self.assertTrue(all(item["participant_zone_rows"] == 3 * item["participants"] for item in summary))

    def test_streamed_bissnp_collapse_uses_only_percent_and_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.bed.gz"
            lines = [
                'track name=fixture type=bedDetail description="CG"',
                "chr1\t10488\t10489\t100.00\t9\t+\t10488\t10489\t0,0,0\t7\t22",
                "chr1\t10489\t10490\t86.36\t22\t-\t10489\t10490\t0,0,0\t8\t9",
                "chr1\t10500\t10501\t25.00\t4\t+\t10500\t10501\t0,0,0\t9\t99",
                "chr2\t200\t201\t50.00\t2\t+\t200\t201\t0,0,0\t10\t0",
                "chr2\t201\t202\t50.00\t2\t-\t201\t202\t0,0,0\t11\t0",
            ]
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write("\n".join(lines) + "\n")
            row = {
                "row_id": "opaque_row",
                "participant_group_id": "opaque_participant",
                "zone": "CV",
                "source_bed_path": str(source),
                "source_bed_bytes": str(source.stat().st_size),
                "source_bed_sha256": file_sha256(source),
                "source_cytosine_rows": "5",
            }
            receipt = collapse_one(row, root / "final", root / "attempt")
            self.assertEqual(receipt["source_rows"], 5)
            self.assertEqual(receipt["canonical_cpgs"], 3)
            self.assertEqual(receipt["reciprocal_aux_pairs"], 1)
            self.assertEqual(receipt["auxiliary_inconsistent_pairs"], 1)
            self.assertEqual(receipt["detail_nonzero_rows"], 5)
            self.assertFalse(receipt["auxiliary_columns_used_for_measurement"])
            with gzip.open(
                root / "final/participants/opaque_row/cpg_counts.hg19.tsv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                observed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                (observed[0]["source_start0"], observed[0]["methylated_reads"], observed[0]["coverage"]),
                ("10488", "28", "31"),
            )
            self.assertEqual(observed[1]["strand_state"], "plus_only")
            self.assertEqual(observed[1]["auxiliary_state"], "unpaired_aux_nonzero")
            self.assertEqual(observed[2]["auxiliary_state"], "auxiliary_inconsistent_ignored")
            self.assertEqual(collapse_one(row, root / "final", root / "replay"), receipt)
            changed = dict(row, source_bed_sha256="0" * 64)
            with self.assertRaisesRegex(GSE105127CollapseError, "replayed RRBS member"):
                collapse_one(changed, root / "final", root / "wrong-replay")

    def test_fastq_audit_binds_record_read_and_base_census(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n")
            self.assertEqual(
                audit_fastq(path, expected_reads=2, expected_bases=8, read_length=4),
                {"read_count": 2, "base_count": 8, "read_length": 4},
            )
            with self.assertRaisesRegex(GSE105127DownloadError, "read/base census"):
                audit_fastq(path, expected_reads=3, expected_bases=12, read_length=4)

    def test_crosswalk_preserves_unmapped_as_distinct_from_roundtrip_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = root / "crosswalk.sqlite"
            connection = sqlite3.connect(database)
            connection.executescript(
                "CREATE TABLE sites (id INTEGER PRIMARY KEY, source_chrom TEXT NOT NULL, source_start INTEGER NOT NULL, source_end INTEGER NOT NULL, state TEXT NOT NULL, target_chrom TEXT, target_start INTEGER, target_end INTEGER);"
                "CREATE TABLE forward (id INTEGER PRIMARY KEY, mappings INTEGER NOT NULL, chrom TEXT, start INTEGER, end INTEGER);"
                "CREATE TABLE reverse (id INTEGER PRIMARY KEY, mappings INTEGER NOT NULL, chrom TEXT, start INTEGER, end INTEGER);"
            )
            connection.executemany(
                "INSERT INTO sites VALUES (?,?,?,?,?,NULL,NULL,NULL)",
                [
                    (1, "chr1", 10, 12, PENDING_CHAIN),
                    (2, "chr1", 20, 22, PENDING_CHAIN),
                    (3, "chr1", 30, 32, PENDING_CHAIN),
                    (4, "chr1", 40, 42, PENDING_CHAIN),
                    (5, "chr1", 50, 52, PENDING_CHAIN),
                    (6, "chr1", 60, 62, PENDING_CHAIN),
                    (7, "chrUn", 70, 72, "source_contig_unrecognized"),
                    (8, "chr1", 80, 82, "source_not_cpg"),
                    (9, "chr1", 90, 92, "strand_pair_inconsistent"),
                    (10, "chr1", 100, 102, PENDING_CHAIN),
                ],
            )
            connection.executemany(
                "INSERT INTO forward VALUES (?,?,?,?,?)",
                [
                    (2, 2, "chr1", 120, 122),
                    (3, 1, "chr1", 130, 133),
                    (4, 1, "chrUn", 140, 142),
                    (5, 1, "chr1", 150, 152),
                    (6, 1, "chr1", 160, 162),
                    (10, 1, "chr1", 200, 202),
                ],
            )
            connection.commit()
            connection.close()
            classify_forward(database, root / "target.bed")
            connection = sqlite3.connect(database)
            connection.execute("UPDATE sites SET state='target_not_cpg' WHERE id=5 AND state=?", (PENDING_TARGET_CHECK,))
            connection.executemany(
                "INSERT INTO reverse VALUES (?,?,?,?,?)",
                [(6, 1, "chr1", 61, 63), (10, 1, "chr1", 100, 102)],
            )
            connection.commit()
            connection.close()
            classify_reverse(database)
            counts = export_crosswalk(database, root / "crosswalk.tsv.gz")
            self.assertEqual(counts["unmapped_chain"], 1)
            self.assertEqual(counts["mapped_nonunique"], 1)
            self.assertEqual(counts["mapped_length_changed"], 1)
            self.assertEqual(counts["mapped_nonprimary"], 1)
            self.assertEqual(counts["target_not_cpg"], 1)
            self.assertEqual(counts["roundtrip_failed"], 1)
            self.assertEqual(counts["mapped_unique_cpg"], 1)
            self.assertEqual(sum(counts.values()), 10)

    def test_rsem_audit_preserves_continuous_values_but_accepts_integer_library(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rsem.genes.results"
            path.write_text(
                "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n"
                "g1\tt1\t100\t80\t1\t250000\t1\n"
                "g2\tt2\t200\t180\t2\t750000\t2\n",
                encoding="utf-8",
            )
            observed = audit_gene_results(path)
            self.assertEqual(observed["fractional_expected_count_rows"], 0)
            path.write_text(path.read_text(encoding="utf-8").replace("750000", "700000"), encoding="utf-8")
            with self.assertRaisesRegex(GSE105127QuantificationError, "TPM census"):
                audit_gene_results(path)

    def test_production_planner_never_opens_historical_outcome_table(self) -> None:
        source = Path("scripts/plan_gse105127_production_dag.py").read_text(encoding="utf-8")
        self.assertNotIn('rrbs_file_qc.tsv', source)
        self.assertIn('rrbs_root / "ARTIFACTS.json"', source)

    def test_slurm_dag_is_cpu_io_only_generic_and_matches_resource_ceiling(self) -> None:
        paths = sorted(Path("slurm").glob("gse105127_stage_*.sbatch"))
        self.assertEqual(len(paths), 9)
        cpu_hours = 0
        memory_gb_hours = 0
        tasks = 0
        for path in paths:
            source = path.read_text(encoding="utf-8")
            header = {
                key: value
                for key, value in re.findall(r"^#SBATCH --([^=]+)=(.+)$", source, re.MULTILINE)
            }
            self.assertIn(header["partition"], {"cpu", "io"})
            self.assertEqual(header["qos"], "nslab")
            self.assertEqual(header["account"], "nslab")
            self.assertRegex(header["job-name"], r"^data-stage-[0-9]{3}$")
            self.assertNotRegex(header["job-name"].lower(), r"masld|gse|model|rna|rrbs")
            count = 8 if header.get("array") == "0-7" else 1
            hours = int(header["time"].split(":", 1)[0])
            cpus = int(header["cpus-per-task"])
            memory = int(header["mem"].removesuffix("G"))
            cpu_hours += count * cpus * hours
            memory_gb_hours += count * memory * hours
            tasks += count
        self.assertEqual(tasks, 30)
        self.assertEqual(cpu_hours, 5_336)
        self.assertEqual(memory_gb_hours, 44_192)
        submit = Path("slurm/submit_gse105127_production_dag.sh").read_text(encoding="utf-8")
        self.assertIn('afterok:${REFERENCE}:${RRBS}', submit)
        self.assertIn('afterok:${RNA_DOWNLOAD}:${RSEM_REFERENCE}', submit)
        self.assertIn('afterok:${CROSSWALK}:${CONSOLIDATE}', submit)


if __name__ == "__main__":
    unittest.main()
