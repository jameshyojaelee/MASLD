#!/usr/bin/env python3
"""Tests for the independent GSE105127 source-DAG audit."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from masld_bench.hashing import sha256_file
from scripts import audit_gse105127_source_dag_unexposed as audit


def production_topology_fixture() -> tuple[
    list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]
]:
    rows = audit.EXPECTED_ROWS
    read_base, read_extra = divmod(audit.EXPECTED_RNA_READS, rows)
    rna_byte_base, rna_byte_extra = divmod(audit.EXPECTED_RNA_BYTES, rows)
    rrbs_byte_base, rrbs_byte_extra = divmod(audit.EXPECTED_RRBS_BYTES, rows)
    rna = []
    rrbs = []
    bundle_totals = {
        str(bundle): {"participants": set(), "participant_zone_rows": 0, "combined_source_bytes": 0}
        for bundle in range(8)
    }
    row_index = 0
    for participant_index in range(19):
        participant = f"participant-{participant_index:02d}"
        bundle = str(participant_index % 8)
        for zone in ("CV", "IZ", "PP"):
            row_id = f"row-{row_index:02d}"
            read_count = read_base + int(row_index < read_extra)
            rna_bytes = rna_byte_base + int(row_index < rna_byte_extra)
            rrbs_bytes = rrbs_byte_base + int(row_index < rrbs_byte_extra)
            common = {
                "bundle_id": bundle,
                "row_id": row_id,
                "participant_group_id": participant,
                "zone": zone,
                "pairing_topology": "adjacent_section",
            }
            rna.append(
                {
                    **common,
                    "sample_accession": f"GSMR{row_index}",
                    "run_accession": f"SRR{row_index}",
                    "library_layout": "SINGLE",
                    "read_length": "76",
                    "read_count": str(read_count),
                    "base_count": str(76 * read_count),
                    "fastq_url": f"https://example.invalid/{row_id}.fastq.gz",
                    "fastq_bytes": str(rna_bytes),
                    "fastq_md5": f"{row_index + 1:032x}",
                    "relative_fastq_path": f"participants/{row_id}/reads.fastq.gz",
                }
            )
            rrbs.append(
                {
                    **common,
                    "sample_accession": f"GSMM{row_index}",
                    "source_bed_path": f"/source/{row_id}.bed.gz",
                    "source_bed_bytes": str(rrbs_bytes),
                    "source_bed_sha256": f"{row_index + 1:064x}",
                    "source_cytosine_rows": "1",
                    "relative_collapsed_path": f"participants/{row_id}/cpg_counts.hg19.tsv.gz",
                }
            )
            bundle_totals[bundle]["participants"].add(participant)
            bundle_totals[bundle]["participant_zone_rows"] += 1
            bundle_totals[bundle]["combined_source_bytes"] += rna_bytes + rrbs_bytes
            row_index += 1
    bundles = [
        {
            "bundle_id": bundle,
            "participants": str(len(values["participants"])),
            "participant_zone_rows": str(values["participant_zone_rows"]),
            "combined_source_bytes": str(values["combined_source_bytes"]),
        }
        for bundle, values in sorted(bundle_totals.items())
    ]
    return rna, rrbs, bundles


def write_submission(path: Path, *, wrong_dependency: bool = False) -> None:
    required = {
        "21083518": ("data-stage-202", "afterok:21083517", "gse105127_stage_plan_cpu.sbatch"),
        "21083520": ("data-stage-203", "afterok:21083518", "gse105127_stage_rrbs_io.sbatch"),
        "21083522": ("data-stage-205", "afterok:21083518", "gse105127_stage_rsem_reference_cpu.sbatch"),
        "21087896": ("data-stage-220", "afterany:21083521", "gse105127_stage_rna_download_resumable_io.sbatch"),
        "21087897": ("data-stage-207", "afterok:21087896", "gse105127_stage_rna_quant_cpu.sbatch"),
        "21087898": ("data-stage-208", "afterok:21087897", "gse105127_stage_consolidate_cpu.sbatch"),
        "21088014": ("data-stage-214", "GSE105127_CAMPAIGN_ID=9675cd48993468d4", "gse105127_stage_reference_legacy_exact_io.sbatch"),
        "21088066": ("data-stage-206", "afterok:21088014", "gse105127_stage_crosswalk_io.sbatch"),
        "21088068": ("data-stage-209", "afterok:21088066:21087898", "gse105127_stage_finalize_cpu.sbatch"),
    }
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("JobID", "JobName", "State", "ExitCode", "SubmitLine"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for job_id, (name, dependency, wrapper) in required.items():
            if wrong_dependency and job_id == "21087897":
                dependency = "afterok:1"
            writer.writerow(
                {
                    "JobID": job_id,
                    "JobName": name,
                    "State": "COMPLETED",
                    "ExitCode": "0:0",
                    "SubmitLine": f"sbatch --dependency={dependency} {wrapper}",
                }
            )


class GSE105127IndependentSourceDagAuditTests(unittest.TestCase):
    def test_exact_participant_zone_topology_passes(self) -> None:
        rna, rrbs, bundles = production_topology_fixture()
        observed = audit.validate_topology(
            audit.RNA_FIELDS, rna, audit.RRBS_FIELDS, rrbs,
            audit.BUNDLE_FIELDS, bundles,
        )
        self.assertEqual(observed["participants"], 19)
        self.assertEqual(observed["pairing_topology"], "adjacent_section")

    def test_participant_bundle_split_fails_closed(self) -> None:
        rna, rrbs, bundles = production_topology_fixture()
        rna[1]["bundle_id"] = "7"
        rrbs[1]["bundle_id"] = "7"
        with self.assertRaisesRegex(
            audit.GSE105127IndependentAuditError, "participant containment"
        ):
            audit.validate_topology(
                audit.RNA_FIELDS, rna, audit.RRBS_FIELDS, rrbs,
                audit.BUNDLE_FIELDS, bundles,
            )

    def test_read_base_length_mismatch_fails_closed(self) -> None:
        rna, rrbs, bundles = production_topology_fixture()
        rna[0]["base_count"] = str(int(rna[0]["base_count"]) + 1)
        with self.assertRaisesRegex(
            audit.GSE105127IndependentAuditError, "source contract"
        ):
            audit.validate_topology(
                audit.RNA_FIELDS, rna, audit.RRBS_FIELDS, rrbs,
                audit.BUNDLE_FIELDS, bundles,
            )

    def test_legacy_member_inventory_separates_tail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "legacy.fa.gz"
            payload = b">1\nACGT\n"
            compressed_member = gzip.compress(payload, mtime=0)
            tail = b"not-a-gzip-member"
            path.write_bytes(compressed_member + tail)
            observed = audit.inspect_legacy_gzip_member(path)
        self.assertEqual(observed["uncompressed_bytes"], len(payload))
        self.assertEqual(observed["trailing_bytes"], len(tail))
        self.assertFalse(observed["trailing_starts_with_gzip_magic"])
        self.assertEqual(observed["gzip_member_bytes"], len(compressed_member))

    def test_second_gzip_member_is_not_admissible_tail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "members.fa.gz"
            path.write_bytes(gzip.compress(b">1\nA\n", mtime=0) + gzip.compress(b">2\nC\n", mtime=0))
            observed = audit.inspect_legacy_gzip_member(path)
        self.assertTrue(observed["trailing_starts_with_gzip_magic"])

    def test_modeling_call_is_rejected_from_source_code(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.py"
            path.write_text("model.fit(values)\n", encoding="utf-8")
            with self.assertRaisesRegex(
                audit.GSE105127IndependentAuditError, "modeling call"
            ):
                audit.audit_python_surface(path, sha256_file(path))

    def test_submission_dependency_contract_passes_and_rejects_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid = root / "valid.tsv"
            invalid = root / "invalid.tsv"
            write_submission(valid)
            write_submission(invalid, wrong_dependency=True)
            self.assertEqual(audit.validate_submission_record(valid), sha256_file(valid))
            with self.assertRaisesRegex(
                audit.GSE105127IndependentAuditError, "dependency differs"
            ):
                audit.validate_submission_record(invalid)


if __name__ == "__main__":
    unittest.main()
