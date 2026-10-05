#!/usr/bin/env python3
"""Synthetic execution evidence tests, without scanning or counting real BAMs."""
import copy
import gzip
import json
from pathlib import Path
import tempfile
import unittest
import subprocess
import sys
from count_verified_fragments import digest, file_stat, validate_count_execution


class ExecutionEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gtf = self.root / "grouped.gtf.gz"
        with gzip.open(self.gtf, "wt") as handle:
            handle.write('chr1\tfixture\texon\t1\t10\t.\t+\t.\tgene_id "GENE1";\n')
        self.correspondence = self.root / "annotation_correspondence.tsv"
        self.correspondence.write_text("counting_gene_id\nGENE1\n")
        self.summary = self.root / "resolution_summary.json"
        self.summary.write_text(json.dumps({"annotation_resolution_complete": True,
            "n_unresolved_annotations": 0, "grouped_counting_annotation": str(self.gtf),
            "grouped_counting_annotation_sha256": digest(self.gtf),
            "annotation_correspondence_sha256": digest(self.correspondence)}))
        self.source = self.root / "SRR1.original.bam"
        self.source.write_bytes(b"synthetic stat fixture; not a BAM")
        self.sorted = self.root / "name-sorted.bam"
        self.sorted.write_bytes(b"synthetic assignment stat fixture; not a BAM")
        self.evidence = self.root / "execution.json"
        self.value = {"status": "complete", "sample_id": "SRR1", "layout": "SINGLE",
            "featurecounts_returncode": 0, "samtools_sort_returncode": 0,
            "annotation_summary_sha256": digest(self.summary),
            "grouped_counting_annotation_sha256": digest(self.gtf),
            "annotation_correspondence_sha256": digest(self.correspondence),
            "source_bam_before": file_stat(self.source), "source_bam_after": file_stat(self.source),
            "assignment_bam": file_stat(self.sorted),
            "featurecounts_argv": ["/fixture/featureCounts", "-T", "1", "-s", "2", "-M", "-R", "BAM",
                "-a", str(self.gtf), "-o", str(self.root / "not_analysis.tsv"), "--Rpath", str(self.root), str(self.source)],
            "samtools_sort_argv": ["/fixture/samtools", "sort", "-N", "-o", str(self.sorted),
                str(self.root / (self.source.name + ".featureCounts.bam"))],
            "versions": {name: "synthetic fixture" for name in ("python", "samtools", "featureCounts")}}

    def check(self, value=None, layout="SINGLE"):
        self.evidence.write_text(json.dumps(self.value if value is None else value))
        return validate_count_execution(self.summary, self.evidence, self.sorted, "SRR1", layout)

    def test_exact_annotation_command_and_source_stats_pass(self):
        self.assertEqual(self.check()["grouped_counting_annotation_sha256"], digest(self.gtf))

    def test_paired_requires_complete_fragment_flags(self):
        value = copy.deepcopy(self.value)
        value["layout"] = "PAIRED"
        with self.assertRaises(ValueError): self.check(value, "PAIRED")
        value["featurecounts_argv"][1:1] = ["-p", "--countReadPairs", "-B"]
        self.check(value, "PAIRED")

    def test_wrong_annotation_argument_rejected(self):
        value = copy.deepcopy(self.value)
        value["featurecounts_argv"][value["featurecounts_argv"].index("-a") + 1] = str(self.correspondence)
        with self.assertRaises(ValueError): self.check(value)

    def test_undeclared_eligibility_flag_rejected(self):
        value = copy.deepcopy(self.value)
        value["featurecounts_argv"][1:1] = ["-Q", "20"]
        with self.assertRaises(ValueError): self.check(value)

    def test_missing_multimapper_or_BAM_assignment_flag_rejected(self):
        for flag in ("-M", "-R"):
            value = copy.deepcopy(self.value)
            value["featurecounts_argv"].remove(flag)
            with self.assertRaises(ValueError): self.check(value)

    def test_source_modified_after_assignment_rejected(self):
        self.source.write_bytes(b"different fixture")
        with self.assertRaises(ValueError): self.check()

    def test_source_before_after_disagreement_rejected(self):
        value = copy.deepcopy(self.value)
        value["source_bam_after"]["mtime_ns"] = "0"
        with self.assertRaises(ValueError): self.check(value)

    def test_correspondence_modified_rejected(self):
        self.correspondence.write_text("counting_gene_id\nGENE2\n")
        with self.assertRaises(ValueError): self.check()

    def test_name_sort_output_must_be_counted_BAM(self):
        value = copy.deepcopy(self.value)
        value["samtools_sort_argv"][value["samtools_sort_argv"].index("-o") + 1] = "/fixture/other.bam"
        with self.assertRaises(ValueError): self.check(value)

    def test_missing_executable_versions_rejected(self):
        value = copy.deepcopy(self.value)
        del value["versions"]["featureCounts"]
        with self.assertRaises(ValueError): self.check(value)

    def test_tiny_reader_and_matrix_certificate_link_execution(self):
        import pysam
        header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
        with pysam.AlignmentFile(str(self.sorted), "wb", header=header) as bam:
            read = pysam.AlignedSegment()
            read.query_name = "fixture1"
            read.query_sequence = "AAAAAAAAAA"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 255
            read.cigar = [(0, 10)]
            read.set_tags([("NH", 1), ("HI", 1), ("XS", "Assigned"), ("XN", 1), ("XT", "GENE1")])
            bam.write(read)
        self.value["assignment_bam"] = file_stat(self.sorted)
        self.check()
        script = Path(__file__).parent / "count_verified_fragments.py"
        count_root = self.root / "counts"
        output = count_root / "SRR1"
        subprocess.run([sys.executable, "-B", str(script), "--annotation-summary", str(self.summary),
            "--assignment-bam", str(self.sorted), "--execution-evidence", str(self.evidence),
            "--sample-id", "SRR1", "--layout", "SINGLE", "--output", str(output)],
            check=True, capture_output=True, text=True)
        certificate = json.loads((output / "count_completion.json").read_text())
        self.assertEqual(certificate["n_assigned_molecules"], 1)
        self.assertEqual(certificate["assignment_execution_sha256"], digest(output / "assignment_execution.json"))
        roster = self.root / "roster.tsv"
        roster.write_text("sample_id\tanalysis_unit_id\tunit_evidence\nSRR1\tfixture_unit\tsource_labels_plus_stable_filtered_genetic_identity_cluster_not_verified_person\n")
        assembled = self.root / "assembled"
        rscript = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
        subprocess.run([rscript, "--vanilla", str(script.parent / "assemble_verified_counts.R"),
            str(self.summary), str(roster), str(count_root), str(assembled)],
            check=True, capture_output=True, text=True)
        completion = json.loads((assembled / "count_completion.json").read_text())
        self.assertEqual(completion["n_samples"], 1)
        self.assertEqual(completion["sample_count_evidence_sha256"], digest(assembled / "sample_count_evidence.tsv"))
        copied = output / "assignment_execution.json"
        altered = json.loads(copied.read_text())
        altered["featurecounts_argv"][altered["featurecounts_argv"].index("-s") + 1] = "1"
        copied.write_text(json.dumps(altered))
        certificate["assignment_execution_sha256"] = digest(copied)
        (output / "count_completion.json").write_text(json.dumps(certificate))
        rejected = subprocess.run([rscript, "--vanilla", str(script.parent / "assemble_verified_counts.R"),
            str(self.summary), str(roster), str(count_root), str(self.root / "rejected")],
            capture_output=True, text=True)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("Actual assignment command/source evidence differs", rejected.stderr)


if __name__ == "__main__": unittest.main()
