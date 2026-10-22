#!/usr/bin/env python3
"""Small complete genetics replays, isolated from all scientific outputs."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE))
import atac_context_v3_lib as contracts


def load_aggregate():
    path = Path(os.environ.get("ATAC_AGGREGATE_SOURCE",
                               PACKAGE / "09_aggregate_genetic_context.py"))
    spec = importlib.util.spec_from_file_location("aggregate_fixture", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_rows(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or tuple(rows[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(root, fasta_path):
    """Interleaved signal families; mapped, excluded and invalid-reference alleles."""
    rows = []
    for study, gene, pair, position, posterior, unique, chrom, allele in (
        ("a", "ENSA", 1, 1, .50, "TRUE", "chr1", "A"),
        ("b", "ENSB", 1, 3, .70, "TRUE", "chr1", "G"),
        ("a", "ENSA", 2, 5, .30, "TRUE", "chr1", "A"),
        ("a", "ENSA", 1, 4, .25, "TRUE", "chr1", "T"),
        ("b", "ENSB", 1, 1, .20, "TRUE", "chrUn", "A"),
        ("a", "ENSA", 2, 6, .40, "TRUE", "chr1", "C"),
        ("a", "ENSA", 1, 8, .25, "FALSE", "chr1", "T"),
        ("b", "ENSB", 1, 7, .10, "TRUE", "chr1", "A"),
        ("a", "ENSA", 2, 8, .30, "FALSE", "chr1", "T"),
    ):
        rows.append(dict(gwas_name=study, gene=gene, ensembl=gene,
                         signal_pair_index=str(pair), snp=f"{study}-{pair}-{position}",
                         **{"SNP.PP.H4": str(posterior)}, unique_liftover=unique,
                         hg38_chrom=chrom, hg38_position_1based=str(position),
                         allele1=allele,
                         allele2="T" if (study == "b" and position == 7)
                         else ("C" if allele != "C" else "T")))
    write_rows(root / "genetics/liftover/variant_liftover_raw.tsv.gz", rows)
    peaks = []
    for lineage in contracts.LINEAGES:
        peaks.extend([
            dict(lineage=lineage, chrom="chr1", start0=0, end=4,
                 supported_gse244832="TRUE", supported_gse281367="FALSE"),
            dict(lineage=lineage, chrom="chr1", start0=2, end=6,
                 supported_gse244832="FALSE", supported_gse281367="TRUE"),
        ])
    write_rows(root / "consensus_peak_manifest.tsv", peaks)
    write_rows(root / "input_manifest.tsv", [dict(
        role="hg38_fasta", relative_or_absolute_path=str(fasta_path),
        sha256=contracts.sha256_file(fasta_path),
    )])
    write_rows(root / "genetics/prepared/replay_plan.tsv", [
        dict(gwas_name=study, ensembl=gene, trait=study, trait_class=trait_class,
             ancestry="EUR")
        for study, gene, trait_class in (("a", "ENSA", "direct_MASLD"),
                                        ("b", "ENSB", "liver_enzyme"))
    ])
    write_rows(root / "genetics/replay_execution/batch_1/exports/study/chr1/gene/signal_pairs.tsv", [
        dict(gwas_name=study, gene=gene, ensembl=gene,
             signal_pair_index=str(pair), **{"PP.H4.abf": str(pp)},
             gwas_signal=signal, eqtl_signal="e1")
        for study, gene, pair, pp, signal in (("a", "ENSA", 1, .9, "g2"),
                                              ("a", "ENSA", 2, .9, "g1"),
                                              ("b", "ENSB", 1, .8, "g1"))
    ])
    return rows


def run_fixture(module, root):
    with patch.object(module, "arguments", return_value=SimpleNamespace(candidate_root=root)), \
         patch.object(module, "candidate_root", side_effect=lambda value: value):
        module.main()


class GeneticAggregationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "candidate"
        self.fasta = Path(self.temporary.name) / "reference.fa"
        self.fasta.write_text(">chr1\nACGTACGTACGT\n")
        # The real FASTA reader is used, including its index and coordinates.
        import pysam
        pysam.faidx(str(self.fasta))
        self.rows = make_fixture(self.root, self.fasta)
        self.module = load_aggregate()

    def test_complete_replay_retains_order_masses_and_primary_tie_break(self):
        run_fixture(self.module, self.root)
        out = self.root / "genetics/context"
        audit = self.module.read_rows(out / "variant_liftover_audit.tsv")
        self.assertEqual([r["hg19_variant_id"] for r in audit],
                         [r["snp"] for r in self.rows])
        self.assertEqual([int(r["posterior_row_id"]) for r in audit], list(range(1, 10)))
        self.assertEqual(audit[4]["mapping_status"], "nonstandard_chromosome")
        self.assertEqual(audit[6]["mapping_status"], "unmapped_or_multimapped")
        self.assertEqual(audit[7]["mapping_status"], "hg38_allele_validation_failed")
        all_pairs = self.module.read_rows(out / "genetic_lineage_context_all_pairs.tsv")
        self.assertEqual(len(all_pairs), 3 * len(contracts.LINEAGES))
        first = all_pairs[0]
        self.assertEqual((first["gwas_name"], first["signal_pair_index"]), ("a", "1"))
        for field, expected in (("mapped_posterior_mass", .75),
                                ("gse244832_accessible_mass", .75),
                                ("gse281367_accessible_mass", .25),
                                ("shared_accessible_mass", .25),
                                ("any_accessible_mass", .75)):
            self.assertEqual(float(first[field]), expected)
        primary = self.module.read_rows(out / "genetic_lineage_context_primary_pairs.tsv")
        self.assertEqual(len(primary), 2 * len(contracts.LINEAGES))
        self.assertEqual({r["signal_pair_index"] for r in primary if r["gwas_name"] == "a"}, {"2"})
        self.assertEqual({r["evidence_state"] for r in primary}, {"untestable"})

    def test_nonadjacent_duplicate_is_rejected_without_publishing(self):
        self.rows.append(dict(self.rows[0]))
        self.assert_rejected("duplicated variant posterior row")

    def test_invalid_probability_is_rejected_without_publishing(self):
        self.rows[-1]["SNP.PP.H4"] = "nan"
        self.assert_rejected("invalid probability")

    def test_incomplete_posterior_mass_is_rejected_without_publishing(self):
        self.rows[-1]["SNP.PP.H4"] = "0.01"
        self.assert_rejected("does not sum to one")

    def test_empty_input_is_rejected_without_publishing(self):
        write_rows(self.root / "genetics/liftover/variant_liftover_raw.tsv.gz",
                   [], tuple(self.rows[0]))
        with self.assertRaisesRegex(contracts.ContractError, "empty"):
            run_fixture(self.module, self.root)
        self.assertFalse((self.root / "genetics/context").exists())

    def test_existing_output_is_preserved(self):
        out = self.root / "genetics/context"
        out.mkdir()
        marker = out / "existing.tsv"
        marker.write_bytes(b"unchanged\n")
        with self.assertRaisesRegex(contracts.ContractError, "refusing to overwrite"):
            run_fixture(self.module, self.root)
        self.assertEqual(marker.read_bytes(), b"unchanged\n")

    def assert_rejected(self, message):
        write_rows(self.root / "genetics/liftover/variant_liftover_raw.tsv.gz", self.rows)
        with self.assertRaisesRegex(contracts.ContractError, message):
            run_fixture(self.module, self.root)
        self.assertFalse((self.root / "genetics/context").exists())


if __name__ == "__main__":
    unittest.main()
