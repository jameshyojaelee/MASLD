#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
import sys
import tempfile
import unittest


TEST_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TEST_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from gene_identity import (  # noqa: E402
    APPROVED_CLASSES,
    CLASS_ANTISENSE_BODY,
    CLASS_ANTISENSE_EXONIC,
    CLASS_COMPLEX,
    CLASS_INTERGENIC,
    CLASS_OPPOSITE_PROMOTER,
    CLASS_SAME_BODY,
    CLASS_SAME_EXONIC,
    GENE_IDENTITY_FIELDS,
    IdentityBuildError,
    RelationEvidence,
    build_release,
    classify_relation,
    parse_attributes,
)


def attrs(gene_id: str, gene_type: str, gene_name: str) -> str:
    return (
        f'gene_id "{gene_id}"; gene_type "{gene_type}"; '
        f'gene_name "{gene_name}"; level 2;'
    )


def feature(
    chromosome: str,
    kind: str,
    start: int,
    end: int,
    strand: str,
    gene_id: str,
    gene_type: str,
    gene_name: str,
) -> str:
    return "\t".join(
        (
            chromosome,
            "TEST",
            kind,
            str(start),
            str(end),
            ".",
            strand,
            ".",
            attrs(gene_id, gene_type, gene_name),
        )
    )


def add_gene(
    lines: list[str],
    chromosome: str,
    start: int,
    end: int,
    strand: str,
    gene_id: str,
    gene_type: str,
    gene_name: str,
    exons: list[tuple[int, int]],
) -> None:
    lines.append(
        feature(
            chromosome,
            "gene",
            start,
            end,
            strand,
            gene_id,
            gene_type,
            gene_name,
        )
    )
    for exon_start, exon_end in exons:
        lines.append(
            feature(
                chromosome,
                "exon",
                exon_start,
                exon_end,
                strand,
                gene_id,
                gene_type,
                gene_name,
            )
        )


class GeneIdentityTests(unittest.TestCase):
    def test_attribute_parser_preserves_versioned_identifiers(self):
        observed = parse_attributes(
            'gene_id "ENSG00000000001.7"; gene_type "lncRNA"; gene_name "LNC_A";'
        )
        self.assertEqual(observed["gene_id"], "ENSG00000000001.7")
        self.assertEqual(observed["gene_type"], "lncRNA")

    def test_priority_and_complex_conflict_are_deterministic(self):
        cases = (
            (
                RelationEvidence(
                    antisense_exonic=frozenset({"a"}),
                    antisense_body=frozenset({"b"}),
                ),
                CLASS_ANTISENSE_EXONIC,
            ),
            (
                RelationEvidence(same_strand_exonic=frozenset({"a"})),
                CLASS_SAME_EXONIC,
            ),
            (
                RelationEvidence(antisense_body=frozenset({"a"})),
                CLASS_ANTISENSE_BODY,
            ),
            (
                RelationEvidence(same_strand_body=frozenset({"a"})),
                CLASS_SAME_BODY,
            ),
            (
                RelationEvidence(opposite_promoter=frozenset({"a"})),
                CLASS_OPPOSITE_PROMOTER,
            ),
            (RelationEvidence(), CLASS_INTERGENIC),
            (
                RelationEvidence(
                    antisense_exonic=frozenset({"a"}),
                    same_strand_exonic=frozenset({"b"}),
                ),
                CLASS_COMPLEX,
            ),
        )
        self.assertEqual({expected for _, expected in cases}, set(APPROVED_CLASSES))
        for evidence, expected in cases:
            self.assertEqual(classify_relation(evidence)[0], expected)

    def test_full_fixture_builds_all_classes_and_ambiguity_flags(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = root / "fixture.gtf"
            output = root / "candidate"
            lines = ["##description: synthetic GENCODE fixture"]

            add_gene(
                lines,
                "chr1",
                100,
                300,
                "+",
                "ENSG00000000001.1",
                "protein_coding",
                "PC_A",
                [(100, 150), (250, 300)],
            )
            add_gene(
                lines,
                "chr1",
                500,
                700,
                "-",
                "ENSG00000000002.2",
                "protein_coding",
                "DUP",
                [(500, 550), (650, 700)],
            )
            add_gene(
                lines,
                "chr1",
                900,
                1100,
                "+",
                "ENSG00000000003.3",
                "protein_coding",
                "DUP",
                [(900, 950), (1050, 1100)],
            )
            add_gene(
                lines,
                "chr1",
                5000,
                5100,
                "+",
                "ENSG00000000004.1",
                "protein_coding",
                "PC_PROM",
                [(5000, 5100)],
            )

            add_gene(
                lines,
                "chr1",
                120,
                140,
                "-",
                "ENSG00000000101.1",
                "lncRNA",
                "L_ANTI_EXON",
                [(120, 140)],
            )
            add_gene(
                lines,
                "chr1",
                260,
                280,
                "+",
                "ENSG00000000102.1",
                "lncRNA",
                "L_SAME_EXON",
                [(260, 280)],
            )
            add_gene(
                lines,
                "chr1",
                180,
                220,
                "-",
                "ENSG00000000103.1",
                "lncRNA",
                "L_ANTI_BODY",
                [(180, 220)],
            )
            add_gene(
                lines,
                "chr1",
                180,
                220,
                "+",
                "ENSG00000000104.1",
                "lncRNA",
                "L_SAME_BODY",
                [(180, 220)],
            )
            add_gene(
                lines,
                "chr1",
                3200,
                3300,
                "-",
                "ENSG00000000105.1",
                "lncRNA",
                "L_PROM",
                [(3200, 3300)],
            )
            add_gene(
                lines,
                "chr1",
                8000,
                8100,
                "+",
                "ENSG00000000106.1",
                "lncRNA",
                "L_INTER",
                [(8000, 8100)],
            )
            add_gene(
                lines,
                "chr1",
                520,
                930,
                "+",
                "ENSG00000000107.1",
                "lncRNA",
                "L_COMPLEX",
                [(520, 530), (920, 930)],
            )
            add_gene(
                lines,
                "chrUn_KI270442v1",
                9000,
                9100,
                "+",
                "ENSG00000000108.1",
                "lncRNA",
                "L_ALT",
                [(9000, 9100)],
            )
            gtf.write_text("\n".join(lines) + "\n", encoding="utf-8")

            os.environ["MASLD_NONCODING_FIXTURE_MODE"] = "1"
            self.addCleanup(os.environ.pop, "MASLD_NONCODING_FIXTURE_MODE", None)
            census = build_release(gtf, output, fixture_mode=True)
            self.assertEqual(census["n_lncrna_genes"], 8)

            with (output / "gene_identity.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                identity = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(tuple(identity[0]), GENE_IDENTITY_FIELDS)
            by_name_identity = {row["gene_name"]: row for row in identity}
            self.assertEqual(
                by_name_identity["PC_A"]["gene_id_versioned"], "ENSG00000000001.1"
            )
            self.assertEqual(
                by_name_identity["PC_A"]["gene_id_base"], "ENSG00000000001"
            )
            self.assertEqual(
                by_name_identity["DUP"]["symbol_mapping_status"], "duplicated"
            )
            self.assertEqual(
                by_name_identity["L_ALT"]["is_canonical_chromosome"], "false"
            )

            with (output / "lncrna_genomic_class.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                classes = list(csv.DictReader(handle, delimiter="\t"))
            by_name = {row["gene_name"]: row for row in classes}
            expected = {
                "L_ANTI_EXON": CLASS_ANTISENSE_EXONIC,
                "L_SAME_EXON": CLASS_SAME_EXONIC,
                "L_ANTI_BODY": CLASS_ANTISENSE_BODY,
                "L_SAME_BODY": CLASS_SAME_BODY,
                "L_PROM": CLASS_OPPOSITE_PROMOTER,
                "L_INTER": CLASS_INTERGENIC,
                "L_COMPLEX": CLASS_COMPLEX,
                "L_ALT": CLASS_INTERGENIC,
            }
            self.assertEqual(
                {name: by_name[name]["lncrna_genomic_class"] for name in expected},
                expected,
            )
            self.assertEqual(
                by_name["L_SAME_EXON"]["mapping_status"], "mapping_ambiguous"
            )
            self.assertEqual(
                by_name["L_COMPLEX"]["mapping_status"], "mapping_ambiguous"
            )
            self.assertEqual(
                by_name["L_ANTI_BODY"]["mapping_status"], "mapping_unambiguous"
            )
            self.assertEqual(by_name["L_ALT"]["main_text_eligible"], "false")

            with (output / "lncrna_mapping_exclusions.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                exclusions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {row["gene_name"] for row in exclusions},
                {"L_SAME_EXON", "L_COMPLEX", "L_ALT"},
            )

            with (output / "checksum_manifest.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                checksums = list(csv.DictReader(handle, delimiter="\t"))
            for row in checksums:
                path = output / row["relative_path"]
                self.assertEqual(int(row["size_bytes"]), path.stat().st_size)
                self.assertEqual(
                    row["sha256"], hashlib.sha256(path.read_bytes()).hexdigest()
                )

    def test_source_gate_and_nonoverwriting_output_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = root / "fixture.gtf"
            gtf.write_text("## empty fixture\n", encoding="utf-8")
            with self.assertRaisesRegex(
                IdentityBuildError, "production annotation release"
            ):
                build_release(
                    gtf,
                    root / "wrong_release",
                    annotation_release="GENCODE v48",
                )
            with self.assertRaisesRegex(IdentityBuildError, "byte count mismatch"):
                build_release(gtf, root / "candidate")
            (root / "existing").mkdir()
            os.environ["MASLD_NONCODING_FIXTURE_MODE"] = "1"
            self.addCleanup(os.environ.pop, "MASLD_NONCODING_FIXTURE_MODE", None)
            with self.assertRaisesRegex(IdentityBuildError, "already exists"):
                build_release(gtf, root / "existing", fixture_mode=True)


if __name__ == "__main__":
    unittest.main()
