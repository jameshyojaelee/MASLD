from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.build_gse256398_gencode49_crosswalk import (
    GSE256398CrosswalkError,
    build_crosswalk,
    parse_gencode_genes,
    stable_gene_id,
)


class GSE256398GENCODE49CrosswalkTests(unittest.TestCase):
    def test_stable_gene_id_removes_only_version(self) -> None:
        self.assertEqual(stable_gene_id("ENSG000001.12"), "ENSG000001")
        self.assertEqual(stable_gene_id("ENSG000001_PAR_Y.3"), "ENSG000001_PAR_Y")

    def test_build_crosswalk_masks_absent_features(self) -> None:
        current = {
            "ENSG1": {
                "stable_id": "ENSG1",
                "versioned_id": "ENSG1.4",
                "gene_name": "A",
                "gene_type": "protein_coding",
                "contig": "chr1",
                "start_1based": "10",
                "end_1based": "20",
                "strand": "+",
            }
        }
        rows, states = build_crosswalk(["ENSG1.1", "ENSG2.2"], ["A", "OLD"], current)
        self.assertEqual(rows[0]["allowed_project_input"], "true")
        self.assertEqual(rows[1]["allowed_project_input"], "false")
        self.assertEqual(states["stable_id_exact_gencode_v49"], 1)
        self.assertEqual(states["absent_or_retired_from_gencode_v49"], 1)

    def test_duplicate_stable_source_ids_fail(self) -> None:
        with self.assertRaisesRegex(GSE256398CrosswalkError, "stable gene IDs"):
            build_crosswalk(["ENSG1.1", "ENSG1.2"], ["A", "A"], {})

    def test_gtf_parser_reads_gene_records_only(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "genes.gtf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("##gencode-version 49\n")
                handle.write(
                    'chr1\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "ENSG1.4"; '
                    'gene_type "protein_coding"; gene_name "A";\n'
                )
                handle.write(
                    'chr1\tHAVANA\ttranscript\t10\t20\t.\t+\t.\tgene_id "ENSG1.4";\n'
                )
            genes = parse_gencode_genes(path)
            self.assertEqual(set(genes), {"ENSG1"})
            self.assertEqual(genes["ENSG1"]["gene_name"], "A")


if __name__ == "__main__":
    unittest.main()
