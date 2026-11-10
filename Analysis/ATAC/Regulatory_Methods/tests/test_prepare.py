import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("prepare_analysis", Path(__file__).resolve().parents[1] / "scripts/prepare_analysis.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CoordinateTests(unittest.TestCase):
    def test_gene_tss_and_promoters_preserve_strand(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "test.gtf.gz"
            with gzip.open(path, "wt") as handle:
                handle.write('chr1\ttest\tgene\t101\t300\t.\t+\t.\tgene_id "ENSG1.2"; gene_name "A";\n')
                handle.write('chr2\ttest\tgene\t4001\t5000\t.\t-\t.\tgene_id "ENSG2.1"; gene_name "B";\n')
            rows = list(module.gene_records(path))
            self.assertEqual([r["tss0"] for r in rows], [100, 4999])
            self.assertEqual((rows[0]["promoter_start0"], rows[0]["promoter_end"]), (0, 2101))
            self.assertEqual((rows[1]["promoter_start0"], rows[1]["promoter_end"]), (2999, 7000))
            self.assertEqual(rows[0]["gene_id"], "ENSG1")


if __name__ == "__main__":
    unittest.main()
