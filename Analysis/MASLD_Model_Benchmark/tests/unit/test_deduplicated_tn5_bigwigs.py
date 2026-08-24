from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

import pyBigWig

from scripts import build_deduplicated_tn5_bigwigs as builder


class DeduplicatedTn5BigWigTests(unittest.TestCase):
    def test_read_support_is_not_used_as_pcr_duplicate_weight(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "fragments.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t0\t10\twell1_AAAA-1\t2\n")
                handle.write("chr1\t4\t10\twell1_CCCC-1\t3\n")
                handle.write("chr2\t5\t20\twell1_GGGG-1\t1\n")
            output = root / "weighted.bw"
            result = builder._build_one(
                {
                    "source": str(source),
                    "output": str(output),
                    "relative_output": "weighted.bw",
                    "chrom_sizes": [(contig, 100) for contig in builder.PRIMARY_CONTIGS],
                    "size_bytes": source.stat().st_size,
                    "sha256": builder.sha256_file(source),
                    "records": "3",
                    "read_support": "6",
                    "donor_id": "1",
                    "outer_fold": "0",
                    "lineage_id": "hepatocyte",
                    "analysis_role": "primary",
                    "nuclei": "3",
                }
            )
            bigwig = pyBigWig.open(str(output))
            try:
                self.assertEqual(
                    bigwig.intervals("chr1"),
                    ((0, 1, 1.0), (4, 5, 1.0), (9, 10, 2.0)),
                )
                self.assertEqual(
                    bigwig.intervals("chr2"), ((5, 6, 1.0), (19, 20, 1.0))
                )
            finally:
                bigwig.close()
            self.assertEqual(result["unique_fragments"], 3)
            self.assertEqual(result["read_support"], 6)
            self.assertEqual(result["tn5_insertions"], 6)
            self.assertEqual(result["nonzero_positions"], 5)
            self.assertEqual(result["max_pending_positions"], 2)

    def test_endpoint_rejects_out_of_bounds_cut_site(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "fragments.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t0\t101\twell1_AAAA-1\t1\n")
            with self.assertRaises(builder.Tn5BigWigError):
                builder._build_one(
                    {
                        "source": str(source),
                        "output": str(root / "bad.bw"),
                        "relative_output": "bad.bw",
                        "chrom_sizes": [
                            (contig, 100) for contig in builder.PRIMARY_CONTIGS
                        ],
                        "size_bytes": source.stat().st_size,
                        "sha256": builder.sha256_file(source),
                        "records": "1",
                        "read_support": "1",
                        "donor_id": "1",
                        "outer_fold": "0",
                        "lineage_id": "hepatocyte",
                        "analysis_role": "primary",
                        "nuclei": "1",
                    }
                )

    def test_coordinate_order_is_required_for_sweep_line_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "fragments.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t10\t20\twell1_AAAA-1\t1\n")
                handle.write("chr1\t5\t15\twell1_CCCC-1\t1\n")
            with self.assertRaises(builder.Tn5BigWigError):
                builder._build_one(
                    {
                        "source": str(source),
                        "output": str(root / "unsorted.bw"),
                        "relative_output": "unsorted.bw",
                        "chrom_sizes": [
                            (contig, 100) for contig in builder.PRIMARY_CONTIGS
                        ],
                        "size_bytes": source.stat().st_size,
                        "sha256": builder.sha256_file(source),
                        "records": "2",
                        "read_support": "2",
                        "donor_id": "1",
                        "outer_fold": "0",
                        "lineage_id": "hepatocyte",
                        "analysis_role": "primary",
                        "nuclei": "2",
                    }
                )


if __name__ == "__main__":
    unittest.main()
