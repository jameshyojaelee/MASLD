from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

import pyBigWig

from scripts import build_training_fold_pseudobulk as builder


class TrainingFoldPseudobulkTests(unittest.TestCase):
    def test_pool_is_sorted_and_uses_bed_exclusive_right_cut(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = []
            for index, rows in enumerate(
                (
                    ("chr1\t5\t10\ta\t3\n",),
                    ("chr1\t1\t8\tb\t2\n", "chr2\t4\t7\tb\t1\n"),
                )
            ):
                path = root / f"source{index}.tsv.gz"
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.writelines(rows)
                sources.append(path)
            fragment = root / "pooled.fragments.tsv.gz"
            bigwig_path = root / "pooled.tn5.bw"
            result = builder.pool_one(
                {
                    "sources": [str(path) for path in sources],
                    "fragment_output": str(fragment),
                    "bigwig_output": str(bigwig_path),
                    "relative_fragment": fragment.name,
                    "relative_bigwig": bigwig_path.name,
                    "chrom_sizes": [
                        (contig, 100) for contig in builder.PRIMARY_CONTIGS
                    ],
                    "donor_test_fold": 0,
                    "donor_valid_fold": 1,
                    "donor_train_folds": "2,3,4",
                    "lineage_id": "hepatocyte",
                    "donors": 2,
                    "nuclei": 2,
                    "expected_records": 3,
                    "expected_read_support": 6,
                }
            )
            with gzip.open(fragment, "rt", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.readlines(),
                    [
                        "chr1\t1\t8\tb\t2\n",
                        "chr1\t5\t10\ta\t3\n",
                        "chr2\t4\t7\tb\t1\n",
                    ],
                )
            bigwig = pyBigWig.open(str(bigwig_path))
            try:
                self.assertEqual(
                    bigwig.intervals("chr1"),
                    ((1, 2, 1.0), (5, 6, 1.0), (7, 8, 1.0), (9, 10, 1.0)),
                )
                self.assertEqual(
                    bigwig.intervals("chr2"), ((4, 5, 1.0), (6, 7, 1.0))
                )
            finally:
                bigwig.close()
            self.assertEqual(result["unique_fragments"], 3)
            self.assertEqual(result["read_support"], 6)
            self.assertEqual(result["tn5_insertions"], 6)
            self.assertEqual(result["max_pending_positions"], 3)


if __name__ == "__main__":
    unittest.main()
