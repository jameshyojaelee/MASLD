from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from scripts import emit_training_tagalign as emitter


class TrainingTagAlignTests(unittest.TestCase):
    def test_held_genomic_fold_is_excluded_without_read_support_weighting(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t2\t8\tcell_a\t4\n")
                handle.write("chr2\t3\t10\tcell_b\t7\n")
                handle.write("chr3\t4\t12\tcell_c\t9\n")
            output = root / "output.tagAlign.gz"
            result = emitter.emit(
                source=source,
                source_size_bytes=source.stat().st_size,
                source_sha256=emitter.sha256_file(source),
                fold_by_contig={"chr1": 0, "chr2": 1, "chr3": 2},
                genomic_test_fold=2,
                genomic_valid_fold=1,
                output=output,
            )
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.readlines(),
                    [
                        "chr1\t2\t3\tcell_a/L\t1000\t+\n",
                        "chr1\t7\t8\tcell_a/R\t1000\t-\n",
                    ],
                )
            self.assertEqual(result["included_unique_fragments"], 1)
            self.assertEqual(result["included_tn5_events"], 2)
            self.assertEqual(result["excluded_genomic_test_fragments"], 1)
            self.assertEqual(result["excluded_genomic_valid_fragments"], 1)
            self.assertEqual(result["allowed_genomic_folds"], [0, 3, 4])
            self.assertFalse(result["read_support_used_as_signal_weight"])

    def test_training_tagaligns_are_partitioned_by_allowed_genomic_fold(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                for fold in range(5):
                    handle.write(
                        f"chr{fold + 1}\t{fold + 2}\t{fold + 8}\tcell_{fold}\t{fold + 4}\n"
                    )
            output = root / "by_fold"
            result = emitter.emit_by_genomic_fold(
                source=source,
                source_size_bytes=source.stat().st_size,
                source_sha256=emitter.sha256_file(source),
                fold_by_contig={f"chr{fold + 1}": fold for fold in range(5)},
                genomic_test_fold=2,
                genomic_valid_fold=1,
                output_directory=output,
            )
            self.assertEqual(result["allowed_genomic_folds"], [0, 3, 4])
            self.assertEqual(
                [row["unique_fragments"] for row in result["outputs"]], [1, 1, 1]
            )
            self.assertEqual(result["included_tn5_events"], 6)
            self.assertEqual(result["excluded_genomic_test_fragments"], 1)
            self.assertEqual(result["excluded_genomic_valid_fragments"], 1)
            for row in result["outputs"]:
                with gzip.open(output / row["path"], "rt", encoding="utf-8") as handle:
                    lines = handle.readlines()
                self.assertEqual(len(lines), 2)
                self.assertEqual(lines[0].split("\t")[0], f"chr{row['genomic_fold'] + 1}")
                self.assertTrue(row["sha256"])


if __name__ == "__main__":
    unittest.main()
