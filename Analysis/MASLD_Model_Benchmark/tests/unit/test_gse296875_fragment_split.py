from __future__ import annotations

import gzip
import io
from pathlib import Path
import tempfile
import unittest

from scripts import split_gse296875_donor_lineage_fragments as splitter


class GSE296875FragmentSplitTests(unittest.TestCase):
    def test_split_well_accepts_cellranger_headers_and_preserves_read_support(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            membership = root / "well1.tsv.gz"
            with gzip.open(membership, "wt", encoding="utf-8") as handle:
                handle.write("\t".join(splitter.MEMBERSHIP_FIELDS) + "\n")
                handle.write(
                    "well1\tAAAA-1\twell1_AAAA-1\t1\tHepatocytes\t"
                    "hepatocyte\tprimary\t0\n"
                )
            fragment = root / "well1.fragments.tsv.gz"
            with gzip.open(fragment, "wt", encoding="utf-8") as handle:
                handle.write("# pipeline_name=cellranger-arc\n")
                handle.write("# pipeline_version=cellranger-arc-2.0.0\n")
                handle.write("# reference_version=2020-A\n")
                handle.write("chr1\t10\t30\tAAAA-1\t3\n")
            result = splitter._split_well(
                {
                    "well": "well1",
                    "fragment": str(fragment),
                    "membership": str(membership),
                    "output": str(root / "split"),
                    "expected_size": str(fragment.stat().st_size),
                    "expected_sha256": splitter.sha256_file(fragment),
                }
            )
            self.assertEqual(result["header_lines"], 3)
            self.assertEqual(result["selected_records"], 1)
            self.assertEqual(result["selected_read_support"], 3)
            output = root / "split" / "donor_1__hepatocyte.tsv.gz"
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.read(), "chr1\t10\t30\twell1_AAAA-1\t3\n"
                )

    def test_hashing_reader_covers_exact_compressed_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t1\t2\tAAAA-1\t1\n")
            raw = path.open("rb")
            hashed = splitter.HashingReader(raw)
            compressed = gzip.GzipFile(fileobj=hashed, mode="rb")
            text = io.TextIOWrapper(compressed, encoding="utf-8")
            self.assertEqual(text.read(), "chr1\t1\t2\tAAAA-1\t1\n")
            text.close()
            raw.close()
            self.assertEqual(hashed.digest.hexdigest(), splitter.sha256_file(path))

    def test_deterministic_gzip_preserves_records_and_read_support(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chunks = []
            rows = (
                ("chr1\t5\t25\twell1_AAAA-1\t2\n",),
                (
                    "chr1\t1\t11\twell2_CCCC-1\t3\n",
                    "chr2\t7\t17\twell2_CCCC-1\t1\n",
                ),
            )
            for index, values in enumerate(rows):
                path = root / f"chunk{index}.tsv.gz"
                text, hashed, raw = splitter._open_deterministic_gzip(path)
                for value in values:
                    text.write(value)
                splitter._close_deterministic_gzip(text, hashed, raw)
                chunks.append(path)
            output = root / "merged.tsv.gz"
            result = splitter._merge_group(
                {
                    "chunks": [str(path) for path in chunks],
                    "output": str(output),
                    "relative_output": "fragments/donor_1/hepatocyte.fragments.tsv.gz",
                    "expected_records": 3,
                    "expected_read_support": 6,
                    "donor_id": "1",
                    "outer_fold": "0",
                    "lineage_id": "hepatocyte",
                    "analysis_role": "primary",
                    "nuclei": "2",
                    "wells": "well1,well2",
                }
            )
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                observed = handle.readlines()
            self.assertEqual([line.split("\t")[:3] for line in observed], [
                ["chr1", "1", "11"],
                ["chr1", "5", "25"],
                ["chr2", "7", "17"],
            ])
            self.assertEqual(result["records"], 3)
            self.assertEqual(result["read_support"], 6)
            self.assertEqual(result["sha256"], splitter.sha256_file(output))

    def test_line_parser_rejects_nonprimary_contig(self) -> None:
        with self.assertRaises(splitter.FragmentSplitError):
            splitter._line_key("chrM\t1\t2\tcell\t1\n")

    def test_cellranger_lexicographic_contig_order_is_registered(self) -> None:
        chr10, _fields = splitter._line_key("chr10\t1\t2\tcell\t1\n")
        chr2, _fields = splitter._line_key("chr2\t1\t2\tcell\t1\n")
        self.assertLess(chr10, chr2)


if __name__ == "__main__":
    unittest.main()
