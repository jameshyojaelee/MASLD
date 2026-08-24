from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.build_corgi_outcome_aligned_tiles import CCRE_FIELDS, build


def write_reference(path: Path, name: str, sequence: str, width: int = 80) -> None:
    with path.open("wb") as handle:
        handle.write(f">{name}\n".encode())
        offset = len(name) + 2
        for start in range(0, len(sequence), width):
            handle.write(sequence[start : start + width].encode() + b"\n")
    lines = (len(sequence) + width - 1) // width
    del lines
    Path(str(path) + ".fai").write_text(
        f"{name}\t{len(sequence)}\t{offset}\t{width}\t{width + 1}\n"
    )


class CorgiOutcomeAlignedTileTests(unittest.TestCase):
    def test_boundary_exclusion_and_greedy_unique_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            split = root / "split"
            split.mkdir()
            (split / "ARTIFACTS.json").write_text("{}\n")
            rows = []
            starts = (10, 70_000, 470_000)
            for fold in range(5):
                for index, start in enumerate(starts):
                    rows.append(
                        {
                            "contig": "chr1",
                            "output_start": str(start),
                            "output_end": str(start + 1_000),
                            "window_id": f"f{fold}w{index}",
                            "genomic_fold": str(fold),
                            "window_class": "encode_ccre",
                            "ccre_class": "PLS",
                            "ccre_id": f"c{fold}_{index}",
                            "ccre_start": str(start + 100),
                            "ccre_end": str(start + 900),
                            "input_start": str(start - 557),
                            "input_end": str(start + 1557),
                            "selection_hash": f"{fold:02d}{index:02d}".ljust(64, "0"),
                        }
                    )
            with (split / "ccre_evaluation_windows.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=CCRE_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            fasta = root / "reference.fa"
            write_reference(fasta, "chr1", "ACGT" * 300_000)
            output = root / "output"
            first = build(
                split_contract=split,
                reference_fasta=fasta,
                reference_fai=Path(str(fasta) + ".fai"),
                output=output,
                expected_windows_per_fold=3,
            )
            self.assertEqual(first["scoreable_windows"], 10)
            self.assertEqual(first["excluded_windows"], 5)
            with (output / "window_map.tsv").open(newline="") as handle:
                mapped = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len({row["window_id"] for row in mapped}), 10)
            self.assertTrue(all(int(row["past_last_output_bin"]) <= 6_144 for row in mapped))
            with gzip.open(output / "tiles.fa.gz", "rt") as handle:
                sequences = [line.strip() for line in handle if not line.startswith(">")]
            self.assertTrue(sequences)
            self.assertEqual(set("".join(sequences)), set("ACGT"))


if __name__ == "__main__":
    unittest.main()
