from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.materialize_gse296875_gse244832_source_window_counts import (
    SourceWindowMaterializationError,
    main,
    read_windows,
)


class MaterializeGSE296875GSE244832SourceWindowCountsTests(unittest.TestCase):
    def _windows(self, path: Path, overlap: bool = False) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("window_index", "window_id", "role", "contig", "start", "end"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            index = 0
            for role, contig in (("valid", "chr1"), ("test", "chr1" if overlap else "chr2")):
                for offset in range(16_000):
                    writer.writerow(
                        {
                            "window_index": index,
                            "window_id": f"{role}_{offset}",
                            "role": role,
                            "contig": contig,
                            "start": offset * 1000,
                            "end": offset * 1000 + 1000,
                        }
                    )
                    index += 1

    def test_exact_whole_contig_axis_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "windows.tsv"
            self._windows(path)
            observed = read_windows(path)
            self.assertEqual(len(observed["valid"]), 16_000)
            self.assertEqual(len(observed["test"]), 16_000)

    def test_cross_rotation_contig_overlap_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "windows.tsv"
            self._windows(path, overlap=True)
            with self.assertRaises(SourceWindowMaterializationError):
                read_windows(path)

    def test_cli_dispatches_config_as_config_path(self) -> None:
        with patch(
            "scripts.materialize_gse296875_gse244832_source_window_counts.materialize",
            return_value={"status": "fixture"},
        ) as mocked, patch(
            "sys.argv",
            [
                "materialize",
                "--root",
                "/tmp/root",
                "--config",
                "/tmp/config.json",
                "--output",
                "/tmp/output",
                "--workers",
                "3",
            ],
        ):
            main()
        self.assertEqual(mocked.call_args.kwargs["config_path"], Path("/tmp/config.json"))
        self.assertNotIn("config", mocked.call_args.kwargs)


if __name__ == "__main__":
    unittest.main()
