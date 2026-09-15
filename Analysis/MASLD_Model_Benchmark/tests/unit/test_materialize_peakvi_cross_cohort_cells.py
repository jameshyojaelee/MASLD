from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

from scripts.materialize_peakvi_cross_cohort_cells import (
    LINEAGES,
    PeakVICellMaterializationError,
    read_windows,
    selected_cells,
    windows_at,
)


class PeakVICrossCohortCellMaterializationTests(unittest.TestCase):
    def test_selection_is_hash_deterministic_and_unit_balanced(self) -> None:
        rows = []
        for donor in range(1, 13):
            for lineage in LINEAGES:
                for cell in range(10):
                    rows.append(
                        {
                            "cell_id": f"Z{donor:02d}_{lineage}_{cell}",
                            "raw_barcode": f"bc{cell}",
                            "well_id": f"Z{donor:02d}",
                            "donor_id": f"Z{donor:02d}",
                            "lineage_id": lineage,
                            "analysis_role": "primary",
                            "outer_fold": str(donor % 5),
                        }
                    )
        first, _ = selected_cells(rows, dataset_id="gse281367", cap=8, seed=17)
        second, _ = selected_cells(list(reversed(rows)), dataset_id="gse281367", cap=8, seed=17)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 12 * 4 * 8)
        self.assertEqual(len({row["cell_hash"] for row in first}), len(first))

    def test_window_roles_are_contig_disjoint_and_overlap_lookup_is_exact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "windows.tsv"
            with path.open("w", encoding="utf-8") as handle:
                handle.write("window_index\twindow_id\trole\tcontig\tstart\tend\n")
                for role, contig, start_index in (("valid", "chr1", 0), ("test", "chr2", 16000)):
                    for offset in range(16_000):
                        start = offset * 2000
                        handle.write(
                            f"{start_index + offset}\t{role}_{offset}\t{role}\t{contig}\t{start}\t{start + 1000}\n"
                        )
            windows = read_windows(path)
            self.assertEqual(list(windows_at(windows["valid"], "chr1", 500)), [0])
            self.assertEqual(list(windows_at(windows["valid"], "chr2", 500)), [])

    def test_target_hash_includes_donor_when_barcodes_recur(self) -> None:
        rows = []
        for donor in range(1, 13):
            for lineage in LINEAGES:
                for cell in range(8):
                    rows.append(
                        {
                            "cell_id": f"shared-{lineage}-{cell}",
                            "raw_barcode": f"shared-{lineage}-{cell}",
                            "well_id": f"Z{donor:02d}",
                            "donor_id": f"Z{donor:02d}",
                            "lineage_id": lineage,
                            "analysis_role": "primary",
                            "outer_fold": str(donor % 5),
                        }
                    )
        selected, _ = selected_cells(
            rows, dataset_id="gse281367", cap=8, seed=17
        )
        self.assertEqual(len(selected), 12 * 4 * 8)
        self.assertEqual(len({row["cell_hash"] for row in selected}), len(selected))

    def test_source_requires_all_39_by_four_units(self) -> None:
        with self.assertRaises(PeakVICellMaterializationError):
            selected_cells([], dataset_id="gse296875", cap=8, seed=0)


if __name__ == "__main__":
    unittest.main()
