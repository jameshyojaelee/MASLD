from __future__ import annotations

import unittest

from scripts.build_corgi_tile_smoke_subset import select_fold


class CorgiTileSmokeSubsetTests(unittest.TestCase):
    def test_selection_is_deterministic_and_covers_strata(self) -> None:
        tiles = []
        windows = []
        for index in range(12):
            tile_id = f"t{index:02d}"
            contig = "chr1" if index < 6 else "chr2"
            cls = "PLS" if index % 2 == 0 else "dELS"
            tiles.append(
                {
                    "tile_id": tile_id,
                    "contig": contig,
                    "central_start": str(index * 1000),
                }
            )
            windows.append({"tile_id": tile_id, "ccre_class": cls})
        first = select_fold(tiles, windows, 8)
        second = select_fold(tiles, windows, 8)
        self.assertEqual(first, second)
        strata = {
            (row["contig"], next(item["ccre_class"] for item in windows if item["tile_id"] == row["tile_id"]))
            for row in first
        }
        self.assertEqual(strata, {("chr1", "PLS"), ("chr1", "dELS"), ("chr2", "PLS"), ("chr2", "dELS")})


if __name__ == "__main__":
    unittest.main()
