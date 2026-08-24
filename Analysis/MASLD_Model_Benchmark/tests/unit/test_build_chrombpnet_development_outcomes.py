from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pyBigWig

from scripts import build_chrombpnet_development_outcomes as outcomes


class ChromBPNetDevelopmentOutcomeTests(unittest.TestCase):
    def test_exact_bigwig_window_sums_are_uint32(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "track.bw"
            bigwig = pyBigWig.open(str(path), "w")
            bigwig.addHeader([("chr1", 100)])
            bigwig.addEntries(
                ["chr1", "chr1", "chr1"],
                [5, 7, 50],
                ends=[6, 8, 51],
                values=[2.0, 3.0, 4.0],
            )
            bigwig.close()
            result = outcomes._extract_track(
                {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": outcomes.sha256_file(path),
                    "donor_id": "1",
                    "outer_fold": 0,
                    "nuclei": 2,
                    "unique_fragments": 4,
                    "windows": [("chr1", 0, 10), ("chr1", 40, 60)],
                }
            )
            np.testing.assert_array_equal(
                result["values"], np.asarray([5, 4], dtype=np.uint32)
            )


if __name__ == "__main__":
    unittest.main()
