from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pyBigWig

from scripts import build_chrombpnet_development_profiles as profiles


class ChromBPNetDevelopmentProfileTests(unittest.TestCase):
    def test_exact_bigwig_profiles_reproduce_window_sums(self) -> None:
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
            result = profiles._extract_track_profile(
                {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": profiles.sha256_file(path),
                    "donor_id": "1",
                    "outer_fold": 0,
                    "nuclei": 2,
                    "width": 20,
                    "windows": [("chr1", 0, 20), ("chr1", 40, 60)],
                }
            )
            values = result["values"]
            self.assertEqual(values.dtype, np.uint32)
            self.assertEqual(values.shape, (2, 20))
            np.testing.assert_array_equal(
                values.sum(axis=1), np.asarray([5, 4], dtype=np.uint64)
            )
            self.assertEqual(values[0, 5], 2)
            self.assertEqual(values[0, 7], 3)
            self.assertEqual(values[1, 10], 4)


if __name__ == "__main__":
    unittest.main()
