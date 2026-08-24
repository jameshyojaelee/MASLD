from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


MODULE = Path(__file__).parents[2] / "scripts" / "enformer_crested_build_track_crosswalk.py"
SPEC = importlib.util.spec_from_file_location("enformer_crested_crosswalk_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
crosswalk = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(crosswalk)


class EnformerCREstedCrosswalkTests(unittest.TestCase):
    def test_assay_families_are_native(self) -> None:
        self.assertEqual(crosswalk._assay_family("DNASE:liver"), "accessibility")
        self.assertEqual(crosswalk._assay_family("ATAC:hepatocyte"), "accessibility")
        self.assertEqual(crosswalk._assay_family("CAGE:liver"), "cage")
        self.assertEqual(crosswalk._assay_family("CHIP:HepG2"), "chip")

    def test_unclassified_description_is_not_silently_accessibility(self) -> None:
        self.assertEqual(crosswalk._assay_family("unknown"), "unknown")


if __name__ == "__main__":
    unittest.main()
