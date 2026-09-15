from __future__ import annotations

import unittest

from scripts.build_gse256398_donor_metadata import fibrosis_metadata


class GSE256398DonorMetadataTests(unittest.TestCase):
    def test_exact_and_interval_source_stages_remain_distinct(self) -> None:
        self.assertEqual(
            fibrosis_metadata("S19, Human, MASLD F0", "masld_f0")[
                "fibrosis_numeric_use"
            ],
            "exact",
        )
        interval = fibrosis_metadata(
            "S33, Human, MASH Fibrosis F2-3", "mash_fibrosis"
        )
        self.assertEqual((interval["fibrosis_min"], interval["fibrosis_max"]), (2, 3))
        self.assertEqual(interval["fibrosis_numeric_use"], "interval")

    def test_cirrhosis_is_not_silently_converted_to_f4(self) -> None:
        observed = fibrosis_metadata(
            "S15, Human, MASH Cirrhosis", "mash_cirrhosis"
        )
        self.assertEqual(observed["fibrosis_source_label"], "cirrhosis")
        self.assertEqual(observed["fibrosis_min"], "")
        self.assertEqual(observed["fibrosis_numeric_use"], "not_numeric_without_source_stage")

    def test_healthy_has_no_invented_histology_stage(self) -> None:
        observed = fibrosis_metadata(
            "S12, Human, Healthy Control", "healthy_control"
        )
        self.assertEqual(observed["fibrosis_state"], "structurally_missing")
        self.assertEqual(observed["fibrosis_source_label"], "")


if __name__ == "__main__":
    unittest.main()
