from __future__ import annotations

import unittest

from scripts.build_gse256398_qc_sensitivity import (
    GSE256398QCSensitivityError,
    capped_retained_row_ids,
    pathological_auto_call,
)


class GSE256398QCSensitivityTests(unittest.TestCase):
    def test_pathological_call_requires_both_floors(self) -> None:
        self.assertTrue(pathological_auto_call(31, 100))
        self.assertFalse(pathological_auto_call(30, 100))
        self.assertFalse(pathological_auto_call(25, 100))
        with self.assertRaises(GSE256398QCSensitivityError):
            pathological_auto_call(2, 1)

    def test_capped_membership_is_score_ranked_and_row_stable(self) -> None:
        rows = [
            {
                "row_id": f"row-{index:02d}",
                "basic_qc_retained": "true",
                "scrublet_score": str(index),
            }
            for index in range(20)
        ]
        retained = capped_retained_row_ids(rows)
        self.assertEqual(len(retained), 18)
        self.assertNotIn("row-19", retained)
        self.assertNotIn("row-18", retained)
        self.assertIn("row-17", retained)


if __name__ == "__main__":
    unittest.main()
