from __future__ import annotations

import unittest

from scripts.diagnose_gse105127_rna_remote_identity import parse_content_range


class GSE105127RemoteIdentityTests(unittest.TestCase):
    def test_content_range_total_is_fail_closed(self) -> None:
        self.assertEqual(parse_content_range("bytes 0-0/1478606522"), 1_478_606_522)
        self.assertIsNone(parse_content_range("bytes 0-0/*"))
        self.assertIsNone(parse_content_range(None))


if __name__ == "__main__":
    unittest.main()
