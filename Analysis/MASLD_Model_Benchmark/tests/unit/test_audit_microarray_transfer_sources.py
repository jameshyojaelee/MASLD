from __future__ import annotations

import gzip
import unittest

from scripts.audit_microarray_transfer_sources import (
    characteristic_key_counts,
    parse_soft,
)


class MicroarrayTransferSourceAuditTests(unittest.TestCase):
    def test_parser_retains_platform_and_repeated_characteristics(self) -> None:
        payload = gzip.compress(
            b"!Series_title = cohort\n"
            b"^SAMPLE = GSM1\n"
            b"!Sample_title = person1 baseline\n"
            b"!Sample_source_name_ch1 = liver\n"
            b"!Sample_organism_ch1 = Homo sapiens\n"
            b"!Sample_platform_id = GPL570\n"
            b"!Sample_characteristics_ch1 = fibrosis: 3\n"
            b"!Sample_characteristics_ch1 = nash: yes\n"
            b"!Sample_relation = BioSample: SAMN1\n"
        )
        series, samples = parse_soft(payload)
        self.assertEqual(series["!Series_title"], ["cohort"])
        self.assertEqual(samples[0]["platform_id"], "GPL570")
        self.assertEqual(
            characteristic_key_counts(samples), {"fibrosis": 1, "nash": 1}
        )

    def test_parser_rejects_missing_platform(self) -> None:
        payload = gzip.compress(
            b"^SAMPLE = GSM1\n"
            b"!Sample_title = person1\n"
            b"!Sample_organism_ch1 = Homo sapiens\n"
        )
        with self.assertRaisesRegex(RuntimeError, "platform"):
            parse_soft(payload)


if __name__ == "__main__":
    unittest.main()
