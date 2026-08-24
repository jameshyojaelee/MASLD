from __future__ import annotations

import gzip
import unittest

from scripts.audit_bulk_expansion_sources import parse_soft


class BulkExpansionSourceAuditTests(unittest.TestCase):
    def test_soft_parser_keeps_series_and_sample_evidence_separate(self) -> None:
        payload = gzip.compress(
            b"!Series_title = cohort\n"
            b"^SAMPLE = GSM1\n"
            b"!Sample_title = person1\n"
            b"!Sample_source_name_ch1 = liver\n"
            b"!Sample_organism_ch1 = Homo sapiens\n"
            b"!Sample_characteristics_ch1 = group: MASH\n"
            b"!Sample_relation = BioSample: SAMN1\n"
        )
        series, samples = parse_soft(payload)
        self.assertEqual(series["!Series_title"], ["cohort"])
        self.assertEqual(samples[0]["accession"], "GSM1")
        self.assertEqual(samples[0]["characteristics"], ["group: MASH"])


if __name__ == "__main__":
    unittest.main()
