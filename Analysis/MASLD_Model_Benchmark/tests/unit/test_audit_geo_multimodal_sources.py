from __future__ import annotations

import gzip
import unittest

from scripts.audit_geo_multimodal_sources import (
    SourceAuditError,
    audit_znf469,
    audit_zonated,
    parse_soft,
)


def _payload(titles: list[str]) -> bytes:
    lines = []
    for index, title in enumerate(titles):
        lines.extend(
            [
                f"^SAMPLE = GSM{index + 1}",
                f"!Sample_title = {title}",
                "!Sample_source_name_ch1 = Liver",
                "!Sample_characteristics_ch1 = tissue: Liver",
            ]
        )
    return gzip.compress(("\n".join(lines) + "\n").encode())


class GeoMultimodalSourceAuditTests(unittest.TestCase):
    def test_parse_soft_preserves_repeated_fields(self) -> None:
        records = parse_soft(_payload(["D1_RNA", "D2_RNA"]))
        self.assertEqual([record["accession"] for record in records], ["GSM1", "GSM2"])
        self.assertEqual(records[0]["characteristics"], ["tissue: Liver"])

    def test_znf_join_is_diagnostic_not_authoritative(self) -> None:
        rna = parse_soft(_payload(["D1", "D2", "D2_rep2"]))
        h3 = parse_soft(_payload(["D1_nor_H3K27ac", "D3_nash_H3K27ac"]))
        audit = audit_znf469(rna, h3)
        self.assertEqual(audit["shared_candidate_ids"], 1)
        self.assertEqual(audit["rna_duplicate_candidate_ids"], {"D2": 2})
        self.assertFalse(audit["candidate_join_is_authoritative"])
        self.assertFalse(audit["admission_ready"])

    def test_zonated_requires_complete_donor_cross_product(self) -> None:
        titles = [
            f"{donor}_{zone}_{assay}"
            for donor in range(1, 20)
            for zone in ("CV", "IZ", "PP")
            for assay in ("RNA", "RRBS")
        ]
        audit = audit_zonated(parse_soft(_payload(titles)))
        self.assertEqual(audit["sample_records"], 114)
        self.assertEqual(audit["complete_candidate_donors"], 19)
        with self.assertRaises(SourceAuditError):
            audit_zonated(parse_soft(_payload(titles[:-1])))


if __name__ == "__main__":
    unittest.main()
