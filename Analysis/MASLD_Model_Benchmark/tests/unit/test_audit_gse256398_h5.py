from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.audit_gse256398_h5 import (
    GSE256398AuditError,
    disease_group,
    parse_soft,
)


class GSE256398H5AuditTests(unittest.TestCase):
    def test_disease_mapping_keeps_alcohol_separate(self) -> None:
        self.assertEqual(
            disease_group("S1, Human, Alcohol-associated Cirrhosis"),
            "alcohol_associated_cirrhosis",
        )
        self.assertEqual(disease_group("S19, Human, MASLD F0"), "masld_f0")
        with self.assertRaises(GSE256398AuditError):
            disease_group("S1, Human, generic liver disease")

    def test_soft_parser_requires_full_source_census(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "small.soft.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("^SAMPLE = GSM1\n")
                handle.write("!Sample_title = S1, Human, Healthy Control\n")
                handle.write("!Sample_organism_ch1 = Homo sapiens\n")
            with self.assertRaisesRegex(GSE256398AuditError, "census"):
                parse_soft(path)


if __name__ == "__main__":
    unittest.main()
