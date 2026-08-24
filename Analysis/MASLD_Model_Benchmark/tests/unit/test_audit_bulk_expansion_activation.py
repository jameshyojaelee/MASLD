from __future__ import annotations

import unittest

from scripts.audit_bulk_expansion_activation import (
    BulkActivationError,
    assign_folds,
    characteristics,
    relation_accessions,
    validate_gse274114_sample_id,
)


class BulkExpansionActivationTests(unittest.TestCase):
    def test_characteristics_are_case_normalized(self) -> None:
        row = {"characteristics_json": '["Sex: Male", "case/control: classic-NAFLD"]'}
        self.assertEqual(characteristics(row), {"sex": "Male", "case/control": "classic-NAFLD"})

    def test_folds_preserve_every_participant_once(self) -> None:
        rows = [
            {"participant_id": f"P{i}", "group": "A" if i < 6 else "B", "fibrosis": "F1"}
            for i in range(11)
        ]
        observed = assign_folds(rows)
        self.assertEqual(len(observed), 11)
        self.assertEqual(set(observed.values()), {0, 1, 2, 3, 4})

    def test_relation_accessions_are_canonical(self) -> None:
        observed = relation_accessions(
            '["BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN40215616", '
            '"SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX23809226"]'
        )
        self.assertEqual(observed, ("SAMN40215616", "SRX23809226"))

    def test_relation_accessions_fail_closed(self) -> None:
        with self.assertRaises(BulkActivationError):
            relation_accessions('["BioSample: https://example.org/SAMN40215616"]')

    def test_gse274114_source_sample_id_families(self) -> None:
        for value in ("FFPE16", "105925-001-001", "106039-001-010"):
            validate_gse274114_sample_id(value)
        with self.assertRaises(BulkActivationError):
            validate_gse274114_sample_id("subject-1")


if __name__ == "__main__":
    unittest.main()
