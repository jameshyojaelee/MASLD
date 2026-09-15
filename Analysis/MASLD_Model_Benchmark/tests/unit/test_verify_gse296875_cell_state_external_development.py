from __future__ import annotations

import unittest

from scripts.verify_gse296875_cell_state_external_development import (
    CLASSES,
    SOURCE_TO_BROAD,
    VerificationError,
    identifier,
    independent_select,
)


class IndependentExternalCellStateVerificationTests(unittest.TestCase):
    def sources(self):
        mapping = []
        labels = []
        for donor in ("a", "b", "c"):
            for source_label in SOURCE_TO_BROAD:
                for replicate in range(3):
                    cell = f"{donor}:{source_label}:{replicate}"
                    mapping.append({
                        "cell_id": cell,
                        "donor_id": donor,
                        "well_id": "well1",
                        "raw_barcode": f"well1_ACGT{replicate}-1",
                    })
                    labels.append({"cell_id": cell, "author_label": source_label})
        return mapping, labels

    def test_independent_selection_is_balanced_and_complete(self):
        mapping, labels = self.sources()
        selected = independent_select(mapping, labels, per_class=6, seed=31)
        self.assertEqual(len(selected), 30)
        for broad in CLASSES:
            rows = [row for row in selected if row["broad_label"] == broad]
            self.assertEqual(len(rows), 6)
            self.assertEqual({row["donor_id"] for row in rows}, {"a", "b", "c"})

    def test_identity_namespaces_are_distinct(self):
        self.assertNotEqual(identifier("row", "1"), identifier("donor", "1"))

    def test_join_mismatch_fails(self):
        mapping, labels = self.sources()
        with self.assertRaisesRegex(VerificationError, "identities differ"):
            independent_select(mapping, labels[:-1], per_class=3)


if __name__ == "__main__":
    unittest.main()
