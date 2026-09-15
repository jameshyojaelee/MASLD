#!/usr/bin/env python3
"""Unit checks for the EpiBERT native-preprocessing adjudicator."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_epibert_native_preprocessing.py"
SPEC = importlib.util.spec_from_file_location("epibert_native", MODULE)
assert SPEC and SPEC.loader
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class EpiBERTNativePreprocessingTest(unittest.TestCase):
    def test_motif_roster_preserves_declaration_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "motifs.meme"
            path.write_text(
                "MEME version 4\nMOTIF AC0001:A first\nMOTIF AC0002:B second\n",
                encoding="utf-8",
            )
            self.assertEqual(
                audit.motif_roster(path),
                [("AC0001:A", "first"), ("AC0002:B", "second")],
            )

    def test_comma_values_requires_numeric_entries(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "values.tsv"
            path.write_text("1.0,2.0,3.0", encoding="utf-8")
            self.assertEqual(audit.comma_values(path), [1.0, 2.0, 3.0])
            path.write_text("1.0,not-a-number", encoding="utf-8")
            with self.assertRaises(ValueError):
                audit.comma_values(path)


if __name__ == "__main__":
    unittest.main()
