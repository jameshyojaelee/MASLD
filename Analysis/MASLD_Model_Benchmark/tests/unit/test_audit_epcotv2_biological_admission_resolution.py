#!/usr/bin/env python3
"""Unit tests for the EPCOTv2 biological-admission resolution audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_epcotv2_biological_admission_resolution.py"
SPEC = importlib.util.spec_from_file_location("epcotv2_resolution", MODULE)
assert SPEC and SPEC.loader
audit_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit_module)


class EPCOTv2BiologicalAdmissionResolutionTest(unittest.TestCase):
    def test_literal_DNA_alphabet_is_resolved_without_import(self) -> None:
        source = b"DNA = ['A', 'C', 'G', 'T']\n"
        self.assertEqual(
            audit_module.literal_assignment(source, "DNA"), ["A", "C", "G", "T"]
        )

    def test_one_hot_neutral_contract_is_source_parsed(self) -> None:
        source = b"""
def one_hot_dna(seq):
    return one_hot(seq, alphabet=DNA, neutral_alphabet=['N'], neutral_value=.25)
"""
        observed = audit_module.verify_one_hot_dna_contract(source)
        self.assertEqual(observed["neutral_alphabet"], ["N"])
        self.assertEqual(observed["neutral_value_per_channel"], 0.25)

    def test_one_hot_contract_fails_without_explicit_neutral_value(self) -> None:
        source = b"""
def one_hot_dna(seq):
    return one_hot(seq, alphabet=DNA, neutral_alphabet=['N'])
"""
        with self.assertRaises(audit_module.EPCOTv2ResolutionError):
            audit_module.verify_one_hot_dna_contract(source)

    def test_terms_do_not_transfer_from_adjacent_code(self) -> None:
        metadata = {
            "sha": "abc",
            "cardData": {"title": "model"},
            "siblings": [{"rfilename": "README.md"}],
        }
        observed = audit_module.api_terms(metadata, "abc")
        self.assertFalse(observed["terms_declared"])
        self.assertIsNone(observed["declared_license"])

    def test_sdist_reader_never_imports_package(self) -> None:
        import io
        import tarfile

        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "source.tar.gz"
            content = b"DNA = ['A', 'C', 'G', 'T']\n"
            info = tarfile.TarInfo("kipoiseq-0.5.2/kipoiseq/utils.py")
            info.size = len(content)
            with tarfile.open(archive, "w:gz") as handle:
                handle.addfile(info, io.BytesIO(content))
            self.assertEqual(
                audit_module.read_unique_tar_suffix(archive, "/kipoiseq/utils.py"),
                content,
            )


if __name__ == "__main__":
    unittest.main()
