#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


validator = load_module("12_validate_protein_atac_adapters.py", "protein_atac_validator_test")
builder = load_module("11_build_protein_atac_adapters.py", "protein_atac_builder_test")


class ProteinAtacFailClosedTests(unittest.TestCase):
    def test_v1_source_and_registry_legacy_ids_are_explicitly_distinct(self):
        row = {"cell_type": "hepatocytes", "module": "20"}
        self.assertEqual(builder._v1_program_id(row), "hepatocytes::20")
        self.assertEqual(builder._legacy_id(row), "hepatocytes:20")

    def test_exact_permutation_effect_pvalue_and_welch_se(self):
        scores = {
            "c1": 0.0,
            "c2": 1.0,
            "c3": 2.0,
            "m1": 3.0,
            "m2": 4.0,
            "m3": 5.0,
        }
        effect, pvalue, std_error = validator.exact_two_group(scores, {"m1", "m2", "m3"})
        self.assertAlmostEqual(effect, 3.0)
        self.assertAlmostEqual(pvalue, 2 / 20)
        self.assertAlmostEqual(std_error, math.sqrt(2 / 3))

    def test_bh_is_complete_two_member_family(self):
        observed = validator.bh_adjust({"program_a": 0.01, "program_b": 0.04})
        self.assertAlmostEqual(observed["program_a"], 0.02)
        self.assertAlmostEqual(observed["program_b"], 0.04)

    def test_family_gate_rejects_missing_duplicate_and_duplicate_expectation(self):
        validator.assert_complete_family((item for item in ("a", "b")), ("a", "b"), "valid")
        with self.assertRaises(Exception):
            validator.assert_complete_family(("a",), ("a", "b"), "missing")
        with self.assertRaises(Exception):
            validator.assert_complete_family(("a", "a", "b"), ("a", "b"), "duplicate")
        with self.assertRaises(Exception):
            validator.assert_complete_family(("a", "b"), ("a", "a", "b"), "bad expectation")

    def test_byte_identity_gate_fails_on_same_length_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            expected = root / "expected.tsv"
            observed = root / "observed.tsv"
            expected.write_bytes(b"gene\tvalue\nA\t1\n")
            observed.write_bytes(expected.read_bytes())
            validator.assert_byte_identical(observed, expected, "fixture")
            observed.write_bytes(b"gene\tvalue\nB\t1\n")
            self.assertEqual(observed.stat().st_size, expected.stat().st_size)
            with self.assertRaises(Exception):
                validator.assert_byte_identical(observed, expected, "fixture")


if __name__ == "__main__":
    unittest.main()
