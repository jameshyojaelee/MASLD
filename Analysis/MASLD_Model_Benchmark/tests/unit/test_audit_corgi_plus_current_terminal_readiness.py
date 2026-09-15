#!/usr/bin/env python3
"""Unit tests for the current Corgi+ terminal-readiness audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_corgi_plus_current_terminal_readiness.py"
SPEC = importlib.util.spec_from_file_location("corgi_plus_current", MODULE)
assert SPEC and SPEC.loader
current = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(current)


class CorgiPlusCurrentTerminalReadinessTest(unittest.TestCase):
    def test_placeholder_methods_must_raise_not_implemented(self) -> None:
        source = """
class corgiplus_pretrained:
    def predict(self):
        raise NotImplementedError()
    def predict_regions(self):
        raise NotImplementedError()
    def predict_regions_with_bigwig(self):
        raise NotImplementedError()
"""
        self.assertTrue(
            current.class_methods_raise_not_implemented(
                source,
                "corgiplus_pretrained",
                ("predict", "predict_regions", "predict_regions_with_bigwig"),
            )
        )

    def test_implemented_method_breaks_placeholder_gate(self) -> None:
        source = """
class corgiplus_pretrained:
    def predict(self):
        return 1
"""
        self.assertFalse(
            current.class_methods_raise_not_implemented(
                source, "corgiplus_pretrained", ("predict",)
            )
        )

    def test_zenodo_record_binds_weight_terms_and_identity(self) -> None:
        record = {
            "id": 18630048,
            "metadata": {"license": {"id": "cc-by-4.0"}},
            "files": [
                {
                    "key": "corgiplus_model.pt",
                    "size": current.CHECKPOINT_SIZE,
                    "checksum": f"md5:{current.CHECKPOINT_MD5}",
                }
            ],
        }
        result = current.zenodo_terms_and_files(record, 18630048)
        self.assertEqual(result["license"], "cc-by-4.0")

    def test_no_project_fixture_guesses_local_rna_or_baseline(self) -> None:
        decisions = current.fixture_decisions()
        self.assertTrue(decisions["gse244832"]["false_cell_pairing_forbidden"])
        self.assertEqual(decisions["gse281367"]["RNA"], "structurally_missing")
        self.assertFalse(
            any(
                row["corgi_plus_family_native_input_complete"]
                for row in decisions.values()
            )
        )

    def test_gse296875_gene_counts_are_not_local_total_rna(self) -> None:
        decision = current.fixture_decisions()["gse296875"]
        self.assertEqual(
            decision["local_stranded_total_RNA"],
            "not_native_to_3prime_single_nucleus_gene_expression",
        )
        self.assertFalse(decision["corgi_plus_family_native_input_complete"])


if __name__ == "__main__":
    unittest.main()
