#!/usr/bin/env python3
"""Fail-closed tests for the GSE105127 exact-reference status repair."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import unittest

from masld_bench.hashing import sha256_file
from scripts import build_gse105127_cpg_crosswalk_status_repair as crosswalk
from scripts import finalize_gse105127_production_activation_status_repair as finalizer
from scripts.gse105127_exact_reference_contract import (
    REFERENCE_RECEIPT_STATUS,
    admits_exact_reference,
)


ROOT = Path(__file__).resolve().parents[2]
LOCKED_SOURCE_SHA256 = {
    "scripts/build_gse105127_cpg_crosswalk.py": "2676dbeb4db530a90f0816d34df480649d66b50755b0978922d95f418d37736e",
    "scripts/finalize_gse105127_production_activation.py": "fb891b529c545972c14f120c341256c993e3d54cbf13f8032499036d7b9d9ac3",
}


def exact_authorities() -> tuple[dict[str, object], dict[str, object]]:
    manifest = {
        "metadata": {
            "artifact_class": "gse105127_reference_bundle",
            "status": "passed_exact_legacy_stream",
        }
    }
    receipt = {
        "schema_version": "masld-bench-gse105127-reference-bundle-v2",
        "status": "passed_exact_legacy_stream",
        "source_assembly": "1000_Genomes_GRCh37_human_g1k_v37",
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    return manifest, receipt


class GSE105127ReferenceStatusRepairTests(unittest.TestCase):
    def test_exact_frozen_reference_contract_is_admitted(self) -> None:
        manifest, receipt = exact_authorities()
        self.assertTrue(admits_exact_reference(manifest, receipt))
        self.assertEqual(REFERENCE_RECEIPT_STATUS, "passed_exact_legacy_stream")

    def test_generic_passed_status_is_not_admitted(self) -> None:
        manifest, receipt = exact_authorities()
        manifest["metadata"]["status"] = "passed"
        receipt["status"] = "passed"
        self.assertFalse(admits_exact_reference(manifest, receipt))

    def test_manifest_receipt_or_firewall_mismatch_fails_closed(self) -> None:
        manifest, receipt = exact_authorities()
        mutations = (
            ("manifest_class", "metadata", "artifact_class", "other"),
            ("manifest_status", "metadata", "status", "other"),
            ("receipt_schema", None, "schema_version", "other"),
            ("source_assembly", None, "source_assembly", "other"),
            ("labels", None, "labels_accessed", True),
            ("fit", None, "fit_or_score_performed", True),
        )
        for name, section, key, value in mutations:
            with self.subTest(name=name):
                current_manifest = deepcopy(manifest)
                current_receipt = deepcopy(receipt)
                target = (
                    current_manifest[section]
                    if section is not None
                    else current_receipt
                )
                target[key] = value
                self.assertFalse(
                    admits_exact_reference(current_manifest, current_receipt)
                )

    def test_source_locked_v2_implementations_remain_byte_identical(self) -> None:
        for relative, expected in LOCKED_SOURCE_SHA256.items():
            with self.subTest(relative=relative):
                self.assertEqual(sha256_file(ROOT / relative), expected)

    def test_repaired_adapters_use_exact_contract_without_outcome_surface(self) -> None:
        for module in (crosswalk, finalizer):
            source = Path(module.__file__).read_text(encoding="utf-8")
            with self.subTest(module=module.__name__):
                self.assertIn("admits_exact_reference", source)
                self.assertIn("REFERENCE_RECEIPT_STATUS", source)
                self.assertNotIn("evaluator_only", source.lower())
                self.assertNotIn("outcomes.tsv", source.lower())
                self.assertNotIn('reference.get("status") != "passed"', source)


if __name__ == "__main__":
    unittest.main()
