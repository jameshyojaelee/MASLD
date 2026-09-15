#!/usr/bin/env python3
"""Unit tests for EPCOTv2's source-defined native contract audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_epcotv2_native_contract.py"
SPEC = importlib.util.spec_from_file_location("epcotv2_native_contract", MODULE)
assert SPEC and SPEC.loader
audit_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit_module)


class EPCOTv2NativeContractTest(unittest.TestCase):
    def test_neighbor_padding_uses_adjacent_bins(self) -> None:
        sequence = np.zeros((3, 4, 1000), dtype=np.int8)
        sequence[0, 0] = 1
        sequence[1, 1] = 1
        sequence[2, 2] = 1
        padded = audit_module.pad_sequence(sequence)
        self.assertEqual(padded.shape, (3, 4, 1600))
        np.testing.assert_array_equal(padded[1, :, :300], sequence[0, :, -300:])
        np.testing.assert_array_equal(padded[1, :, 300:1300], sequence[1])
        np.testing.assert_array_equal(padded[1, :, 1300:], sequence[2, :, :300])

    def test_raw_atac_locus_propagates_to_three_model_bins(self) -> None:
        raw_mask = np.zeros((5, 1000), dtype=np.float32)
        raw_mask[2] = 1
        propagated = audit_module.pad_signal(raw_mask) > 0
        affected = np.flatnonzero(propagated.any(axis=1)).tolist()
        self.assertEqual(affected, [1, 2, 3])
        self.assertEqual(int(propagated.sum()), 1600)

    def test_terms_require_card_or_license_file(self) -> None:
        absent = {"sha": "abc", "cardData": {}, "siblings": [{"rfilename": "README.md"}]}
        self.assertFalse(audit_module.api_terms(absent, "abc")["terms_declared"])
        card = {"sha": "abc", "cardData": {"license": "mit"}, "siblings": []}
        self.assertTrue(audit_module.api_terms(card, "abc")["terms_declared"])
        file_record = {
            "sha": "abc",
            "cardData": {},
            "siblings": [{"rfilename": "LICENSE.txt"}],
        }
        self.assertTrue(audit_module.api_terms(file_record, "abc")["terms_declared"])

    def test_source_defined_fixture_is_biological_free(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            receipt = audit_module.build_source_defined_fixture(output)
            self.assertFalse(receipt["base_to_channel_mapping_claimed"])
            self.assertEqual(receipt["atac_mask_propagated_model_bins"], [299, 300, 301])
            with np.load(output / receipt["path"], allow_pickle=False) as fixture:
                self.assertEqual(
                    fixture["sequence_plus_observed_atac"].shape,
                    (1, 600, 5, 1600),
                )
                self.assertEqual(
                    int(fixture["propagated_raw_atac_locus_mask"].sum()), 1600
                )

    def test_wrong_shapes_fail_closed(self) -> None:
        with self.assertRaises(audit_module.EPCOTv2AuditError):
            audit_module.pad_sequence(np.zeros((2, 4, 999), dtype=np.int8))
        with self.assertRaises(audit_module.EPCOTv2AuditError):
            audit_module.pad_signal(np.zeros((2, 999), dtype=np.float32))


if __name__ == "__main__":
    unittest.main()
