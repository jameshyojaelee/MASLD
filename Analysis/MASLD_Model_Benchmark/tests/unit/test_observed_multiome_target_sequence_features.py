from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.materialize_gse296875_observed_multiome_target_sequence_features import (
    IndexedFasta,
    TargetSequenceFeatureError,
    feature_names,
    reverse_complement,
    sequence_features,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_target_sequence_features_20260825.json"


class ObservedMultiomeTargetSequenceFeatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_config_passes_without_outcome_access(self) -> None:
        resolved = validate_config(ROOT, deepcopy(self.config))
        self.assertTrue(resolved["fasta"].is_file())
        self.assertFalse(self.config["feature_contract"]["fit_on_outcomes"])

    def test_feature_axis_and_reverse_complement_invariance(self) -> None:
        names = feature_names()
        self.assertEqual(len(names), 85)
        sequence = "ACGTN" * 205
        sequence = sequence[:1024]
        forward = sequence_features(sequence, 231)
        reverse = sequence_features(reverse_complement(sequence), 231)
        np.testing.assert_allclose(forward, reverse, rtol=0, atol=1e-7)

    def test_indexed_fasta_fetches_lines_and_boundary_padding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "tiny.fa"
            fai = root / "tiny.fa.fai"
            fasta.write_bytes(b">chr1\nACGT\nTGCA\n")
            fai.write_text("chr1\t8\t6\t4\t5\n", encoding="utf-8")
            with IndexedFasta(fasta, fai) as indexed:
                self.assertEqual(indexed.fetch("chr1", 2, 7), "GTTGC")
                self.assertEqual(indexed.fetch("chr1", -2, 3), "NNACG")
                self.assertEqual(indexed.fetch("chr1", 6, 10), "CANN")

    def test_outcome_fitted_feature_contract_is_rejected(self) -> None:
        mutated = deepcopy(self.config)
        mutated["feature_contract"]["fit_on_outcomes"] = True
        with self.assertRaises(TargetSequenceFeatureError):
            validate_config(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
