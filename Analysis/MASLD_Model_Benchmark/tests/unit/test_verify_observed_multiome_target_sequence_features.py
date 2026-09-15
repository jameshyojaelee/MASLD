from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

import numpy as np

from scripts.verify_gse296875_observed_multiome_target_sequence_features import (
    TargetSequenceVerificationError,
    _features,
    _reverse_complement,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_target_sequence_feature_verification_20260825.json"


class TargetSequenceFeatureVerificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_parent_contracts_pass(self) -> None:
        resolved = validate_config(ROOT, deepcopy(self.config))
        self.assertTrue(resolved["fasta"].is_file())

    def test_independent_features_are_reverse_complement_invariant(self) -> None:
        sequence = ("ACGTN" * 205)[:1024]
        np.testing.assert_allclose(_features(sequence, 217), _features(_reverse_complement(sequence), 217), rtol=0, atol=1e-7)

    def test_open_firewall_is_rejected_before_parent_reads(self) -> None:
        mutated = deepcopy(self.config)
        mutated["firewall"]["model_fit"] = True
        with self.assertRaises(TargetSequenceVerificationError):
            validate_config(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()
