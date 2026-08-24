from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "epibert_probe_released_checkpoint.py"
SPEC = importlib.util.spec_from_file_location("epibert_checkpoint_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
probe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(probe)


class EpiBERTReleasedCheckpointProbeTests(unittest.TestCase):
    def test_pretraining_constructor_matches_released_checkpoint_audit(self) -> None:
        contract = probe.constructor_contract("pretrained_atac")
        self.assertEqual(contract["filter_list_seq"], [512, 640, 640, 768, 896, 1024])
        self.assertEqual(contract["num_heads"], 8)
        self.assertEqual(contract["seed"], 21)
        self.assertEqual(contract["final_output_length"], 4092)
        self.assertNotIn("predict_atac", contract)

    def test_rampage_constructor_matches_released_command(self) -> None:
        contract = probe.constructor_contract("fine_tuned_rampage")
        self.assertEqual(contract["num_heads"], 8)
        self.assertEqual(contract["seed"], 19)
        self.assertEqual(contract["final_output_length"], 896)
        self.assertTrue(contract["predict_atac"])
        self.assertEqual(contract["num_heads"], 8)

    def test_unknown_model_kind_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            probe.EpiBERTCheckpointProbeError, "model kind"
        ):
            probe.constructor_contract("sequence_only")

    def test_motif_model_input_matches_released_deserializer_batching(self) -> None:
        fixture = np.zeros((1, 693), dtype=np.float32)
        model_input = probe._motif_model_input(fixture)
        self.assertEqual(model_input.shape, (1, 1, 693))
        self.assertTrue(np.shares_memory(fixture, model_input))

    def test_motif_model_input_rejects_documented_but_unexecutable_rank(self) -> None:
        with self.assertRaisesRegex(
            probe.EpiBERTCheckpointProbeError, "motif fixture"
        ):
            probe._motif_model_input(np.zeros((1, 1, 693), dtype=np.float32))


if __name__ == "__main__":
    unittest.main()
