from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from scripts.probe_scooby_epicardioids_checkpoint_decoder import (
    ScoobyCheckpointDecoderProbeError,
    validate_task_boundary,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/scooby_epicardioids_checkpoint_decoder_probe.json"
SCRIPT = ROOT / "scripts/probe_scooby_epicardioids_checkpoint_decoder.py"


class ScoobyCheckpointDecoderProbeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_probe_is_decoder_only_and_not_rna_conditioned_atac(self) -> None:
        validate_task_boundary(self.config)
        self.assertFalse(self.config["probe_policy"]["full_sequence_forward_allowed"])
        self.assertFalse(self.config["task_scope"]["rna_conditioned_atac_eligible"])

    def test_probe_rejects_native_context_claim(self) -> None:
        mutated = copy.deepcopy(self.config)
        mutated["probe_policy"]["native_cell_context_claim_allowed"] = True
        with self.assertRaises(ScoobyCheckpointDecoderProbeError):
            validate_task_boundary(mutated)

    def test_source_uses_only_safetensors_loader(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("from safetensors.torch import load_file", source)
        self.assertNotIn("torch.load(", source)
        self.assertNotIn("pickle", source)


if __name__ == "__main__":
    unittest.main()
