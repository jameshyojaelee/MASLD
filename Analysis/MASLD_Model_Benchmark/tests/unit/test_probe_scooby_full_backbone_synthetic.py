from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from scripts.probe_scooby_full_backbone_synthetic import (
    ScoobyFullBackboneProbeError,
    validate_boundary,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/scooby_epicardioids_full_backbone_synthetic_probe.json"
SCRIPT = ROOT / "scripts/probe_scooby_full_backbone_synthetic.py"


class ScoobyFullBackboneSyntheticTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_probe_is_synthetic_and_not_biological(self) -> None:
        validate_boundary(self.config)
        self.assertTrue(self.config["task_scope"]["synthetic_architecture_admission_only"])
        self.assertFalse(self.config["task_scope"]["biological_task_executed"])
        self.assertFalse(self.config["task_scope"]["native_cell_context_used"])

    def test_probe_rejects_rna_conditioned_claim(self) -> None:
        mutated = copy.deepcopy(self.config)
        mutated["task_scope"]["rna_conditioned_atac_eligible"] = True
        with self.assertRaises(ScoobyFullBackboneProbeError):
            validate_boundary(mutated)

    def test_source_uses_safetensors_and_no_unrestricted_loader(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("from safetensors.torch import load_file", source)
        self.assertNotIn("torch.load(", source)
        self.assertNotIn("pickle", source)


if __name__ == "__main__":
    unittest.main()
