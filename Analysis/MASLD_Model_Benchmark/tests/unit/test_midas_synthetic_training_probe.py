from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "midas_synthetic_training_probe.py"
SPEC = importlib.util.spec_from_file_location("midas_probe_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
probe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(probe)


class MIDASSyntheticProbeTests(unittest.TestCase):
    def test_small_config_disables_worker_persistence(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn('"num_workers": 0', source)
        self.assertIn('"persistent_workers": False', source)

    def test_atac_is_split_and_rna_only_translation_is_required(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn('dims_x={"rna": [20], "atac": [10, 20]}', source)
        self.assertIn('"x": {"rna": batch["x"]["rna"]}', source)
        self.assertIn('if "atac" not in translated:', source)

    def test_only_weights_only_state_loading_is_used(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertIn("weights_only=True", source)
        self.assertNotIn("load_from_checkpoint", source)
        self.assertNotIn("MIDAS.load(", source)


if __name__ == "__main__":
    unittest.main()
