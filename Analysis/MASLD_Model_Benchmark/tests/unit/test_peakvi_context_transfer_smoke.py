from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).parents[2]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
MODULE = SCRIPTS / "peakvi_context_transfer_smoke.py"
SPEC = importlib.util.spec_from_file_location("peakvi_context_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
context = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(context)


class PeakVIContextTransferTests(unittest.TestCase):
    def test_model_identity_is_not_rna_conditioned(self) -> None:
        self.assertEqual(context.MODEL_ID, "peakvi_training_context")
        self.assertNotIn("rna", context.MODEL_ID)

    def test_lineage_roster_is_exact(self) -> None:
        self.assertEqual(
            context.LINEAGES,
            ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"),
        )


if __name__ == "__main__":
    unittest.main()
