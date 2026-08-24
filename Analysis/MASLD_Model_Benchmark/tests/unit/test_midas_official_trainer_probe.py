from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "midas_official_trainer_probe.py"


class MIDASOfficialTrainerProbeSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = MODULE.read_text(encoding="utf-8")

    def test_official_train_path_is_called(self) -> None:
        self.assertIn("model.train(\n        max_epochs=1,", self.source)
        self.assertIn("limit_train_batches=2", self.source)

    def test_checkpoint_read_is_tensor_only(self) -> None:
        self.assertIn("weights_only=True", self.source)
        self.assertNotIn("load_from_checkpoint", self.source)


if __name__ == "__main__":
    unittest.main()
