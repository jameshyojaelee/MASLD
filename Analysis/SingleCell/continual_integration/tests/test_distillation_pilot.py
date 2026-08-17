from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.distillation_pilot import load_distillation_policy


class TestDistillationPilot(unittest.TestCase):
    def setUp(self):
        self.root = Path(__file__).resolve().parents[1]
        self.config = load_config(self.root / "config_v1.json")
        self.policy = self.root / "reference" / "latent_distillation_policy_v3.json"

    def test_source_controlled_policy_is_bound_to_failed_v2(self):
        value = load_distillation_policy(self.config, self.policy)
        self.assertEqual(value["alpha_grid"], [100.0, 300.0, 1000.0, 3000.0, 10000.0])
        self.assertTrue(value["case_stage_program_hero_gene_and_cas13_outcomes_locked"])

    def test_policy_copy_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / self.policy.name
            copied.write_text(self.policy.read_text())
            with self.assertRaisesRegex(ContractError, "source-controlled"):
                load_distillation_policy(self.config, copied)


if __name__ == "__main__":
    unittest.main()
