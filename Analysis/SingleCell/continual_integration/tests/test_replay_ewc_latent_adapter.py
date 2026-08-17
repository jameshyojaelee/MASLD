import json
import unittest
from pathlib import Path

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.replay_ewc_latent_adapter import _load_policy


ROOT = Path(__file__).resolve().parents[1]


class TestReplayEWCLatentAdapter(unittest.TestCase):
    def test_policy_is_bound_to_original_paper_setting_and_v26(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy = _load_policy(
            config, ROOT / "reference" / "replay_ewc_latent_adapter_policy_v27.json"
        )
        self.assertEqual(policy["parent"]["required_method"], "continual_learning")
        self.assertEqual(policy["parent"]["required_ewc_lambda"], 100.0)
        self.assertEqual(policy["parent"]["required_replay_fraction"], 0.2)
        self.assertFalse(policy["diagnostic_basis"]["query_labels_may_select_v27"])

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            _load_policy(config, ROOT / "reference" / "reference_policy_v1.json")


if __name__ == "__main__":
    unittest.main()
