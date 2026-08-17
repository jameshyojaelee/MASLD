import unittest
from pathlib import Path

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.production_stress import _load_policy


ROOT = Path(__file__).resolve().parents[1]


class TestProductionStress(unittest.TestCase):
    def test_policy_cannot_mislabel_development_as_confirmation(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy, _ = _load_policy(
            config, ROOT / "reference" / "production_stress_policy_v36.json"
        )
        self.assertTrue(policy["status"]["GSE136103_used_during_method_development"])
        self.assertTrue(policy["status"]["results_cannot_establish_independent_secondary_confirmation"])

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            _load_policy(config, ROOT / "reference" / "reference_policy_v1.json")


if __name__ == "__main__":
    unittest.main()
