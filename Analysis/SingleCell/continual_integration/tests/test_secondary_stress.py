import unittest
from pathlib import Path

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.secondary_stress import _load_policy


ROOT = Path(__file__).resolve().parents[1]


class TestSecondaryStress(unittest.TestCase):
    def test_policy_preserves_donor_minimum_and_mapping_only_boundary(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy, _ = _load_policy(
            config, ROOT / "reference" / "secondary_stress_policy_v34.json"
        )
        self.assertEqual(policy["minimum_donors_per_group"], 3)
        self.assertEqual(policy["mapping_only_study"], "GSE189600")
        self.assertTrue(policy["firewall"]["donors_are_the_inferential_unit"])

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            _load_policy(config, ROOT / "reference" / "reference_policy_v1.json")


if __name__ == "__main__":
    unittest.main()
