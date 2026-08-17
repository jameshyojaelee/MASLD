import unittest
from pathlib import Path

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.production_expansion import _load_policy


ROOT = Path(__file__).resolve().parents[1]


class TestProductionExpansion(unittest.TestCase):
    def test_policy_forbids_condition_and_descriptive_fit(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy, _ = _load_policy(
            config, ROOT / "reference" / "production_expansion_policy_v35.json"
        )
        self.assertIn("query_control", policy["adaptation"]["forbidden_inputs"])
        self.assertTrue(policy["firewall"]["descriptive_only_cells_are_unavailable_to_fit"])
        self.assertEqual(
            {spec["lineage_label"] for name, spec in policy["models"].items()
             if name != "all_lineage"},
            set(config["lineages"]),
        )

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            _load_policy(config, ROOT / "reference" / "reference_policy_v1.json")


if __name__ == "__main__":
    unittest.main()
