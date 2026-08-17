import unittest
from pathlib import Path

from masld_cl.class_routed_adapter import _load_policy
from masld_cl.config import load_config
from masld_cl.contracts import ContractError


ROOT = Path(__file__).resolve().parents[1]


class TestClassRoutedAdapter(unittest.TestCase):
    def test_policy_forbids_query_labels(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy = _load_policy(
            config, ROOT / "reference" / "class_routed_replay_ewc_adapter_policy_v28.json"
        )
        self.assertFalse(
            policy["routing_classifier"]["query_labels_available_to_fit_or_route"]
        )
        self.assertTrue(
            policy["firewall"]["query_audit_labels_not_used_to_fit_route_construct_or_select"]
        )

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            _load_policy(config, ROOT / "reference" / "reference_policy_v1.json")

    def test_v29_compartments_cover_reference_vocabulary_once(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy = _load_policy(
            config,
            ROOT / "reference" / "compartment_routed_replay_ewc_adapter_policy_v29.json",
        )
        labels = [label for values in policy["routing_groups"].values() for label in values]
        self.assertEqual(len(labels), len(set(labels)))
        self.assertEqual(len(policy["routing_groups"]), 6)

    def test_v30_merges_all_immune_labels(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy = _load_policy(
            config,
            ROOT / "reference" / "immune_compartment_replay_ewc_adapter_policy_v30.json",
        )
        self.assertEqual(len(policy["routing_groups"]), 5)
        self.assertIn("T cells", policy["routing_groups"]["immune"])
        self.assertIn("Macrophages", policy["routing_groups"]["immune"])


if __name__ == "__main__":
    unittest.main()
