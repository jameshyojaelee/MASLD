import unittest
from pathlib import Path

import numpy as np

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.gpu_tolerance import _arguments_without_output, _load_policy, _numeric_tolerance


ROOT = Path(__file__).resolve().parents[1]


class TestGPUTolerance(unittest.TestCase):
    def test_policy_is_report_only_and_exact_setting(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy, _ = _load_policy(
            config, ROOT / "reference" / "gpu_tolerance_policy_v41.json",
            verify_hashes=False,
        )
        self.assertTrue(policy["report_only_no_selection_threshold"])
        self.assertEqual(policy["required_setting"]["ewc_lambda"], 100.0)
        self.assertEqual(policy["required_setting"]["replay_fraction"], 0.2)

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            _load_policy(
                config, ROOT / "reference" / "reference_policy_v1.json",
                verify_hashes=False,
            )

    def test_numeric_tolerance_reports_exact_and_perturbed(self):
        left = np.arange(12, dtype=np.float32).reshape(3, 4)
        exact = _numeric_tolerance(left, left.copy())
        self.assertTrue(exact["bitwise_identical"])
        self.assertEqual(exact["maximum_absolute_difference"], 0.0)
        right = left.copy()
        right[1, 2] += 1e-4
        changed = _numeric_tolerance(left, right)
        self.assertFalse(changed["bitwise_identical"])
        self.assertEqual(changed["nonzero_elements"], 1)

    def test_output_argument_is_the_only_allowed_difference(self):
        left = ["--seed", "17", "--output", "/a"]
        right = ["--seed", "17", "--output", "/b"]
        self.assertEqual(_arguments_without_output(left), _arguments_without_output(right))


if __name__ == "__main__":
    unittest.main()
