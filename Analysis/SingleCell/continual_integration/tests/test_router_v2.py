import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from masld_cl.config import load_config
from masld_cl.router_v2 import (
    RouterV2Error,
    _validate_counts,
    collapse_router_labels,
    load_router_v2_policy,
)
from masld_cl.router_v2_evaluation import (
    RouterV2EvaluationError,
    paired_donor_bootstrap,
)


class TestRouterV2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.config = load_config(cls.root / "config_v1.json")
        cls.policy_path = cls.root / "reference" / "router_v2_policy_v43.json"

    def test_policy_is_locked_before_external_labels(self):
        _, policy, _, _ = load_router_v2_policy(self.config, self.policy_path)
        self.assertTrue(policy["selection_status"]["gse212837_is_diagnostic_only"])
        self.assertFalse(policy["selection_status"]["gse212837_may_change_method_or_thresholds"])
        self.assertFalse(policy["selection_status"]["gse296875_author_labels_opened_before_lock"])
        self.assertEqual(policy["validation"]["minimum_donor_balanced_macro_f1"], 0.70)

    def test_policy_copy_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            copied = Path(tmp) / "policy.json"
            copied.write_bytes(self.policy_path.read_bytes())
            with self.assertRaises(RouterV2Error):
                load_router_v2_policy(self.config, copied)

    def test_locked_label_collapse(self):
        _, policy, _, _ = load_router_v2_policy(self.config, self.policy_path)
        observed = collapse_router_labels(
            ["T cells", "Resident NK", "Mono+mono derived cells", "Hepatocytes"],
            policy["collapsed_labels"],
        )
        np.testing.assert_array_equal(observed, ["NK-T", "NK-T", "Macrophages", "Hepatocytes"])
        with self.assertRaises(RouterV2Error):
            collapse_router_labels(["new label"], policy["collapsed_labels"])

    def test_counts_fail_closed(self):
        _validate_counts(csr_matrix([[1, 0], [0, 2]], dtype=float))
        for values in ([[1.5, 0]], [[-1, 2]], [[0, 0]]):
            with self.assertRaises(RouterV2Error):
                _validate_counts(csr_matrix(values, dtype=float))

    def test_paired_donor_bootstrap_detects_improvement(self):
        observed = paired_donor_bootstrap(
            {"a": 0.8, "b": 0.9, "c": 0.85},
            {"a": 0.2, "b": 0.3, "c": 0.25}, 2000, 17,
        )
        self.assertGreater(observed["ci_low"], 0)
        with self.assertRaises(RouterV2EvaluationError):
            paired_donor_bootstrap(
                {"a": 0.8, "b": 0.9, "c": 0.85},
                {"a": 0.2, "b": 0.3}, 100, 17,
            )


if __name__ == "__main__":
    unittest.main()
