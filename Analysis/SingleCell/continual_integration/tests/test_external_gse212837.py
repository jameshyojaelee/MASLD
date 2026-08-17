import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.external_gse212837 import (
    _balanced_positions,
    _normalize_log_common,
    load_external_gse212837_policy,
)


ROOT = Path(__file__).resolve().parents[1]


class TestExternalGSE212837(unittest.TestCase):
    def test_policy_is_independent_and_condition_blind(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy, _ = load_external_gse212837_policy(
            config, ROOT / "reference" / "external_gse212837_policy_v40.json",
            verify_source_hashes=False,
        )
        self.assertFalse(policy["status"]["cohort_used_during_method_development"])
        self.assertIn("celltype_pred", policy["mapping"]["forbidden_mapping_inputs"])
        self.assertEqual(
            policy["evaluation_only"]["gates"]["required_evaluable_models"],
            ["all_lineage", "hepatocytes", "fibroblasts", "cholangiocytes"],
        )
        self.assertFalse(policy["models"]["macrophages"]["independently_evaluable"])
        self.assertFalse(policy["models"]["t_cells"]["independently_evaluable"])

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            load_external_gse212837_policy(
                config, ROOT / "reference" / "reference_policy_v1.json",
                verify_source_hashes=False,
            )

    def test_common_normalization_is_finite_and_fixed_sum_before_log(self):
        counts = sparse.csr_matrix(np.asarray([[1, 1, 0], [0, 3, 1]], dtype=np.int32))
        observed = _normalize_log_common(counts)
        expected = np.log1p(np.asarray([[5000, 5000, 0], [0, 7500, 2500]], dtype=float))
        np.testing.assert_allclose(observed.toarray(), expected, rtol=1e-6)

    def test_balanced_sampler_is_exact_and_deterministic(self):
        keys = ["b"] * 10 + ["a"] * 7
        left = _balanced_positions(keys, cap=3, seed=17)
        right = _balanced_positions(keys, cap=3, seed=17)
        np.testing.assert_array_equal(left, right)
        self.assertEqual(sum(np.asarray(keys)[left] == "a"), 3)
        self.assertEqual(sum(np.asarray(keys)[left] == "b"), 3)

    def test_uniform_translation_preserves_centered_geometry(self):
        rng = np.random.default_rng(17)
        before = rng.normal(size=(50, 30)).astype(np.float32)
        offset = rng.normal(size=30)
        after = (before.astype(np.float64) + offset).astype(np.float32)
        np.testing.assert_allclose(
            after - after.mean(axis=0), before - before.mean(axis=0), atol=2e-6
        )


if __name__ == "__main__":
    unittest.main()
