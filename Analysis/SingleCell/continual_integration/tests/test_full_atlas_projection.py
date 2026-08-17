import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from masld_cl.config import load_config
from masld_cl.contracts import ContractError
from masld_cl.full_atlas_projection import (
    _project_frozen_bridge,
    load_full_atlas_projection_policy,
)


ROOT = Path(__file__).resolve().parents[1]


class TestFullAtlasProjection(unittest.TestCase):
    def test_policy_is_bound_to_passing_v32_and_exact_census(self):
        config = load_config(ROOT / "config_v1.json")
        _, policy, _ = load_full_atlas_projection_policy(
            config, ROOT / "reference" / "full_atlas_projection_policy_v33.json"
        )
        self.assertEqual(policy["expected_contract"]["descriptive_cells"], 1_232_318)
        self.assertEqual(policy["expected_contract"]["descriptive_only_cells"], 33)
        self.assertFalse(policy["projection"]["secondary_query_labels_available_to_fit"])

    def test_policy_copy_fails_closed(self):
        config = load_config(ROOT / "config_v1.json")
        with self.assertRaises(ContractError):
            load_full_atlas_projection_policy(
                config, ROOT / "reference" / "reference_policy_v1.json"
            )

    def test_frozen_bridge_matches_explicit_formula(self):
        raw = np.asarray([[1.0, 2.0], [3.0, 5.0]], dtype=np.float32)
        bridge = {"bridge": {
            "source_mean": [1.0, 1.0], "target_mean": [0.5, -0.5],
            "rotation": [[0.0, 1.0], [1.0, 0.0]], "scale": 0.5,
        }}
        with tempfile.TemporaryDirectory() as tmp:
            output = np.lib.format.open_memmap(
                Path(tmp) / "x.npy", mode="w+", dtype=np.float32, shape=raw.shape
            )
            _project_frozen_bridge(raw, bridge, output, chunk_size=1)
            expected = ((raw.astype(np.float64) - [1.0, 1.0])
                        @ np.asarray([[0.0, 1.0], [1.0, 0.0]])
                        * 0.5 + [0.5, -0.5]).astype(np.float32)
            np.testing.assert_array_equal(output, expected)


if __name__ == "__main__":
    unittest.main()
