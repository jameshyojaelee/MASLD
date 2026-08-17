import unittest

import numpy as np
import pandas as pd
import torch

from masld_cl.external_reference_common import (
    ExternalReferenceError, _enable_inference_without_optimization,
    common_gene_mapping, external26_roster_mask,
)


class TestExternalReferenceCommon(unittest.TestCase):
    def test_common_gene_mapping_prefers_unique_id_and_excludes_collisions(self):
        result = common_gene_mapping(
            ["A", "B", "C", "D"],
            ["A", "x", "y", "z"],
            ["other", "B", "C", "C"],
        )
        self.assertTrue(np.array_equal(result["canonical_indices"], np.asarray([0, 1])))
        self.assertTrue(np.array_equal(result["external_indices"], np.asarray([0, 1])))
        self.assertEqual(result["mapping_methods"], {
            "external_feature_name": 1,
            "external_var_id": 1,
        })

    def test_common_gene_mapping_rejects_duplicate_canonical_ids(self):
        with self.assertRaisesRegex(ExternalReferenceError, "duplicated"):
            common_gene_mapping(["A", "A"], ["A"], ["A"])

    def test_external26_rosters_are_nested(self):
        obs = pd.DataFrame({
            "strict_reference": [True, False, False],
            "external_reference_candidate": [False, True, False],
        })
        self.assertTrue(np.array_equal(
            external26_roster_mask(obs, "common_strict7"), [True, False, False]
        ))
        self.assertTrue(np.array_equal(
            external26_roster_mask(obs, "external_clean26"), [True, True, False]
        ))
        with self.assertRaisesRegex(ExternalReferenceError, "unsupported"):
            external26_roster_mask(obs, "bad")

    def test_inference_flag_does_not_change_loaded_mapper_state(self):
        class Mapper:
            def __init__(self):
                self.module = torch.nn.Linear(2, 2)
                self.is_trained = False

        mapper = Mapper()
        before = {
            name: value.detach().clone()
            for name, value in mapper.module.state_dict().items()
        }
        result = _enable_inference_without_optimization(mapper)
        self.assertTrue(mapper.is_trained)
        self.assertFalse(result["load_query_data_weights_optimized"])
        for name, value in mapper.module.state_dict().items():
            self.assertTrue(torch.equal(value, before[name]))


if __name__ == "__main__":
    unittest.main()
