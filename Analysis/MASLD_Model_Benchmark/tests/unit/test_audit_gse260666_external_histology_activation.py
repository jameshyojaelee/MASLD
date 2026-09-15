from __future__ import annotations

import unittest

from scripts.audit_gse260666_external_histology_activation import lock_modal_recipe


class ExternalHistologyActivationTests(unittest.TestCase):
    def test_modal_primary_recipe_is_locked_without_secondary_endpoints(self) -> None:
        modal = {
            "feature_request": 1000,
            "pca_components": 5,
            "c": 0.001,
            "l1_ratio": 0.0,
            "neighbor_count": 0,
            "temperature": 1.0,
        }
        alternative = {**modal, "pca_components": 20, "c": 0.01}
        rows = [dict(modal) for _ in range(13)] + [
            dict(alternative) for _ in range(12)
        ]
        result = lock_modal_recipe(rows)
        self.assertEqual(result["parameters"], modal)
        self.assertEqual(result["selection_count"], 13)
        self.assertEqual(result["distinct_recipe_count"], 2)

    def test_modal_tie_prefers_lower_effective_complexity(self) -> None:
        simple = {
            "feature_request": 1000,
            "pca_components": 5,
            "c": 0.001,
            "l1_ratio": 0.0,
            "neighbor_count": 0,
            "temperature": 1.0,
        }
        complex_recipe = {**simple, "pca_components": 20, "c": 0.1}
        third = {**simple, "pca_components": 10, "c": 0.01}
        rows = [dict(simple) for _ in range(10)]
        rows += [dict(complex_recipe) for _ in range(10)]
        rows += [dict(third) for _ in range(5)]
        result = lock_modal_recipe(rows)
        self.assertEqual(result["parameters"], simple)
        self.assertEqual(result["tied_modal_recipe_count"], 2)


if __name__ == "__main__":
    unittest.main()
