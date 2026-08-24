from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

try:
    import numpy as np
except ModuleNotFoundError:  # Control-plane Python intentionally has no NumPy.
    np = None  # type: ignore[assignment]

try:
    import torch
except ModuleNotFoundError:  # Control-plane Python intentionally has no PyTorch.
    torch = None  # type: ignore[assignment]


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ADAPTER_PATH = (
    PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_scpair.py"
)
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_scpair_unit", ADAPTER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load standalone scientific adapter: {ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)


class RNAATACScPairTests(unittest.TestCase):
    def _rows(self, donors: int = 20) -> list[dict[str, str]]:
        return [
            {
                "cell_id": f"cell-{donor:02d}-{lineage}",
                "donor_id": f"donor-{donor:02d}",
                "lineage": lineage,
                "well_id": f"well{donor % 8 + 1}",
            }
            for donor in range(donors)
            for lineage in adapter.LINEAGES
        ]

    def test_inner_validation_is_deterministic_and_donor_disjoint(self) -> None:
        rows = self._rows()
        first = adapter.inner_validation_indices(
            rows, seed=1103, validation_fraction=0.2
        )
        second = adapter.inner_validation_indices(
            rows, seed=1103, validation_fraction=0.2
        )
        self.assertEqual(first, second)
        training, validation = first
        self.assertFalse(
            {rows[index]["donor_id"] for index in training}
            & {rows[index]["donor_id"] for index in validation}
        )
        self.assertEqual(
            {rows[index]["lineage"] for index in training}, set(adapter.LINEAGES)
        )
        self.assertEqual(
            {rows[index]["lineage"] for index in validation}, set(adapter.LINEAGES)
        )

    @unittest.skipIf(np is None, "requires the locked scientific NumPy environment")
    def test_profile_normalization_is_finite_and_depth_free(self) -> None:
        values = np.asarray([[0.2, 0.3, 0.5], [2.0, 3.0, 5.0]])
        profiles = adapter.profile_normalize(values, 1e-8)
        np.testing.assert_allclose(profiles.sum(axis=1), 1.0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(profiles[0], profiles[1], rtol=0, atol=1e-8)

    @unittest.skipIf(torch is None, "requires the locked scientific PyTorch environment")
    def test_pinned_architecture_emits_distinct_valid_output_px(self) -> None:
        torch.manual_seed(1103)
        model = adapter._make_model(
            input_dim=12, output_dim=7, hidden_layers=[900, 40], dropout=0.1
        )
        model.eval()
        values = torch.arange(36, dtype=torch.float32).reshape(3, 12)
        output_result, output_px, embedding = model(values)
        self.assertEqual(tuple(output_result.shape), (3, 7))
        self.assertEqual(tuple(output_px.shape), (3, 7))
        self.assertEqual(tuple(embedding.shape), (3, 40))
        self.assertTrue(bool(torch.all(output_px >= 0)))
        self.assertTrue(bool(torch.all(output_px <= 1)))
        self.assertTrue(bool(torch.all(output_result <= output_px)))
        self.assertFalse(torch.equal(output_result, output_px))


if __name__ == "__main__":
    unittest.main()
