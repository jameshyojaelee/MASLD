from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

try:
    import numpy as np
except ModuleNotFoundError:  # Control-plane Python intentionally has no NumPy.
    np = None  # type: ignore[assignment]

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ADAPTER_PATH = (
    PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_regularized.py"
)
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_regularized_unit", ADAPTER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load standalone scientific adapter: {ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)

LINEAGES = adapter.LINEAGES
choose_shrinkage_weight = adapter.choose_shrinkage_weight
context_permutation_indices = adapter.context_permutation_indices


class RNAATACRegularizedTests(unittest.TestCase):
    def _rows(self, donors: int = 9) -> list[dict[str, str]]:
        return [
            {"donor_id": f"donor-{donor:02d}", "lineage": lineage}
            for donor in range(donors)
            for lineage in LINEAGES
        ]

    @unittest.skipIf(np is None, "requires the locked scientific NumPy environment")
    def test_context_permutation_is_deterministic_lineage_stratified_derangement(
        self,
    ) -> None:
        rows = self._rows()
        first = context_permutation_indices(rows, seed=20260824)
        second = context_permutation_indices(rows, seed=20260824)
        self.assertEqual(first, second)
        self.assertEqual(sorted(first), list(range(len(rows))))
        for destination, source in enumerate(first):
            self.assertNotEqual(destination, source)
            self.assertEqual(rows[destination]["lineage"], rows[source]["lineage"])

    @unittest.skipIf(np is None, "requires the locked scientific NumPy environment")
    def test_shrinkage_selection_uses_grouped_training_profiles(self) -> None:
        rows = self._rows()
        common = np.linspace(1.0, 2.0, 12)
        common /= common.sum()
        profiles = np.vstack([common for _ in rows])
        selected, records = choose_shrinkage_weight(
            profiles, rows, pseudocount=1e-8
        )
        self.assertEqual(selected, 0.0)
        self.assertEqual(len(records), 6)
        self.assertTrue(all(int(row["validation_profiles"]) == len(rows) for row in records))


if __name__ == "__main__":
    unittest.main()
