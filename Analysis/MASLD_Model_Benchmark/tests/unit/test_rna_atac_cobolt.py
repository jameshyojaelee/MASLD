from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ADAPTER_PATH = PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_cobolt.py"
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_cobolt_unit", ADAPTER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load standalone adapter: {ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)


class RNAATACCoboltUnitTests(unittest.TestCase):
    def test_inner_split_is_donor_disjoint_and_deterministic(self) -> None:
        rows = [
            {"cell_id": f"cell-{donor}-{lineage}", "donor_id": donor, "lineage": lineage}
            for donor in [f"donor-{index}" for index in range(10)]
            for lineage in adapter.LINEAGES
        ]
        first = adapter.inner_validation_indices(rows, seed=1103, validation_fraction=0.2)
        second = adapter.inner_validation_indices(rows, seed=1103, validation_fraction=0.2)
        self.assertEqual(first, second)
        training, validation = first
        training_donors = {rows[index]["donor_id"] for index in training}
        validation_donors = {rows[index]["donor_id"] for index in validation}
        self.assertFalse(training_donors & validation_donors)

    def test_parameter_contract_rejects_dataset_adjustments(self) -> None:
        parameters = {
            "alpha": 5.0,
            "annealing_epochs": 30,
            "batch_size": 128,
            "early_stopping_patience": 20,
            "hidden_dims": [128, 64],
            "inference_batch_size": 1,
            "intercept_adjustment": False,
            "join_namespace": "namespace",
            "learning_rate": 0.005,
            "max_epochs": 100,
            "n_hvg": 2000,
            "n_latent": 10,
            "n_smoke_peaks": 2000,
            "outer_folds": 5,
            "slope_adjustment": False,
            "split_seed": 20260821,
            "validation_fraction": 0.2,
        }
        adapter._validate_parameters(parameters)
        parameters["slope_adjustment"] = True
        with self.assertRaises(adapter.RNAATACCoboltError):
            adapter._validate_parameters(parameters)


if __name__ == "__main__":
    unittest.main()
