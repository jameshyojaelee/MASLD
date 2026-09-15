from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
import tempfile
import unittest

import h5py
from scipy import sparse

from scripts import materialize_gse296875_observed_multiome_evaluator_outcomes as materializer


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_evaluator_materialization_20260825.json"


class ObservedMultiomeEvaluatorMaterializationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_config_passes_and_firewall_mutation_fails(self) -> None:
        resolved = materializer.validate_config(ROOT, deepcopy(self.config))
        self.assertEqual(set(resolved), {"materialization_contract", "mask_plan", "verified_model_input_control", "input_h5"})
        mutated = deepcopy(self.config)
        mutated["firewall"]["observed_atac_input_count_values_read"] = True
        with self.assertRaises(materializer.EvaluatorMaterializationError):
            materializer.validate_config(ROOT, mutated)

    def test_evaluator_h5_contains_only_rows_and_target_counts(self) -> None:
        planned = [{"donor_hash": "a" * 64, "lineage": "Hepatocyte", "donor_fold": "0", "nuclei": "2"}]
        targets = [{"target_hash": "b" * 64, "chromosome": "chr1", "bed_start_0based": "10", "bed_end_half_open": "20"}]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evaluator.h5"
            materializer._write_fold_h5(path, 0, planned, ["c" * 64], sparse.csr_matrix([[3]]), targets)
            with h5py.File(path, "r") as handle:
                self.assertEqual(set(handle), {"rows", "target_atac"})
                self.assertTrue(bool(handle.attrs["evaluator_only"]))
                self.assertNotIn("rna", handle)
                self.assertNotIn("observed_atac_input", handle)
                self.assertIn("counts_csr", handle["target_atac"])

    def test_source_reads_only_selected_target_atac_values(self) -> None:
        source = inspect.getsource(materializer)
        self.assertNotIn("read_full_csr", source)
        self.assertIn("selectively_read_csr_columns(atac_group, target_columns)", source)
        self.assertNotIn("selectively_read_csr_columns(atac_group, input_columns)", source)
        self.assertNotIn('handle["rna/', source)
        self.assertNotIn('atac_group["data"][:]', source)


if __name__ == "__main__":
    unittest.main()
