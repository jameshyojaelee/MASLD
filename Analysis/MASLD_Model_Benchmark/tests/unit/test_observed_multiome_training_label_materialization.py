from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
from scipy import sparse

from scripts import materialize_gse296875_observed_multiome_training_labels as materializer


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_training_label_materialization_20260825.json"


class ObservedMultiomeTrainingLabelMaterializationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_config_passes_and_held_row_mutation_fails(self) -> None:
        resolved = materializer.validate_config(ROOT, deepcopy(self.config))
        self.assertEqual(set(resolved), {"supervised_training_contract", "training_target_mask_plan", "base_mask_plan", "input_h5"})
        mutated = deepcopy(self.config)
        mutated["artifact"]["held_donor_rows"] = "allowed"
        with self.assertRaises(materializer.TrainingLabelMaterializationError):
            materializer.validate_config(ROOT, mutated)

    def test_row_and_column_selective_read_excludes_held_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "matrix.h5"
            with h5py.File(path, "w") as handle:
                group = handle.create_group("csr")
                group.create_dataset("data", data=np.asarray([1, 2, 3, 4, 5, 6]))
                group.create_dataset("indices", data=np.asarray([0, 2, 1, 2, 0, 3]))
                group.create_dataset("indptr", data=np.asarray([0, 2, 4, 6]))
                group.create_dataset("shape", data=np.asarray([3, 4]))
            with h5py.File(path, "r") as handle:
                matrix, positions, source_rows = materializer.selectively_read_csr_rows_columns(handle["csr"], [0, 2], [2, 0])
            np.testing.assert_array_equal(matrix.toarray(), [[2, 1], [0, 5]])
            np.testing.assert_array_equal(positions, [0, 1, 4])
            np.testing.assert_array_equal(source_rows, [0, 0, 2])

    def test_training_h5_excludes_held_rows_and_nontraining_modalities(self) -> None:
        planned = [{"donor_hash": "a" * 64, "lineage": "Hepatocyte", "donor_fold": "1", "nuclei": "2"}]
        targets = [{"training_target_hash": "b" * 64, "chromosome": "chr2", "bed_start_0based": "10", "bed_end_half_open": "20", "source_genomic_fold": "1"}]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "training.h5"
            materializer._write_child(path, 0, 0, planned, ["c" * 64], sparse.csr_matrix([[3]]), targets)
            with h5py.File(path, "r") as handle:
                self.assertEqual(set(handle), {"rows", "training_target_atac"})
                self.assertTrue(bool(handle.attrs["training_only"]))
                self.assertNotIn("rna", handle)
                self.assertNotIn("observed_atac_input", handle)
                self.assertNotIn("evaluator", handle)
                self.assertNotIn(0, handle["rows/donor_fold"][:])

    def test_source_never_loads_full_atac_or_other_modality_values(self) -> None:
        source = inspect.getsource(materializer)
        self.assertNotIn('group["data"][:]', source)
        self.assertNotIn('handle["rna/', source)
        self.assertNotIn("evaluator_outcomes.h5", source)
        self.assertNotIn("model_input.h5", source)


if __name__ == "__main__":
    unittest.main()
