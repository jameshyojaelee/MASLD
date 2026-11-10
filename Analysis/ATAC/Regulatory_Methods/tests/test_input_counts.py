"""Small identity/count invariants; run inside the existing compute allocation."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/input_counts.py"
SPEC = importlib.util.spec_from_file_location("input_counts", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class CountInvariants(unittest.TestCase):
    def test_fragment_and_detected_cell_counts_are_distinct(self):
        matrix = sp.csr_matrix([[3, 0, 1], [2, 1, 0], [0, 4, 0]])
        raw, detected = module.aggregate_counts(matrix, np.array([0, 0, 1]), 2)
        np.testing.assert_array_equal(raw.toarray(), [[5, 1, 1], [0, 4, 0]])
        np.testing.assert_array_equal(detected.toarray(), [[2, 1, 1], [0, 1, 0]])

    def test_chunk_aggregation_conserves_all_donor_peak_counts(self):
        matrix = sp.csr_matrix([[1, 0], [3, 2], [0, 4], [7, 0]])
        codes = np.array([0, 1, 0, 1])
        whole = module.aggregate_counts(matrix, codes, 2)
        left = module.aggregate_counts(matrix[:2], codes[:2], 2)
        right = module.aggregate_counts(matrix[2:], codes[2:], 2)
        for total, first, second in zip(whole, left, right):
            self.assertEqual((total != first + second).nnz, 0)

    def test_noninteger_counts_are_rejected(self):
        with self.assertRaises(ValueError):
            module.aggregate_counts(sp.csr_matrix([[0.5]]), np.array([0]), 1)
        with self.assertRaises(ValueError):
            module.aggregate_counts(sp.csr_matrix([[-1]]), np.array([0]), 1)

    def test_duplicate_sparse_entries_are_one_detected_cell(self):
        matrix = sp.csr_matrix((np.array([2, 3]), np.array([0, 0]), np.array([0, 2])), shape=(1, 1))
        raw, detected = module.aggregate_counts(matrix, np.array([0]), 1)
        self.assertEqual(int(raw[0, 0]), 5)
        self.assertEqual(int(detected[0, 0]), 1)

    def test_pilot_retains_100_cells_per_selected_donor(self):
        donors = np.repeat([f"D{x:02}" for x in range(18)], 175)
        first = module.pilot_indices(donors, 1000, 42)
        second = module.pilot_indices(donors, 1000, 42)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(first), 1000)
        self.assertEqual(len(np.unique(first)), 1000)
        _, counts = np.unique(donors[first], return_counts=True)
        self.assertEqual(len(counts), 10)
        self.assertTrue(np.all(counts >= 100))

    def test_pilot_excludes_underpopulated_donors(self):
        donors = np.array(["small"] * 20 + ["large"] * 150)
        selected = module.pilot_indices(donors, 1000, 42)
        self.assertEqual(set(donors[selected]), {"large"})
        self.assertEqual(len(selected), 150)

    def test_balancing_samples_whole_cells_and_preserves_primary_roster(self):
        roster = pd.DataFrame({"donor_id": ["D1"] * 25 + ["D2"] * 40 + ["MASL"] * 30,
                               "da_eligible": [True] * 65 + [False] * 30})
        mask, target = module.balanced_cell_mask(roster, 42)
        self.assertEqual(target, 25)
        self.assertEqual(roster.loc[mask].groupby("donor_id").size().to_dict(), {"D1": 25, "D2": 25})
        self.assertEqual(mask.sum(), 50)
        np.testing.assert_array_equal(mask, module.balanced_cell_mask(roster, 42)[0])


class IdentityInvariants(unittest.TestCase):
    def test_rna_duplicate_library_barcode_is_not_independent_cell(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "donors.csv"
            metadata.write_text("donor_id,condition,rna_srrs,rna_gsm\nD01,NORMAL,SRR1;SRR2,GSM1\n")
            h5ad = root / "rna.h5ad"
            with h5py.File(h5ad, "w") as handle:
                obs = handle.create_group("obs")
                dtype = h5py.string_dtype()
                obs.create_dataset("_index", data=np.array(["SRR1_AAAC-1", "SRR2_AAAC-1"], dtype=object), dtype=dtype)
                obs.create_dataset("sample", data=np.array(["SRR1", "SRR2"], dtype=object), dtype=dtype)
                obs.create_dataset("condition", data=np.array(["MASLD", "MASLD"], dtype=object), dtype=dtype)
            table, summary = module.audit_rna({"project_root": str(root), "rna_reference": str(h5ad), "cohorts": {"GSE244832": {"metadata": str(metadata)}}})
            self.assertEqual(summary["duplicate_library_barcode_rows"], 2)
            self.assertEqual(summary["unique_resolved_library_barcodes"], 1)
            self.assertFalse(summary["rna_model_eligible"])
            self.assertFalse(summary["expression_matrix_read"])
            self.assertEqual(table["donor_condition"].tolist(), ["NORMAL", "NORMAL"])

    def test_missing_categorical_identity_is_empty(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "metadata.h5"
            with h5py.File(path, "w") as handle:
                group = handle.create_group("obs")
                value = group.create_group("donor")
                value.create_dataset("categories", data=np.array(["D01"], dtype=object), dtype=h5py.string_dtype())
                value.create_dataset("codes", data=[0, -1])
                self.assertEqual(module.obs_column(group, "donor").tolist(), ["D01", ""])

    def test_json_outputs_refuse_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "diagnostics.json"
            module.write_json(path, {"count": 2})
            with self.assertRaises(FileExistsError):
                module.write_json(path, {"count": 3})
            self.assertEqual(json.loads(path.read_text()), {"count": 2})

    def test_bounded_rna_vector_check_does_not_claim_rawness(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "rna.h5ad"
            with h5py.File(source, "w") as handle:
                matrix = handle.create_group("X")
                matrix.attrs["encoding-type"] = "csr_matrix"
                matrix.attrs["shape"] = [2, 2]
                matrix.create_dataset("data", data=[3.0, 3.0])
                matrix.create_dataset("indices", data=[0, 0])
                matrix.create_dataset("indptr", data=[0, 1, 2])
            pd.DataFrame({"rna_library_gsm": ["GSM1", "GSM1"], "raw_barcode": ["AA-1", "AA-1"],
                          "rna_run": ["SRR1", "SRR2"], "duplicate_library_barcode": [True, True]}).to_csv(root / "GSE244832.rna_identity.tsv.gz", sep="\t", index=False)
            with patch.dict("os.environ", {"SLURM_JOB_ID": "1", "SLURM_JOB_NODELIST": "fixture"}):
                module.audit_rna_vectors({"rna_reference": str(source)}, root)
            result = json.loads((root / "rna_duplicate_vector_check.json").read_text())
            self.assertEqual(result["groups_checked"], 1)
            self.assertTrue(result["all_groups_identical"])
            self.assertFalse(result["rna_model_eligible"])
            self.assertFalse(result["rows_deleted_or_summed"])
            self.assertFalse(result["raw_counts_verified"])


if __name__ == "__main__":
    unittest.main()
