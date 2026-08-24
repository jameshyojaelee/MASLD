import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

from masld_cl.config import load_config
from masld_cl.gse296875_router import (
    GSE296875RouterError,
    load_10x_canonical_subset,
    load_10x_raw_subset,
    load_gse296875_mapping_policy,
    load_gse296875_source_lock,
)


class TestGSE296875Router(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.config = load_config(cls.root / "config_v1.json")
        cls.policy_path = cls.root / "reference" / "gse296875_source_lock_v44.json"
        cls.mapping_policy_path = cls.root / "reference" / "gse296875_mapping_policy_v45.json"

    def _fixture(self, directory: str) -> Path:
        path = Path(directory) / "raw.h5"
        dense = np.asarray([
            [1, 2],
            [3, 4],
            [5, 0],
            [7, 8],
        ], dtype=np.int32)
        matrix = sparse.csc_matrix(dense)
        with h5py.File(path, "w") as handle:
            group = handle.create_group("matrix")
            group.create_dataset("barcodes", data=np.asarray([b"a", b"b"]))
            group.create_dataset("data", data=matrix.data)
            group.create_dataset("indices", data=matrix.indices)
            group.create_dataset("indptr", data=matrix.indptr)
            group.create_dataset("shape", data=np.asarray(matrix.shape))
            features = group.create_group("features")
            features.create_dataset("id", data=np.asarray([b"i1", b"i2", b"i3", b"p1"]))
            features.create_dataset("name", data=np.asarray([b"G1", b"G1", b"G2", b"P1"]))
            features.create_dataset(
                "feature_type",
                data=np.asarray([b"Gene Expression", b"Gene Expression", b"Gene Expression", b"Peaks"]),
            )
        return path

    def test_subset_order_and_duplicate_gene_aggregation(self):
        with tempfile.TemporaryDirectory() as directory:
            result = load_10x_raw_subset(
                self._fixture(directory), ["b", "a"], ["G2", "G1", "ABSENT"],
            )
        np.testing.assert_array_equal(result["genes"], ["G2", "G1"])
        np.testing.assert_array_equal(result["counts"].toarray(), [[0, 6], [5, 4]])
        self.assertEqual(result["source_rna_features"], 3)

    def test_feature_id_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            result = load_10x_raw_subset(
                self._fixture(directory), ["a"], ["i3"], match_feature_ids=True,
            )
        np.testing.assert_array_equal(result["genes"], ["i3"])
        np.testing.assert_array_equal(result["counts"].toarray(), [[5]])

    def test_v40_one_to_one_mapping_excludes_symbol_collision(self):
        with tempfile.TemporaryDirectory() as directory:
            result = load_10x_canonical_subset(
                self._fixture(directory), ["a"], ["G1", "i3"],
            )
        np.testing.assert_array_equal(result["genes"], ["i3"])
        np.testing.assert_array_equal(result["counts"].toarray(), [[5]])
        self.assertEqual(result["mapping_methods"], {"external_var_id": 1})

    def test_missing_barcode_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(GSE296875RouterError):
                load_10x_raw_subset(self._fixture(directory), ["missing"], ["G1"])

    def test_external_source_lock_and_copy_failure(self):
        _, policy, processed, raw = load_gse296875_source_lock(
            self.config, self.policy_path,
        )
        self.assertEqual(processed.stat().st_size, policy["processed_source"]["bytes"])
        self.assertEqual(set(raw), {f"well{index}" for index in range(1, 9)})
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / "policy.json"
            copied.write_bytes(self.policy_path.read_bytes())
            with self.assertRaises(GSE296875RouterError):
                load_gse296875_source_lock(self.config, copied)

    def test_mapping_policy_is_locked_before_author_values(self):
        _, policy, _, _, _ = load_gse296875_mapping_policy(
            self.config, self.mapping_policy_path,
        )
        self.assertFalse(policy["structure"]["metadata_values_opened"])
        self.assertEqual(
            policy["author_label_evaluation"]["primary_column"],
            "celltype_dXRNA_res.0.25_harmony",
        )
        self.assertEqual(policy["validation"]["minimum_donor_balanced_macro_f1"], 0.70)


if __name__ == "__main__":
    unittest.main()
