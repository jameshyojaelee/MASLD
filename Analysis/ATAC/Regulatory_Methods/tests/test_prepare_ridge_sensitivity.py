"""Scientific identity checks for read-only ridge input reconstruction."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

source = Path(__file__).resolve().parents[1] / "scripts/prepare_ridge_sensitivity.py"
spec = importlib.util.spec_from_file_location("prepare_ridge", source)
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)


def table_digest(frame):
    return hashlib.sha256(frame.to_csv(sep="\t", index=False).encode()).hexdigest()


def matrix_digest(matrix):
    value = matrix.tocsr()
    return hashlib.sha256(value.toarray().tobytes()).hexdigest()


class ReconstructionIdentity(unittest.TestCase):
    def test_contract_rejects_each_changed_scientific_axis(self):
        core = SimpleNamespace(table_digest=table_digest, matrix_digest=matrix_digest)
        counts = sparse.csr_matrix([[1., 0], [2, 1]])
        metadata = pd.DataFrame(dict(donor_id=["a", "b"], metacell_id=["a0", "b0"]))
        peaks = pd.DataFrame(dict(peak_id=["chr1:1-2", "chr1:4-5"]))
        genes = pd.DataFrame(dict(gene_id=["g1", "g2"]))
        rna = sparse.csr_matrix([[3., 4], [5, 6]])
        folds = metadata.assign(outer_fold=[0, 1])
        objects = [counts, metadata, peaks, genes, rna, genes.copy(), folds]
        keys = ["counts_sha256", "metadata_sha256", "peaks_sha256", "genes_sha256",
                "rna_sha256", "rna_genes_sha256", "fold_assignments_sha256"]
        functions = [matrix_digest, table_digest, table_digest, table_digest,
                     matrix_digest, table_digest, table_digest]
        contract = {key: function(value) for key, function, value in zip(keys, functions, objects)}
        self.assertEqual(prep.verify_contract(core, contract, *objects), contract)
        for index, key in enumerate(keys):
            altered = list(objects)
            altered[index] = objects[index][::-1].copy() if sparse.issparse(objects[index]) else objects[index].iloc[::-1].copy()
            with self.subTest(axis=key), self.assertRaisesRegex(ValueError, key):
                prep.verify_contract(core, contract, *altered)

    def test_exact_membership_and_pairing_order(self):
        actual = pd.DataFrame(dict(donor_id=["a", "a"], cell_id=["c1", "c2"], source_index=[0, 1]))
        prep.same_frame(actual, actual.astype(str), list(actual), "membership")
        with self.assertRaisesRegex(ValueError, "membership"):
            prep.same_frame(actual.iloc[::-1], actual, list(actual), "membership")
        changed = actual.copy()
        changed.loc[0, "donor_id"] = "b"
        with self.assertRaisesRegex(ValueError, "membership"):
            prep.same_frame(changed, actual, list(actual), "membership")

    def test_refuses_existing_output_before_loading_or_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "protected"
            marker.write_text("preserve")
            before = {str(path): path.read_bytes() for path in root.iterdir()}
            with self.assertRaises(FileExistsError):
                prep.prepare(root, root, root, root, root)
            self.assertEqual(before, {str(path): path.read_bytes() for path in root.iterdir()})

    def test_changed_snapshot_rejected_without_import(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "source_snapshot").mkdir()
            snapshot = root / "source_snapshot/network_compare.py"
            snapshot.write_text("raise AssertionError('must not import changed source')\n")
            (root / "completion.json").write_text(json.dumps(dict(status="completed", pilot=False)))
            (root / "source_sha256.json").write_text(json.dumps({"/original/network_compare.py": "incorrect"}))
            with self.assertRaisesRegex(ValueError, "snapshot hash"):
                prep.load_core(root)

    def test_rejects_output_under_original_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, "separate from original"):
                prep.prepare(root, root, root, root, root / "new_sensitivity")
            self.assertFalse((root / "new_sensitivity").exists())


if __name__ == "__main__":
    unittest.main()
