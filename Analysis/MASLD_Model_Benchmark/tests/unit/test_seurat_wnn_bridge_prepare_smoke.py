from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.sparse import csr_matrix, save_npz


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "seurat_wnn_bridge_prepare_smoke.py"
SPEC = importlib.util.spec_from_file_location("seurat_wnn_bridge_prepare_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
prepare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prepare)


class SeuratWNNBridgePrepareTests(unittest.TestCase):
    def test_load_matrix_rejects_negative_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.npz"
            save_npz(path, csr_matrix(np.array([[1, -1]], dtype=np.int32)))
            with self.assertRaisesRegex(
                prepare.SeuratWNNPreparationError, "matrix values"
            ):
                prepare.load_matrix(path)

    def test_load_matrix_preserves_axes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.npz"
            save_npz(path, csr_matrix(np.eye(3, dtype=np.int32)))
            observed = prepare.load_matrix(path)
            self.assertEqual(observed.shape, (3, 3))
            self.assertEqual(observed.nnz, 3)


if __name__ == "__main__":
    unittest.main()
