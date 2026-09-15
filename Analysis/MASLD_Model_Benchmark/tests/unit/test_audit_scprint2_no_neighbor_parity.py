from __future__ import annotations

import unittest

import numpy as np
from scipy import sparse

from scripts.audit_scprint2_no_neighbor_parity import row_totals


class Scprint2NoNeighborParityAuditTests(unittest.TestCase):
    def test_row_totals_preserve_source_depth_unit(self) -> None:
        matrix = sparse.csr_matrix(np.asarray([[1, 2, 0], [0, 3, 4]]))
        np.testing.assert_array_equal(row_totals(matrix), np.asarray([3.0, 7.0]))


if __name__ == "__main__":
    unittest.main()
