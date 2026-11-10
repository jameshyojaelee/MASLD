"""Read-only column lookup preserving original CSR summation order."""
import numpy as np
from scipy import sparse


class PromoterLookup:
    """Index original CSR data positions by column without changing their order."""
    def __init__(self, matrix):
        if not sparse.isspmatrix_csr(matrix):
            raise TypeError("Original CSR matrix required")
        self.matrix = matrix
        positions = sparse.csr_matrix(
            (np.arange(matrix.nnz, dtype=np.int64), matrix.indices, matrix.indptr),
            shape=matrix.shape, copy=False).tocsc()
        self.column_positions = positions

    def columns(self, columns):
        columns = np.asarray(columns, dtype=np.int64)
        if columns.ndim != 1 or len(set(columns.tolist())) != len(columns):
            raise ValueError("Column selection must be one-dimensional and unique")
        if np.any(columns < 0) or np.any(columns >= self.matrix.shape[1]):
            raise IndexError("Column outside original matrix")
        selected = self.column_positions[:, columns]
        positions = selected.data
        column_ids = np.repeat(np.arange(len(columns), dtype=self.matrix.indices.dtype), np.diff(selected.indptr))
        order = np.argsort(positions, kind="stable")
        positions, column_ids = positions[order], column_ids[order]
        rows = np.searchsorted(self.matrix.indptr, positions, side="right") - 1
        indptr = np.empty(self.matrix.shape[0] + 1, dtype=self.matrix.indptr.dtype)
        indptr[0] = 0
        np.cumsum(np.bincount(rows, minlength=self.matrix.shape[0]), out=indptr[1:])
        return sparse.csr_matrix((self.matrix.data[positions], column_ids, indptr),
                                 shape=(self.matrix.shape[0], len(columns)), copy=False)

    def promoter(self, columns):
        subset = self.columns(columns)
        return np.log1p(np.asarray(subset.sum(axis=1)).ravel())

    @property
    def index_bytes(self):
        return sum(value.nbytes for value in [self.column_positions.data,
                                              self.column_positions.indices,
                                              self.column_positions.indptr])
