"""Sparse, donor-first helpers for isolated observed-multiome output files."""

from __future__ import annotations

from hashlib import sha256
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


class ObservedMultiomeMaterializationError(ValueError):
    """Raised when sparse selection, aggregation, or identifiers drift."""


def salted_hash(namespace: str, *parts: object) -> str:
    if not namespace or not parts or any(str(part) == "" for part in parts):
        raise ObservedMultiomeMaterializationError("invalid hash input")
    return sha256(
        "\0".join((namespace, *(str(part) for part in parts))).encode("utf-8")
    ).hexdigest()


def decode_strings(values: Iterable[Any], *, label: str) -> list[str]:
    result: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        text = str(value)
        if not text:
            raise ObservedMultiomeMaterializationError(f"empty {label}")
        result.append(text)
    return result


def validate_csr_structure(
    indices: np.ndarray, indptr: np.ndarray, shape: Sequence[int]
) -> tuple[int, int]:
    rows, columns = map(int, shape)
    indices = np.asarray(indices)
    indptr = np.asarray(indptr)
    if (
        rows < 1
        or columns < 1
        or indices.ndim != 1
        or indptr.shape != (rows + 1,)
        or int(indptr[0]) != 0
        or int(indptr[-1]) != indices.size
        or np.any(indptr[1:] < indptr[:-1])
        or np.any(indices < 0)
        or np.any(indices >= columns)
    ):
        raise ObservedMultiomeMaterializationError("CSR structure differs")
    return rows, columns


def selected_value_positions(
    indices: np.ndarray, selected_columns: Sequence[int], *, total_columns: int
) -> np.ndarray:
    selected = np.asarray(selected_columns, dtype=np.int64)
    if (
        selected.ndim != 1
        or selected.size < 1
        or len(set(map(int, selected))) != selected.size
        or np.any(selected < 0)
        or np.any(selected >= total_columns)
    ):
        raise ObservedMultiomeMaterializationError("selected CSR columns differ")
    membership = np.zeros(total_columns, dtype=bool)
    membership[selected] = True
    return np.flatnonzero(membership[np.asarray(indices, dtype=np.int64)])


def selectively_read_csr_columns(
    group: Any, selected_columns: Sequence[int]
) -> tuple[Any, np.ndarray]:
    """Read values only at sparse positions belonging to selected columns."""

    from scipy import sparse

    shape = tuple(int(value) for value in group["shape"][:])
    indices = np.asarray(group["indices"][:], dtype=np.int64)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    rows, columns = validate_csr_structure(indices, indptr, shape)
    selected = np.asarray(selected_columns, dtype=np.int64)
    positions = selected_value_positions(indices, selected, total_columns=columns)
    remap = np.full(columns, -1, dtype=np.int64)
    remap[selected] = np.arange(selected.size, dtype=np.int64)
    selected_indices = remap[indices[positions]]
    keep = np.zeros(indices.size, dtype=np.int64)
    keep[positions] = 1
    prefix = np.concatenate((np.zeros(1, dtype=np.int64), np.cumsum(keep)))
    row_counts = prefix[indptr[1:]] - prefix[indptr[:-1]]
    selected_indptr = np.concatenate(
        (np.zeros(1, dtype=np.int64), np.cumsum(row_counts))
    )
    data = np.asarray(group["data"][positions])
    matrix = sparse.csr_matrix(
        (data, selected_indices, selected_indptr),
        shape=(rows, selected.size),
    )
    if matrix.nnz != positions.size or np.any(matrix.data < 0):
        raise ObservedMultiomeMaterializationError("selected CSR values differ")
    return matrix, positions


def read_full_csr(group: Any) -> Any:
    from scipy import sparse

    shape = tuple(int(value) for value in group["shape"][:])
    indices = np.asarray(group["indices"][:], dtype=np.int64)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    validate_csr_structure(indices, indptr, shape)
    data = np.asarray(group["data"][:])
    matrix = sparse.csr_matrix((data, indices, indptr), shape=shape)
    if np.any(matrix.data < 0):
        raise ObservedMultiomeMaterializationError("CSR values are negative")
    return matrix


def donor_lineage_aggregation(
    *,
    cell_unit_keys: Sequence[tuple[str, str]],
    planned_units: Sequence[Mapping[str, str]],
) -> tuple[Any, list[str]]:
    from scipy import sparse

    unit_keys = [
        (str(row["donor_hash"]), str(row["lineage"])) for row in planned_units
    ]
    if len(unit_keys) != len(set(unit_keys)) or not unit_keys:
        raise ObservedMultiomeMaterializationError("planned units differ")
    by_key = {key: index for index, key in enumerate(unit_keys)}
    row_indices: list[int] = []
    for key in cell_unit_keys:
        if key not in by_key:
            raise ObservedMultiomeMaterializationError("cell unit is absent from plan")
        row_indices.append(by_key[key])
    aggregation = sparse.csr_matrix(
        (
            np.ones(len(row_indices), dtype=np.int64),
            (np.asarray(row_indices), np.arange(len(row_indices))),
        ),
        shape=(len(unit_keys), len(row_indices)),
    )
    counts = np.asarray(aggregation.sum(axis=1)).ravel().astype(int)
    expected = np.asarray([int(row["nuclei"]) for row in planned_units])
    if not np.array_equal(counts, expected) or np.any(counts < 1):
        raise ObservedMultiomeMaterializationError("unit nucleus census differs")
    return aggregation, [f"{donor_hash}\0{lineage}" for donor_hash, lineage in unit_keys]


def write_string_dataset(group: Any, name: str, values: Sequence[str]) -> None:
    import h5py

    group.create_dataset(
        name,
        data=list(values),
        dtype=h5py.string_dtype(encoding="utf-8"),
        compression="gzip",
    )


def write_csr(group: Any, matrix: Any) -> dict[str, Any]:
    value = matrix.tocsr()
    group.create_dataset("data", data=value.data, compression="gzip", shuffle=True)
    group.create_dataset(
        "indices", data=value.indices.astype(np.int64), compression="gzip", shuffle=True
    )
    group.create_dataset(
        "indptr", data=value.indptr.astype(np.int64), compression="gzip", shuffle=True
    )
    group.create_dataset("shape", data=np.asarray(value.shape, dtype=np.int64))
    return {
        "shape": list(value.shape),
        "nnz": int(value.nnz),
        "sum": int(value.sum()),
        "data_sha256": sha256(np.asarray(value.data).tobytes()).hexdigest(),
        "indices_sha256": sha256(value.indices.astype(np.int64).tobytes()).hexdigest(),
        "indptr_sha256": sha256(value.indptr.astype(np.int64).tobytes()).hexdigest(),
    }
