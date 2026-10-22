"""Numerical helpers for the fixed H48 knockout-signature comparison."""

from pathlib import Path

import numpy as np


def empirical_two_sided(observed, null):
    """Add-one tail probability for a finite statistic and a complete null."""
    null = np.asarray(null, dtype=float)
    if null.ndim != 1 or null.size == 0:
        raise ValueError("The permutation null must be a nonempty vector")
    if not np.isfinite(observed) or not np.isfinite(null).all():
        return np.nan
    return (1 + np.count_nonzero(np.abs(null) >= abs(observed))) / (null.size + 1)


def fraction_below(values, threshold):
    """An incomplete p-value family cannot establish its rejection fraction."""
    values = np.asarray(values, dtype=float)
    if values.size == 0 or not np.isfinite(values).all():
        return np.nan
    return float(np.mean(values < threshold))


def read_null_array(path, shape, index=None):
    """Map a float32 null, optionally restricting resident data to one knockout."""
    expected = int(np.prod(shape, dtype=np.int64)) * np.dtype(np.float32).itemsize
    if Path(path).stat().st_size != expected:
        raise ValueError(f"Null array size does not match its declared dimensions: {path}")
    offset = 0
    if index is not None:
        if not 0 <= index < shape[0]:
            raise IndexError("Null knockout index is out of range")
        offset = int(index) * int(np.prod(shape[1:], dtype=np.int64)) * np.dtype(np.float32).itemsize
        shape = shape[1:]
    return np.memmap(path, mode="r", dtype=np.float32, offset=offset, shape=tuple(shape))
