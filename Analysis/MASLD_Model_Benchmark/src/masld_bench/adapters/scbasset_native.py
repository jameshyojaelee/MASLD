"""Pinned scBasset architecture helpers for project-local training."""

from __future__ import annotations

from hashlib import sha256
import importlib.util
from pathlib import Path
import sys
from types import ModuleType
from typing import Any


SOURCE_REVISION = "aed3a6f713091fd988196b297e9c06e092ff1d22"
UTILS_SHA256 = "7b653646fc83867cf38c51ae14fdcdf523f0dc99685d1b9446ea5940de71543c"
BASENJI_SHA256 = "165fdc3c2dd85afd300bf71e27199ae907b005e325c33fa0cc1e9b96be9a92b2"
SEQUENCE_LENGTH = 1344
BOTTLENECK_DIMENSION = 32


class ScBassetNativeError(RuntimeError):
    """Raised when the pinned native architecture or tensors differ."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_basenji_helpers(source: Path) -> Any:
    source = source.resolve(strict=True)
    utils_path = source / "scbasset" / "utils.py"
    basenji_path = source / "scbasset" / "basenji_utils.py"
    if sha256_file(utils_path) != UTILS_SHA256 or sha256_file(basenji_path) != BASENJI_SHA256:
        raise ScBassetNativeError("pinned scBasset source hashes differ")
    # pysam is imported only for optional FASTA extraction. Project sequence
    # materialization is separate and no production dependency is added here.
    if "pysam" not in sys.modules:
        sys.modules["pysam"] = ModuleType("pysam")
    spec = importlib.util.spec_from_file_location(
        "masld_scbasset_basenji_pinned", basenji_path
    )
    if spec is None or spec.loader is None:
        raise ScBassetNativeError("cannot load pinned scBasset Basenji helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_model(source: Path, n_training_cells: int) -> Any:
    import tensorflow as tf

    if isinstance(n_training_cells, bool) or n_training_cells < 1:
        raise ScBassetNativeError("training-cell output count must be positive")
    basenji = load_basenji_helpers(source)
    sequence = tf.keras.Input(shape=(SEQUENCE_LENGTH, 4), name="sequence")
    current, reverse_bool = basenji.StochasticReverseComplement()(sequence)
    current = basenji.StochasticShift(3)(current)
    current = basenji.conv_block(
        current, filters=288, kernel_size=17, pool_size=3
    )
    current = basenji.conv_tower(
        current,
        filters_init=288,
        filters_mult=1.122,
        repeat=6,
        kernel_size=5,
        pool_size=2,
    )
    current = basenji.conv_block(current, filters=256, kernel_size=1)
    current = basenji.dense_block(
        current, flatten=True, units=BOTTLENECK_DIMENSION, dropout=0.2
    )
    current = basenji.GELU()(current)
    current = basenji.final(
        current, units=n_training_cells, activation="sigmoid"
    )
    current = basenji.SwitchReverse()([current, reverse_bool])
    current = tf.keras.layers.Flatten()(current)
    model = tf.keras.Model(inputs=sequence, outputs=current)
    if tuple(model.input_shape) != (None, SEQUENCE_LENGTH, 4) or tuple(
        model.output_shape
    ) != (None, n_training_cells):
        raise ScBassetNativeError("native scBasset model shape differs")
    return model


def one_hot_base_codes(codes: Any) -> Any:
    import numpy as np

    values = np.asarray(codes)
    if values.ndim != 2 or values.shape[1] != SEQUENCE_LENGTH:
        raise ScBassetNativeError("base-code tensor shape differs")
    if not np.issubdtype(values.dtype, np.integer) or np.any(values < 0) or np.any(values > 3):
        raise ScBassetNativeError("base codes must be integer A/C/G/T indices")
    return np.eye(4, dtype=np.float32)[values]


def hdf5_tensor_manifest(path: Path) -> list[dict[str, Any]]:
    import h5py
    import numpy as np

    records: list[dict[str, Any]] = []
    with h5py.File(path, "r") as handle:
        def record(name: str, value: object) -> None:
            if isinstance(value, h5py.Dataset):
                array = np.asarray(value[()])
                records.append(
                    {
                        "name": name,
                        "shape": list(array.shape),
                        "dtype": str(array.dtype),
                        "sha256": sha256(array.tobytes(order="C")).hexdigest(),
                    }
                )
        handle.visititems(record)
    if not records:
        raise ScBassetNativeError("HDF5 checkpoint contains no tensor members")
    return records
