"""Streaming readers for a pretraining corpus that does not fit in memory.

The corpus this method was developed against is 581,503 x 20,010 float32 -- 46 GB if
materialised. Nothing in this file ever holds more than `chunk_rows` rows at a time, so
peak memory is set by the chunk size and the query, never by the corpus. That is not an
optimisation; it is the difference between the tool running at all and not.

Four containers are supported, all row-major over corpus samples:

  .h5ad / .h5 / .hdf5   an HDF5 dataset, dense or AnnData CSR. Gene identifiers are read
                        from `var/` when present, including AnnData's categorical
                        encoding (a `categories` + `codes` group rather than a string
                        array).
  .npy                  memory-mapped, so a row slice touches only those rows.
  .npz                  decompressed on open. Convenient for fixtures, wrong for a real
                        corpus; the reader says so rather than silently using 46 GB.
  .tsv / .csv (+.gz)    read line by line, never wholly.

Every reader exposes the same three things: `n_rows`, `gene_ids` (or None), and
`iter_chunks()`. `take_rows()` gathers a scattered set of rows for the background probes
and for block confirmation, and is the only method allowed to read out of order.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np

__all__ = ["open_corpus", "CorpusReader"]


def _decode(a) -> list[str]:
    """Bytes, numpy string scalars and objects all become plain str."""
    out = []
    for v in a:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8"))
        elif isinstance(v, np.bytes_):
            out.append(v.tobytes().decode("utf-8"))
        else:
            out.append(str(v))
    return out


class CorpusReader:
    """Row-streaming view of a corpus matrix. Subclasses fill in the three hooks."""

    n_rows: int
    n_cols: int
    gene_ids: list[str] | None
    source: str

    def iter_chunks(self, chunk_rows: int) -> Iterator[tuple[int, int, np.ndarray]]:
        raise NotImplementedError

    def take_rows(self, rows: Sequence[int]) -> np.ndarray:
        raise NotImplementedError

    def close(self) -> None:
        pass

    def __enter__(self) -> "CorpusReader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def describe(self) -> dict:
        return {
            "source": self.source,
            "n_rows": int(self.n_rows),
            "n_cols": int(self.n_cols),
            "carries_gene_ids": self.gene_ids is not None,
        }


class _DenseArrayReader(CorpusReader):
    """Anything that supports `A[start:stop, :]` and `A[list_of_rows, :]`."""

    def __init__(self, arr, source: str, gene_ids: list[str] | None):
        self._a = arr
        self.n_rows, self.n_cols = int(arr.shape[0]), int(arr.shape[1])
        self.gene_ids = gene_ids
        self.source = source

    def iter_chunks(self, chunk_rows: int) -> Iterator[tuple[int, int, np.ndarray]]:
        for s in range(0, self.n_rows, chunk_rows):
            e = min(s + chunk_rows, self.n_rows)
            yield s, e, np.asarray(self._a[s:e, :], dtype=np.float32)

    def take_rows(self, rows: Sequence[int]) -> np.ndarray:
        rows = np.asarray(rows, dtype=np.int64)
        # h5py needs a sorted, increasing selection; restore the caller's order after.
        order = np.argsort(rows, kind="stable")
        got = np.asarray(self._a[rows[order], :], dtype=np.float32)
        back = np.empty_like(got)
        back[order] = got
        return back


class _H5Reader(_DenseArrayReader):
    def __init__(self, path: Path, dataset: str, gene_key: str | None):
        import h5py  # imported lazily: numpy/scipy alone are enough for .npy and .tsv

        self._f = h5py.File(str(path), "r")
        if dataset not in self._f:
            raise KeyError(
                f"{path}: no dataset '{dataset}'. Top-level keys: {list(self._f.keys())}"
            )
        node = self._f[dataset]
        if isinstance(node, h5py.Group):
            raise TypeError(
                f"{path}:{dataset} is a group, not a dense dataset. If it is AnnData CSR, "
                f"use _H5SparseReader (open_corpus dispatches to it automatically)."
            )
        gene_ids = _read_h5_var_ids(self._f, gene_key)
        super().__init__(node, f"{path}::{dataset}", gene_ids)

    def close(self) -> None:
        self._f.close()


class _H5SparseReader(CorpusReader):
    """AnnData CSR stored as a group of data/indices/indptr.

    Densified one chunk at a time. A corpus written this way is usually mostly zeros, so
    this costs chunk_rows x n_cols floats and nothing more.
    """

    def __init__(self, path: Path, dataset: str, gene_key: str | None):
        import h5py

        self._f = h5py.File(str(path), "r")
        g = self._f[dataset]
        self._data, self._indices, self._indptr = g["data"], g["indices"], g["indptr"]
        shape = g.attrs.get("shape", g.attrs.get("h5sparse_shape"))
        if shape is None:
            raise KeyError(f"{path}:{dataset} is CSR but carries no shape attribute")
        self.n_rows, self.n_cols = int(shape[0]), int(shape[1])
        self.gene_ids = _read_h5_var_ids(self._f, gene_key)
        self.source = f"{path}::{dataset} (CSR)"

    def _densify(self, s: int, e: int) -> np.ndarray:
        ptr = np.asarray(self._indptr[s:e + 1])
        lo, hi = int(ptr[0]), int(ptr[-1])
        cols = np.asarray(self._indices[lo:hi], dtype=np.int64)
        vals = np.asarray(self._data[lo:hi], dtype=np.float32)
        out = np.zeros((e - s, self.n_cols), dtype=np.float32)
        counts = np.diff(ptr)
        rows = np.repeat(np.arange(e - s, dtype=np.int64), counts)
        out[rows, cols] = vals
        return out

    def iter_chunks(self, chunk_rows: int) -> Iterator[tuple[int, int, np.ndarray]]:
        for s in range(0, self.n_rows, chunk_rows):
            e = min(s + chunk_rows, self.n_rows)
            yield s, e, self._densify(s, e)

    def take_rows(self, rows: Sequence[int]) -> np.ndarray:
        out = np.empty((len(rows), self.n_cols), dtype=np.float32)
        for k, r in enumerate(rows):
            out[k] = self._densify(int(r), int(r) + 1)[0]
        return out

    def close(self) -> None:
        self._f.close()


def _read_h5_var_ids(f, gene_key: str | None) -> list[str] | None:
    """Gene identifiers from an AnnData `var`, if the deposit carries any.

    Tries the caller's key first, then the conventional ones. Handles both a plain string
    array and AnnData's categorical encoding. Returns None rather than guessing when the
    deposit has no usable identifier column -- the caller then needs --corpus-genes or
    --assume-aligned, and is told so.
    """
    import h5py

    if "var" not in f:
        return None
    var = f["var"]
    keys = ([gene_key] if gene_key else []) + [
        "_index", "index", "ensg_id", "gene_ids", "gene_id", "gene_symbol", "feature_id",
    ]
    for k in keys:
        if not k or k not in var:
            continue
        node = var[k]
        if isinstance(node, h5py.Group):
            if "categories" in node and "codes" in node:
                cats = _decode(np.asarray(node["categories"]))
                codes = np.asarray(node["codes"], dtype=np.int64)
                return [cats[c] if c >= 0 else "" for c in codes]
            continue
        vals = _decode(np.asarray(node))
        # An integer row index is not a gene vocabulary; refuse it rather than align on it.
        if len(set(vals)) == len(vals) and not all(v.isdigit() for v in vals):
            return vals
    return None


class _NpyReader(_DenseArrayReader):
    def __init__(self, path: Path, gene_ids: list[str] | None):
        arr = np.load(str(path), mmap_mode="r")
        if arr.ndim != 2:
            raise ValueError(f"{path}: expected a 2-D corpus, got shape {arr.shape}")
        super().__init__(arr, str(path), gene_ids)


class _NpzReader(_DenseArrayReader):
    def __init__(self, path: Path, key: str | None, gene_ids: list[str] | None):
        z = np.load(str(path), allow_pickle=False)
        names = list(z.files)
        if key is None:
            two_d = [n for n in names if z[n].ndim == 2]
            if len(two_d) != 1:
                raise ValueError(
                    f"{path}: pass --corpus-dataset; 2-D arrays present: {two_d or names}"
                )
            key = two_d[0]
        arr = z[key]
        if os.environ.get("PRETRAINING_EXPOSURE_QUIET") != "1":
            print(
                f"  NOTE: .npz is decompressed whole ({arr.nbytes / 1e9:.2f} GB). "
                f"Convert a production corpus to .npy or .h5ad to stream it.",
                flush=True,
            )
        super().__init__(arr, f"{path}::{key}", gene_ids)


class _TextReader(CorpusReader):
    """Row-per-corpus-sample text matrix, read line by line and never wholly.

    The first column may be a row label; `--corpus-text-has-rowname` says so. A header row
    of gene identifiers is used when `--corpus-text-header` is set.
    """

    def __init__(self, path: Path, sep: str, header: bool, rowname: bool,
                 gene_ids: list[str] | None):
        self._path, self._sep = path, sep
        self._header, self._rowname = header, rowname
        self._open = (lambda: gzip.open(str(path), "rt")) if str(path).endswith(".gz") \
            else (lambda: open(str(path), "r"))
        n, ncol, hdr = 0, None, None
        with self._open() as fh:
            for i, line in enumerate(fh):
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split(sep)
                if i == 0 and header:
                    hdr = parts[1:] if rowname else parts
                    continue
                if ncol is None:
                    ncol = len(parts) - (1 if rowname else 0)
                n += 1
        self.n_rows, self.n_cols = n, int(ncol or 0)
        self.gene_ids = gene_ids if gene_ids is not None else hdr
        if self.gene_ids is not None and len(self.gene_ids) != self.n_cols:
            raise ValueError(
                f"{path}: {len(self.gene_ids)} gene ids for {self.n_cols} columns"
            )
        self.source = str(path)

    def _rows(self) -> Iterator[np.ndarray]:
        with self._open() as fh:
            for i, line in enumerate(fh):
                if not line.strip():
                    continue
                if i == 0 and self._header:
                    continue
                parts = line.rstrip("\n").split(self._sep)
                if self._rowname:
                    parts = parts[1:]
                yield np.asarray(parts, dtype=np.float32)

    def iter_chunks(self, chunk_rows: int) -> Iterator[tuple[int, int, np.ndarray]]:
        buf, s = [], 0
        for row in self._rows():
            buf.append(row)
            if len(buf) == chunk_rows:
                yield s, s + len(buf), np.stack(buf)
                s += len(buf)
                buf = []
        if buf:
            yield s, s + len(buf), np.stack(buf)

    def take_rows(self, rows: Sequence[int]) -> np.ndarray:
        want = {int(r): k for k, r in enumerate(rows)}
        out = np.empty((len(rows), self.n_cols), dtype=np.float32)
        for i, row in enumerate(self._rows()):
            if i in want:
                out[want[i]] = row
        return out


def open_corpus(path, dataset: str = "X", gene_key: str | None = None,
                gene_ids: list[str] | None = None, sep: str = "\t",
                text_header: bool = True, text_rowname: bool = True) -> CorpusReader:
    """Dispatch on suffix. `gene_ids`, when given, overrides whatever the file carries."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"corpus not found: {p}")
    suf = "".join(p.suffixes[-2:]) if p.name.endswith(".gz") else p.suffix
    if p.suffix in {".h5ad", ".h5", ".hdf5"}:
        import h5py

        with h5py.File(str(p), "r") as f:
            is_group = dataset in f and isinstance(f[dataset], h5py.Group)
        r = _H5SparseReader(p, dataset, gene_key) if is_group \
            else _H5Reader(p, dataset, gene_key)
    elif p.suffix == ".npy":
        r = _NpyReader(p, None)
    elif p.suffix == ".npz":
        r = _NpzReader(p, dataset if dataset != "X" else None, None)
    elif suf in {".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz", ".txt.gz"}:
        s = "," if ".csv" in suf else sep
        r = _TextReader(p, s, text_header, text_rowname, None)
    else:
        raise ValueError(
            f"{p}: unsupported container '{p.suffix}'. "
            f"Supported: .h5ad/.h5/.hdf5, .npy, .npz, .tsv/.csv(.gz)"
        )
    if gene_ids is not None:
        if len(gene_ids) != r.n_cols:
            raise ValueError(
                f"--corpus-genes has {len(gene_ids)} ids but the corpus has {r.n_cols} columns"
            )
        r.gene_ids = list(gene_ids)
    return r
