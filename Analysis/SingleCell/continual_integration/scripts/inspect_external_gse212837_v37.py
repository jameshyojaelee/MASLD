#!/usr/bin/env python3
"""Inspect a downloaded H5AD without materializing its expression matrix."""

from __future__ import annotations

import argparse
import json
from collections import Counter

import h5py
import numpy as np


def _text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def _attrs(node):
    return {str(key): _text(value) for key, value in node.attrs.items()}


def _matrix(node):
    if isinstance(node, h5py.Dataset):
        return {"storage": "dense", "shape": list(node.shape), "dtype": str(node.dtype), "attrs": _attrs(node)}
    result = {"storage": "group", "keys": sorted(node.keys()), "attrs": _attrs(node)}
    if "shape" in node.attrs:
        result["shape"] = [int(value) for value in node.attrs["shape"]]
    if "data" in node:
        result["data_dtype"] = str(node["data"].dtype)
        result["data_length"] = int(len(node["data"]))
        if len(node["data"]):
            sample = node["data"][: min(1_000_000, len(node["data"]))]
            result["data_sample_min"] = float(np.min(sample))
            result["data_sample_max"] = float(np.max(sample))
            result["data_sample_integer_valued"] = bool(np.equal(sample, np.floor(sample)).all())
    return result


def _column(node):
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = [_text(value) for value in node["categories"][:]]
        codes = node["codes"][:]
        counts = Counter(int(value) for value in codes if int(value) >= 0)
        return {
            "kind": "categorical", "n": int(len(codes)),
            "categories": categories,
            "counts": {categories[index]: int(counts[index]) for index in range(len(categories))},
            "missing": int(np.sum(codes < 0)),
        }
    if not isinstance(node, h5py.Dataset):
        return {"kind": "group", "keys": sorted(node.keys()), "attrs": _attrs(node)}
    values = node[:]
    result = {"kind": "dataset", "n": int(len(values)), "dtype": str(values.dtype)}
    if values.dtype.kind in {"S", "O", "U"}:
        counts = Counter(_text(value) for value in values)
        result["unique"] = len(counts)
        if len(counts) <= 100:
            result["counts"] = dict(sorted(counts.items()))
    elif len(values):
        result.update({
            "min": float(np.nanmin(values)), "max": float(np.nanmax(values)),
            "missing": int(np.sum(~np.isfinite(values))) if values.dtype.kind == "f" else 0,
        })
    return result


def inspect(path: str):
    with h5py.File(path, "r") as handle:
        result = {
            "schema_version": "masld-cl-external-h5ad-inspection-v37",
            "path": path,
            "root_keys": sorted(handle.keys()),
            "root_attrs": _attrs(handle),
            "matrices": {}, "obs": {}, "var": {},
        }
        if "X" in handle:
            result["matrices"]["X"] = _matrix(handle["X"])
        if "raw" in handle and "X" in handle["raw"]:
            result["matrices"]["raw/X"] = _matrix(handle["raw/X"])
        if "layers" in handle:
            for key in sorted(handle["layers"]):
                result["matrices"][f"layers/{key}"] = _matrix(handle["layers"][key])
        for frame in ("obs", "var"):
            if frame not in handle:
                continue
            result[frame]["attrs"] = _attrs(handle[frame])
            result[frame]["columns"] = {
                key: _column(handle[frame][key])
                for key in sorted(handle[frame])
                if key != "_index"
            }
            if "_index" in handle[frame]:
                result[frame]["index"] = _column(handle[frame]["_index"])
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    args = parser.parse_args()
    print(json.dumps(inspect(args.h5ad), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
