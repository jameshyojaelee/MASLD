#!/usr/bin/env python3
"""Compute-node fixture for the SnapATAC2 merge boundary schema."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import polars as pl
import snapatac2 as snap


PACKAGE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "atac_context_v3_finalize", PACKAGE / "02c_finalize_consensus_peaks.py"
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("failed to load consensus finalizer")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def native(start: int, pvalue: float) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": ["chr1"],
            "start0": [start],
            "end": [start + 100],
            "name": ["."],
            "score": [100],
            "strand": ["."],
            "signal_value": [10.0],
            "p_value": [pvalue],
            "q_value": [pvalue - 1.0],
            "peak": [50],
        }
    )


left = native(1000, 20.0)
right = native(1100, 10.0)
raw = snap.tl.merge_peaks(
    {
        "left": MODULE.merge_input_frame(left),
        "right": MODULE.merge_input_frame(right),
    },
    snap.genome.hg38,
    half_width=250,
)
assert raw.columns == ["Peaks", "right", "left"] or raw.columns == ["Peaks", "left", "right"]
assert raw.height == 1, raw
merged = MODULE.normalize_merged_peak_frame(raw)
assert merged.height == 1, merged
assert merged[0, "chrom"] == "chr1", merged
assert merged[0, "start0"] == 800, merged
assert merged[0, "end"] == 1300, merged
assert merged[0, "end"] - merged[0, "start0"] == 500, merged
duplicated_fixture = MODULE.normalize_merged_peak_frame(pl.concat([raw, raw]))
assert duplicated_fixture.height == 2, duplicated_fixture
print("SnapATAC2 merge-boundary fixture passed")
