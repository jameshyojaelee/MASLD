#!/usr/bin/env python3
"""Corrected working entry point for the unchanged v1.2 severity model.

Use the original CLI arguments, for example::

    python scripts/score_severity_v12_corrected.py --counts counts.npz \
        --orientation samples_x_genes --out scores.tsv

Only the default sample IDs for sample-by-gene NPZ input are corrected. Model
weights, transformations, thresholds and scoring are loaded from the original
release directory. This script does not create or adopt a new model release.
The v1.3 candidate accepts TSV only and does not expose this NPZ option.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


SOURCE_PATH = (
    Path(__file__).resolve().parents[1]
    / "release" / "masld-severity-v1.2" / "score.py"
)
_spec = importlib.util.spec_from_file_location("_severity_v12_source", SOURCE_PATH)
if _spec is None or _spec.loader is None:
    raise ImportError(f"cannot load severity scorer: {SOURCE_PATH}")
source = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(source)
_original_read_counts = source.read_counts


def read_counts(
    counts_path, genes_path=None, samples_path=None, orientation="genes_x_samples"
):
    """Read original formats, correcting only missing sample IDs after transpose."""
    if (
        not counts_path.endswith(".npz")
        or orientation != "samples_x_genes"
        or samples_path
    ):
        return _original_read_counts(counts_path, genes_path, samples_path, orientation)

    np = source.np
    with np.load(counts_path, allow_pickle=True) as bundle:
        if "samples" in bundle.files or "counts" not in bundle.files:
            return _original_read_counts(counts_path, genes_path, samples_path, orientation)
        mat = np.asarray(bundle["counts"], dtype=np.float64)
        if mat.ndim != 2:
            return _original_read_counts(counts_path, genes_path, samples_path, orientation)
        if genes_path:
            with source._open(genes_path) as handle:
                genes = np.array([line.strip() for line in handle if line.strip()])
        elif "genes" in bundle.files:
            genes = np.asarray(bundle["genes"]).astype(str)
        else:
            raise source.ScoreError("npz has no 'genes' array and --genes was not given")

    mat = mat.T
    samples = np.array([f"sample_{i+1}" for i in range(mat.shape[1])])
    if mat.shape[0] != len(genes):
        raise source.ScoreError(
            f"counts has {mat.shape[0]} rows and {len(genes)} gene ids. If your matrix is "
            f"samples-by-genes, pass --orientation samples_x_genes."
        )
    if not np.isfinite(mat).all():
        raise source.ScoreError("counts contain NaN or inf")
    if (mat < 0).any():
        raise source.ScoreError(
            "counts contain negative values; this is a counts matrix, not a log ratio"
        )
    return mat, genes, samples


def main(argv=None):
    source.read_counts = read_counts
    return source.main(argv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except source.ScoreError as error:
        print(f"score.py: {error}", file=sys.stderr)
        sys.exit(2)
