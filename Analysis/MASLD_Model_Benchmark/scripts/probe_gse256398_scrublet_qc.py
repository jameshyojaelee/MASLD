#!/usr/bin/env python3
"""Probe one frozen GSE256398 donor through author-aligned Scrublet and QC."""

from __future__ import annotations

import argparse
import importlib.metadata
import inspect
import json
from pathlib import Path
import time

import h5py
import numpy as np
from scipy import sparse
import scrublet


class GSE256398QCProbeError(RuntimeError):
    """Raised when the installed QC runtime or source matrix differs."""


def decode_many(values: np.ndarray) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def read_10x_counts(path: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][...])
        counts = sparse.csc_matrix(
            (
                np.asarray(matrix["data"][...]),
                np.asarray(matrix["indices"][...]),
                np.asarray(matrix["indptr"][...]),
            ),
            shape=shape,
        ).T.tocsr()
        barcodes = decode_many(matrix["barcodes"][...])
        names = decode_many(matrix["features/name"][...])
    if counts.shape != (len(barcodes), len(names)) or counts.shape[1] != 36_601:
        raise GSE256398QCProbeError("10x count matrix axis differs")
    if counts.nnz == 0 or np.any(counts.data <= 0) or len(set(barcodes)) != len(barcodes):
        raise GSE256398QCProbeError("10x count values or barcodes differ")
    return counts, barcodes, names


def basic_qc(counts: sparse.csr_matrix, names: list[str]) -> dict[str, np.ndarray]:
    total = np.asarray(counts.sum(axis=1)).ravel()
    genes = np.diff(counts.indptr)
    mitochondrial = np.asarray([name.startswith("MT-") for name in names], dtype=bool)
    mt_counts = np.asarray(counts[:, mitochondrial].sum(axis=1)).ravel()
    percent_mt = np.divide(
        mt_counts * 100.0,
        total,
        out=np.zeros_like(mt_counts, dtype=np.float64),
        where=total > 0,
    )
    retained = (
        (genes > 200)
        & (genes < 6_500)
        & (total < 40_000)
        & (percent_mt < 20.0)
    )
    return {
        "total_counts": total,
        "detected_genes": genes,
        "percent_mt": percent_mt,
        "basic_qc_retained": retained,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if args.output.exists():
        raise GSE256398QCProbeError("probe output exists")
    started = time.monotonic()
    counts, barcodes, names = read_10x_counts(args.h5)
    qc = basic_qc(counts, names)
    model = scrublet.Scrublet(counts, random_state=args.seed)
    scores, predicted = model.scrub_doublets()
    scores = np.asarray(scores)
    predicted = np.asarray(predicted)
    if (
        scores.shape != (counts.shape[0],)
        or predicted.shape != (counts.shape[0],)
        or not np.all(np.isfinite(scores))
        or predicted.dtype != np.bool_
        or model.threshold_ is None
    ):
        raise GSE256398QCProbeError("Scrublet result differs")
    final = qc["basic_qc_retained"] & ~predicted
    receipt = {
        "schema_version": "masld-bench-gse256398-scrublet-qc-probe-v1",
        "status": "pass",
        "source_h5": args.h5.name,
        "nuclei_input": counts.shape[0],
        "features_input": counts.shape[1],
        "nonzero_entries": counts.nnz,
        "basic_qc_retained": int(qc["basic_qc_retained"].sum()),
        "scrublet_predicted_doublets": int(predicted.sum()),
        "final_retained": int(final.sum()),
        "scrublet_threshold": float(model.threshold_),
        "scrublet_score_min": float(scores.min()),
        "scrublet_score_max": float(scores.max()),
        "seed": args.seed,
        "scrublet_constructor_signature": str(inspect.signature(scrublet.Scrublet)),
        "scrublet_method_signature": str(inspect.signature(scrublet.Scrublet.scrub_doublets)),
        "package_versions": {
            name: importlib.metadata.version(name)
            for name in ("h5py", "numpy", "scipy", "scrublet")
        },
        "elapsed_seconds": time.monotonic() - started,
        "author_order": "Scrublet_per_donor_then_Seurat_threshold_equivalent_basic_QC",
        "outcomes_or_disease_labels_read": False,
        "model_training_activated": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=False)
    args.output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
