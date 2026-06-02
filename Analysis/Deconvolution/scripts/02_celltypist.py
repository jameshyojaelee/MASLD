#!/usr/bin/env python3
"""Annotate cells using CellTypist and write labels into AnnData."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc

try:
    import celltypist
except ImportError as exc:
    raise SystemExit("celltypist is not installed in this environment") from exc


def normalize_if_needed(adata, force: bool) -> None:
    if force:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)


def extract_labels(result, adata) -> pd.DataFrame:
    obs = adata.obs.copy()
    # Try to use to_adata() if available
    if hasattr(result, "to_adata"):
        annotated = result.to_adata()
        if "cell_type" in annotated.obs.columns:
            obs["celltypist_label"] = annotated.obs["cell_type"].values
        elif "majority_voting" in annotated.obs.columns:
            obs["celltypist_label"] = annotated.obs["majority_voting"].values
        elif "label" in annotated.obs.columns:
            obs["celltypist_label"] = annotated.obs["label"].values
        else:
            # Fallback to first column if present
            if annotated.obs.shape[1] > 0:
                obs["celltypist_label"] = annotated.obs.iloc[:, 0].values
    else:
        labels = getattr(result, "predicted_labels", None)
        if labels is not None:
            if isinstance(labels, pd.DataFrame):
                if "majority_voting" in labels.columns:
                    obs["celltypist_label"] = labels["majority_voting"].values
                elif "label" in labels.columns:
                    obs["celltypist_label"] = labels["label"].values
                else:
                    obs["celltypist_label"] = labels.iloc[:, 0].values
            else:
                obs["celltypist_label"] = labels
    # Add confidence if possible
    prob = getattr(result, "probability_matrix", None)
    if prob is not None:
        obs["celltypist_confidence"] = np.max(prob, axis=1)
    return obs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input .h5ad")
    parser.add_argument("--output", required=True, help="Output .h5ad")
    parser.add_argument("--model", required=True, help="CellTypist model name or path")
    parser.add_argument("--majority-voting", action="store_true")
    parser.add_argument("--normalize", action="store_true", help="Force normalize/log1p")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(str(input_path))
    normalize_if_needed(adata, args.normalize)

    # Ensure model availability if a named model is used.
    if not Path(args.model).exists():
        try:
            celltypist.models.download_model(args.model)
        except Exception:
            pass

    result = celltypist.annotate(
        adata,
        model=args.model,
        majority_voting=args.majority_voting,
    )

    obs = extract_labels(result, adata)
    if "celltypist_label" not in obs.columns:
        raise SystemExit(
            "CellTypist did not return labels. Please check the model or output format."
        )
    adata.obs = obs
    adata.write_h5ad(str(output_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
