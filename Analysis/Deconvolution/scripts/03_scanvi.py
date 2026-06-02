#!/usr/bin/env python3
"""Train scANVI using CellTypist labels and export refined labels."""

from __future__ import annotations

import argparse
import os
import random
from pathlib import Path

import numpy as np
import scanpy as sc
import scvi

# --- Seed pinning (added 2026-04-22 per T0.8) -------------------------------
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
import torch
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
try:
    scvi.settings.seed = 42
except Exception:
    pass
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input .h5ad with celltypist labels")
    parser.add_argument("--output", required=True, help="Output .h5ad")
    parser.add_argument("--label-key", default="celltypist_label")
    parser.add_argument("--confidence-key", default="celltypist_confidence")
    parser.add_argument("--min-confidence", type=float, default=0.5)
    parser.add_argument("--unknown-label", default="Unknown")
    parser.add_argument("--batch-key", default="sample")
    parser.add_argument("--use-gpu", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(str(input_path))

    if args.label_key not in adata.obs:
        raise SystemExit(f"Missing label column: {args.label_key}")

    if args.confidence_key in adata.obs:
        low_conf = adata.obs[args.confidence_key].astype(float) < args.min_confidence
        # Ensure the label column is plain string for safe assignment
        adata.obs[args.label_key] = adata.obs[args.label_key].astype(str)
        adata.obs.loc[low_conf, args.label_key] = args.unknown_label

    layer = "counts" if "counts" in adata.layers else None
    scvi.model.SCVI.setup_anndata(adata, batch_key=args.batch_key, layer=layer)
    scvi_model = scvi.model.SCVI(adata)
    accelerator = "gpu" if args.use_gpu else "cpu"
    devices = 1 if args.use_gpu else 1
    scvi_model.train(accelerator=accelerator, devices=devices)

    scvi_latent = scvi_model.get_latent_representation()
    adata.obsm["X_scVI"] = scvi_latent

    scvi.model.SCANVI.setup_anndata(
        adata,
        batch_key=args.batch_key,
        labels_key=args.label_key,
        unlabeled_category=args.unknown_label,
        layer=layer,
    )
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        labels_key=args.label_key,
        unlabeled_category=args.unknown_label,
    )
    scanvi_model.train(accelerator=accelerator, devices=devices)

    preds = scanvi_model.predict(soft=True)
    adata.obs["scanvi_label"] = scanvi_model.predict()
    adata.obs["scanvi_label_prob"] = np.max(preds, axis=1)

    adata.write_h5ad(str(output_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
