#!/usr/bin/env python3
"""Subsample ScaleSC annotated atlas for deconvolution reference.

MuSiC/BayesPrism become very slow with >100K cells. This script
subsamples proportionally per cell type to ~100K cells total.
"""
import argparse
import logging
import numpy as np
import anndata as ad

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Annotated h5ad with cell_type + raw counts")
    parser.add_argument("--output", required=True, help="Output subsampled h5ad")
    parser.add_argument("--max-cells", type=int, default=100000)
    parser.add_argument("--min-per-type", type=int, default=50,
                        help="Minimum cells per cell type (rare types kept fully)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    np.random.seed(args.seed)

    log.info(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input)
    log.info(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Use raw counts for deconvolution reference
    if adata.raw is not None:
        log.info("Using adata.raw for raw counts")
        adata_raw = adata.raw.to_adata()
        adata_raw.obs = adata.obs.copy()
    else:
        log.info("No adata.raw; using X directly")
        adata_raw = adata

    ct_counts = adata_raw.obs["cell_type"].value_counts()
    total = ct_counts.sum()
    log.info(f"Cell types: {len(ct_counts)}")

    # Proportional sampling with minimum per type
    target = args.max_cells
    indices = []
    for ct, n in ct_counts.items():
        prop_target = max(args.min_per_type, int(target * n / total))
        n_sample = min(n, prop_target)
        ct_indices = np.where(adata_raw.obs["cell_type"] == ct)[0]
        sampled = np.random.choice(ct_indices, size=n_sample, replace=False)
        indices.extend(sampled)
        log.info(f"  {ct}: {n:,} -> {n_sample:,}")

    indices = sorted(indices)
    adata_sub = adata_raw[indices].copy()
    log.info(f"Subsampled: {adata_sub.n_obs:,} cells x {adata_sub.n_vars:,} genes")

    adata_sub.write_h5ad(args.output)
    log.info(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
