#!/usr/bin/env python3
"""Combine CellTypist and scANVI labels into a consensus annotation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input .h5ad with labels")
    parser.add_argument("--output", required=True, help="Output .h5ad")
    parser.add_argument("--celltypist-key", default="celltypist_label")
    parser.add_argument("--celltypist-confidence-key", default="celltypist_confidence")
    parser.add_argument("--scanvi-key", default="scanvi_label")
    parser.add_argument("--scanvi-prob-key", default="scanvi_label_prob")
    parser.add_argument("--min-scanvi-prob", type=float, default=0.6)
    parser.add_argument("--min-celltypist-prob", type=float, default=0.5)
    parser.add_argument("--unknown-label", default="Unknown")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(str(input_path))

    ct = adata.obs.get(args.celltypist_key)
    sv = adata.obs.get(args.scanvi_key)

    ct_prob = adata.obs.get(args.celltypist_confidence_key)
    sv_prob = adata.obs.get(args.scanvi_prob_key)

    consensus = []
    source = []
    for i in range(adata.n_obs):
        sv_ok = sv is not None and sv_prob is not None and float(sv_prob.iloc[i]) >= args.min_scanvi_prob
        ct_ok = ct is not None and ct_prob is not None and float(ct_prob.iloc[i]) >= args.min_celltypist_prob

        if sv_ok:
            consensus.append(str(sv.iloc[i]))
            source.append("scanvi")
        elif ct_ok:
            consensus.append(str(ct.iloc[i]))
            source.append("celltypist")
        else:
            if sv is not None:
                consensus.append(str(sv.iloc[i]))
                source.append("scanvi_lowconf")
            elif ct is not None:
                consensus.append(str(ct.iloc[i]))
                source.append("celltypist_lowconf")
            else:
                consensus.append(args.unknown_label)
                source.append("unknown")

    adata.obs["consensus_label"] = consensus
    adata.obs["consensus_source"] = source

    adata.write_h5ad(str(output_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
