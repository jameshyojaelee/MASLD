#!/usr/bin/env python3
"""Write label summary tables for annotated references."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def summarize(path: Path, label_key: str) -> pd.DataFrame:
    adata = sc.read_h5ad(str(path))
    counts = adata.obs[label_key].astype(str).value_counts().rename_axis("label").reset_index(name="count")
    counts["fraction"] = counts["count"] / counts["count"].sum()
    counts.insert(0, "species", path.stem.split("_")[1])
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mouse", required=True, help="Mouse consensus h5ad")
    parser.add_argument("--human", required=True, help="Human consensus h5ad")
    parser.add_argument("--label-key", default="consensus_label")
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    mouse_df = summarize(Path(args.mouse), args.label_key)
    human_df = summarize(Path(args.human), args.label_key)
    out = pd.concat([mouse_df, human_df], ignore_index=True)
    out.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
