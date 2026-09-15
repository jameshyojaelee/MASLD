#!/usr/bin/env python3
"""C1 endpoint 2, step 4a: pick the Borzoi reproduction set.

The enlarged Borzoi column mixes two runs: the archived src/79 panel (1,236 Currin leads) and
this package's run (everything else). A mixed-provenance column needs a reproduction check, so
200 leads that the archived panel already scored are rescored here with this package's loader and
the two values are compared. Deterministic pick: sort by canon, take every k-th row.
"""
from __future__ import annotations

import argparse
import os

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DIRFEAT = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc/direction_features_caqtl.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=200)
    args = ap.parse_args()

    d = pd.read_csv(DIRFEAT, sep="\t")
    d = d[d["borzoi_atac_delta_labelframe"].notna()].copy()
    d["chr"] = d["canon"].str.split(":").str[0]
    d["pos_hg38"] = d["canon"].str.split(":").str[1].astype(int)
    d["ref"] = d["canon"].str.split(":").str[2]
    d["alt"] = d["canon"].str.split(":").str[3]
    d = d.sort_values("canon", kind="mergesort").reset_index(drop=True)
    k = max(1, len(d) // args.n)
    pick = d.iloc[::k].head(args.n).copy()
    pick["priority"] = 0
    pick["target"] = "borzoi_reproduction"
    pick[["chr", "pos_hg38", "ref", "alt", "canon", "priority", "target"]].to_csv(
        args.out, sep="\t", index=False)
    print(f"archived Borzoi leads {len(d)}; stride {k}; reproduction targets {len(pick)}")


if __name__ == "__main__":
    main()
