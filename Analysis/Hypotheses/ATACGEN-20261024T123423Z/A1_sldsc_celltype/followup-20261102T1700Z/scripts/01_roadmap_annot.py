#!/usr/bin/env python3
"""Follow-up step 1: Roadmap H3K27ac annotations for one chromosome.

Writes ldsc/annot_l2/rm.<CHR>.annot.gz (thin; RMU, RM_E122, RM_E100, RM_E108, RM_E066; rows =
1000G EUR bim) and ldsc/annot_stats/rm.<CHR>.tsv (common-SNP counts, 0.05 < FRQ < 0.95).
A SNP is in a peak if start < pos <= end (0-based BED, hg19).
Usage: 01_roadmap_annot.py <F> <CHR> <BFILE_PREFIX> <FRQ_PREFIX>
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

TISSUES = ["E122", "E100", "E108", "E066"]


def merged(df):
    df = df.sort_values("s")
    s, e = df["s"].to_numpy(), df["e"].to_numpy()
    out_s, out_e = [], []
    for a, b in zip(s, e):
        if out_e and a <= out_e[-1]:
            out_e[-1] = max(out_e[-1], b)
        else:
            out_s.append(a)
            out_e.append(b)
    return np.array(out_s), np.array(out_e)


def in_peaks(pos, s, e):
    i = np.searchsorted(s, pos, side="left") - 1  # last interval with s < pos
    ok = i >= 0
    res = np.zeros(len(pos), bool)
    res[ok] = pos[ok] <= e[i[ok]]
    return res


def main(f, c, bfile, frq_prefix):
    bim = pd.read_csv(f"{bfile}{c}.bim", sep=r"\s+", header=None, usecols=[1, 3], names=["SNP", "BP"])
    frq = pd.read_csv(f"{frq_prefix}{c}.frq", sep=r"\s+", usecols=["SNP", "MAF"])
    if not np.array_equal(bim["SNP"].values, frq["SNP"].values):
        raise SystemExit(f"chr{c}: bim and frq SNP order differ")
    common = ((frq["MAF"] > 0.05) & (frq["MAF"] < 0.95)).to_numpy()
    pos = bim["BP"].to_numpy()
    cols = {}
    for t in TISSUES:
        pk = pd.read_csv(os.path.join(f, "external/roadmap", f"{t}-H3K27ac.narrowPeak.gz"), sep="\t",
                         header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
        s, e = merged(pk[pk["chr"] == f"chr{c}"])
        cols["RM_" + t] = in_peaks(pos, s, e)
    A = pd.DataFrame({"RMU": np.any(np.column_stack(list(cols.values())), axis=1), **cols}).astype(np.int8)
    os.makedirs(os.path.join(f, "ldsc/annot_l2"), exist_ok=True)
    os.makedirs(os.path.join(f, "ldsc/annot_stats"), exist_ok=True)
    with gzip.open(os.path.join(f, "ldsc/annot_l2", f"rm.{c}.annot.gz"), "wt", compresslevel=1) as fh:
        A.to_csv(fh, sep="\t", index=False)
    pd.DataFrame({"annotation": A.columns, "n_snp_in": A.sum().to_numpy(),
                  "n_common_in": A[common].sum().to_numpy(), "n_snp_total": len(A),
                  "n_common_total": int(common.sum())}).to_csv(
        os.path.join(f, "ldsc/annot_stats", f"rm.{c}.tsv"), sep="\t", index=False)
    print(f"chr{c}: {len(A)} SNPs; " + ", ".join(f"{k} {int(v)}" for k, v in A.sum().items()))


if __name__ == "__main__":
    main(*sys.argv[1:5])
