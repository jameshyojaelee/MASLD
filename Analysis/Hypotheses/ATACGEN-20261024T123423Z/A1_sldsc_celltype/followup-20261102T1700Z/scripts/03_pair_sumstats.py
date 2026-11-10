#!/usr/bin/env python3
"""Follow-up step 3 (Q7): restrict an AST and an ALT munged file to one identical SNP set.

Kept: SNPs with non-missing Z in both files and chi^2 <= max(80, 0.001 * max N) in both, so
ldsc removes no further SNP for either trait and the jackknife blocks coincide.
Usage: 03_pair_sumstats.py <AST_SUMSTATS> <ALT_SUMSTATS> <OUT_PREFIX>
Writes <OUT_PREFIX>_ast.sumstats.gz, <OUT_PREFIX>_alt.sumstats.gz
"""
import sys

import numpy as np
import pandas as pd


def main(ast_path, alt_path, out):
    a = pd.read_csv(ast_path, sep="\t").dropna(subset=["Z"])
    b = pd.read_csv(alt_path, sep="\t").dropna(subset=["Z"])
    snps = np.intersect1d(a["SNP"].values, b["SNP"].values)
    a, b = a[a["SNP"].isin(snps)], b[b["SNP"].isin(snps)]
    keep = set(snps)
    for d in (a, b):
        cap = max(80.0, 0.001 * d["N"].max())
        keep &= set(d.loc[d["Z"] ** 2 <= cap, "SNP"])
    for d, name in ((a, "ast"), (b, "alt")):
        d[d["SNP"].isin(keep)].to_csv(f"{out}_{name}.sumstats.gz", sep="\t", index=False,
                                      compression={"method": "gzip", "compresslevel": 1})
    print(f"{out}: {len(snps)} shared SNPs with Z, {len(keep)} kept after the chi^2 cap")


if __name__ == "__main__":
    main(*sys.argv[1:4])
