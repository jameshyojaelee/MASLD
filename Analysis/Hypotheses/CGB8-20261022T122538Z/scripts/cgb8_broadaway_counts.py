#!/usr/bin/env python
"""Per-gene Broadaway liver eQTL meta counts (the eQTL source used by COLOC r2): min p and n variants at p<1e-5."""
import argparse
import glob

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--eqtl-dir", default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl")
ap.add_argument("--out", required=True)
a = ap.parse_args()
parts = []
for f in sorted(glob.glob(f"{a.eqtl_dir}/chr*_marginal_summary_results.tsv")):
    for ch in pd.read_csv(f, sep="\t", usecols=["ENSG", "PVAL"], chunksize=5_000_000):
        ch["sig"] = ch.PVAL < 1e-5
        parts.append(ch.groupby("ENSG").agg(min_p=("PVAL", "min"), n_p1e5=("sig", "sum")))
    print(f, flush=True)
d = pd.concat(parts).groupby(level=0).agg(min_p=("min_p", "min"), n_p1e5=("n_p1e5", "sum"))
d.index = d.index.str.split(".").str[0]
d.rename_axis("g").reset_index().to_csv(a.out, sep="\t", index=False)
print(len(d), "genes")
