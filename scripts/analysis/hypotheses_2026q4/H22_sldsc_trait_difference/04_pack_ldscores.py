#!/usr/bin/env python3
"""H22 step 4: after `ldsc.py --l2` on the testable annotations of one chromosome,
(i) check that the LD score SNP column equals baselineLD's, (ii) store the
matrix as float32 .npy for on-the-fly null models, (iii) write the per-model
LD score and M files for the non-null models.

Usage: 04_pack_ldscores.py <OUT> <CHR> <BASELINE_PREFIX>
  BASELINE_PREFIX e.g. GWAS/ldsc_ref/sldsc/baselineLD_v2.2/baselineLD.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd


def read_m(path):
    with open(path) as fh:
        return [float(x) for x in fh.readline().split()]


def main(out, c, base_prefix):
    c = int(c)
    l2 = pd.read_csv(os.path.join(out, f"ldsc/l2_testable/all.{c}.l2.ldscore.gz"), sep="\t")
    base = pd.read_csv(f"{base_prefix}{c}.l2.ldscore.gz", sep="\t", usecols=["CHR", "SNP", "BP"])
    if not np.array_equal(l2["SNP"].values, base["SNP"].values):
        raise SystemExit(f"chr{c}: LD score SNP column differs from baselineLD")
    names = [x[:-2] if x.endswith("L2") else x for x in l2.columns[3:]]
    meta = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t")
    if names != list(meta.loc[meta["testable"], "annotation"]):
        raise SystemExit(f"chr{c}: LD score columns differ from the testable rows of annotations.tsv")
    d = os.path.join(out, "ldsc/packed")
    os.makedirs(d, exist_ok=True)
    # annotation-major so one annotation is a contiguous read
    np.save(os.path.join(d, f"l2.{c}.npy"), np.ascontiguousarray(l2.iloc[:, 3:].to_numpy(dtype=np.float32).T))
    l2[["CHR", "SNP", "BP"]].to_csv(os.path.join(d, f"snps.{c}.tsv.gz"), sep="\t", index=False)
    m = read_m(os.path.join(out, f"ldsc/l2_testable/all.{c}.l2.M"))
    m550 = read_m(os.path.join(out, f"ldsc/l2_testable/all.{c}.l2.M_5_50"))
    pd.DataFrame({"annotation": names, "M": m, "M_5_50": m550}).to_csv(
        os.path.join(d, f"M.{c}.tsv"), sep="\t", index=False)

    idx = {n: i for i, n in enumerate(names)}
    models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", keep_default_na=False)  # "null" is a kind, not NA
    for _, r in models[models["kind"] != "null"].iterrows():
        cols = r["columns"].split(",")
        ii = [idx[x] for x in cols]
        md = os.path.join(out, "ldsc/models", r["model"])
        os.makedirs(md, exist_ok=True)
        pre = os.path.join(md, f"{r['model']}.{c}")
        sub = l2.iloc[:, [0, 1, 2] + [3 + i for i in ii]].copy()
        sub.columns = ["CHR", "SNP", "BP"] + [x + "L2" for x in cols]
        with gzip.open(pre + ".l2.ldscore.gz", "wt", compresslevel=1) as fh:
            sub.to_csv(fh, sep="\t", index=False, float_format="%.3f")
        for suf, vec in ((".l2.M", m), (".l2.M_5_50", m550)):
            with open(pre + suf, "w") as fh:
                fh.write("\t".join(repr(vec[i]) for i in ii) + "\n")
    print(f"chr{c}: {len(l2)} SNPs x {len(names)} annotations packed")


if __name__ == "__main__":
    main(*sys.argv[1:4])
