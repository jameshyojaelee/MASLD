#!/usr/bin/env python3
"""Follow-up step 2: per-chromosome model files (LD scores, M, M_5_50, thin .annot.gz).

J7  = UNION + SPEC_<L> x 7 (H22 packed LD scores and annot, read-only)
RM_<T> = RMU + RM_<T> (T in E122, E100, E108, E066)
RMJ = RMU + RM_E122 + RM_E066 + UNION + SPEC_endothelial_cell + SPEC_hepatocyte
Usage: 02_models.py <F> <CHR> <H22> <BASELINE_PREFIX>
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

LINEAGES = ["b_cell", "cholangiocyte", "endothelial_cell", "fibroblast", "hepatocyte", "macrophage", "t_cell"]
SPEC = ["SPEC_" + l for l in LINEAGES]
MODELS = {"J7": ["UNION"] + SPEC}
for t in ["E122", "E100", "E108", "E066"]:
    MODELS["RM_" + t] = ["RMU", "RM_" + t]
MODELS["RMJ"] = ["RMU", "RM_E122", "RM_E066", "UNION", "SPEC_endothelial_cell", "SPEC_hepatocyte"]


def read_m(path):
    with open(path) as fh:
        return [float(x) for x in fh.readline().split()]


def main(f, c, h22, base_prefix):
    c = int(c)
    base = pd.read_csv(f"{base_prefix}{c}.l2.ldscore.gz", sep="\t", usecols=["CHR", "SNP", "BP"])
    # H22 packed LD scores: rows of l2.<c>.npy follow M.<c>.tsv
    hm = pd.read_csv(os.path.join(h22, f"ldsc/packed/M.{c}.tsv"), sep="\t")
    hsnp = pd.read_csv(os.path.join(h22, f"ldsc/packed/snps.{c}.tsv.gz"), sep="\t")
    if not np.array_equal(hsnp["SNP"].values, base["SNP"].values):
        raise SystemExit(f"chr{c}: H22 packed SNPs differ from baselineLD")
    hidx = {n: i for i, n in enumerate(hm["annotation"])}
    hcols = ["UNION"] + SPEC
    harr = np.load(os.path.join(h22, f"ldsc/packed/l2.{c}.npy"), mmap_mode="r")
    L = {n: np.asarray(harr[hidx[n]]) for n in hcols}
    M = {n: (hm.loc[hidx[n], "M"], hm.loc[hidx[n], "M_5_50"]) for n in hcols}
    hann = pd.read_csv(os.path.join(h22, f"ldsc/annot_l2/all.{c}.annot.gz"), sep=r"\s+", usecols=hcols,
                       dtype=np.float32)
    # Roadmap LD scores from step 1
    rl = pd.read_csv(os.path.join(f, f"ldsc/l2/rm.{c}.l2.ldscore.gz"), sep="\t")
    if not np.array_equal(rl["SNP"].values, base["SNP"].values):
        raise SystemExit(f"chr{c}: Roadmap LD score SNPs differ from baselineLD")
    rnames = [x[:-2] for x in rl.columns[3:]]
    rm, rm550 = read_m(os.path.join(f, f"ldsc/l2/rm.{c}.l2.M")), read_m(os.path.join(f, f"ldsc/l2/rm.{c}.l2.M_5_50"))
    for i, n in enumerate(rnames):
        L[n] = rl.iloc[:, 3 + i].to_numpy()
        M[n] = (rm[i], rm550[i])
    rann = pd.read_csv(os.path.join(f, f"ldsc/annot_l2/rm.{c}.annot.gz"), sep="\t", dtype=np.float32)
    if list(rann.columns) != rnames or len(rann) != len(hann):
        raise SystemExit(f"chr{c}: Roadmap annot columns or rows differ")
    ann = pd.concat([hann, rann], axis=1)

    for model, cols in MODELS.items():
        md = os.path.join(f, "ldsc/models", model)
        os.makedirs(md, exist_ok=True)
        pre = os.path.join(md, f"{model}.{c}")
        sub = base.copy()
        for n in cols:
            sub[n + "L2"] = L[n]
        with gzip.open(pre + ".l2.ldscore.gz", "wt", compresslevel=1) as fh:
            sub.to_csv(fh, sep="\t", index=False, float_format="%.4f")
        for suf, k in ((".l2.M", 0), (".l2.M_5_50", 1)):
            with open(pre + suf, "w") as fh:
                fh.write("\t".join(repr(float(M[n][k])) for n in cols) + "\n")
        with gzip.open(pre + ".annot.gz", "wt", compresslevel=1) as fh:
            ann[cols].to_csv(fh, sep="\t", index=False, float_format="%g")
    print(f"chr{c}: {len(base)} SNPs, {len(MODELS)} models written")


if __name__ == "__main__":
    main(*sys.argv[1:5])
