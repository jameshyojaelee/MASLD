#!/usr/bin/env python3
"""A4 step 5b: S-LDSC annotations and LD-score packing for the zonation peak sets.

Adapted copy of H22 03_make_annot.py (snp_annot, write_thin, annot_stats) and
04_pack_ldscores.py; H22 files are not edited. Peak sets come from
work/ldsc/sets.tsv.gz (hg19, 1-based closed intervals as in H22) and
work/ldsc/annotations.tsv (testable flag).

Modes:
  annot <A4> <CHR> <BFILE_PREFIX> <FRQ_PREFIX>
      ldsc/annot_l2/all.<CHR>.annot.gz  thin annot, testable columns (--l2 input)
      ldsc/annot_stats/<CHR>.tsv         SNP counts over common SNPs (0.05 < MAF < 0.95)
  pack <A4> <CHR> <BASELINE_PREFIX>
      checks the LD score SNP column against baselineLD, stores float32 .npy, M files,
      and per-model LD scores and thin .annot.gz (read by --overlap-annot) for MAIN
      (ALL_HEP, PC, PP) and SEC (ALL_HEP, PCq, PPq).
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

MODELS = {"MAIN": ["ALL_HEP", "PC", "PP"], "SEC": ["ALL_HEP", "PCq", "PPq"]}


def ld(a4, *p):
    return os.path.join(a4, "work", "ldsc", *p)


def testable(a4):
    meta = pd.read_csv(ld(a4, "annotations.tsv"), sep="\t")
    return list(meta.loc[meta["testable"], "annotation"])


def write_thin(path, names, arr):
    n, m = arr.shape
    buf = np.empty((n, 2 * m), dtype=np.uint8)
    buf[:, 0::2] = arr + ord("0")
    buf[:, 1::2] = ord("\t")
    buf[:, -1] = ord("\n")
    with gzip.open(path, "wb", compresslevel=1) as fh:
        fh.write(("\t".join(names) + "\n").encode())
        fh.write(buf.tobytes())


def annot(a4, c, bfile, frq_prefix):
    c = int(c)
    names = testable(a4)
    sets = pd.read_csv(ld(a4, "sets.tsv.gz"), sep="\t")
    sets = sets[sets["chr"] == c]
    bim = pd.read_csv(f"{bfile}{c}.bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    pos = bim["BP"].values
    assert np.all(np.diff(pos) >= 0), "bim not sorted by position"
    lo = np.searchsorted(pos, sets["start"].values, side="left")
    hi = np.searchsorted(pos, sets["end"].values, side="right")
    cnt = hi - lo
    rows = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)]) if cnt.sum() else np.array([], int)
    cols = np.repeat(np.arange(len(sets)), cnt)
    S = sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(bim), len(sets)))
    A = ((S @ sets[names].to_numpy(np.float32)) > 0).astype(np.uint8)
    os.makedirs(ld(a4, "annot_l2"), exist_ok=True)
    write_thin(ld(a4, "annot_l2", f"all.{c}.annot.gz"), names, A)
    frq = pd.read_csv(f"{frq_prefix}{c}.frq", sep=r"\s+")
    assert np.array_equal(frq["SNP"].values, bim["SNP"].values), "frq and bim SNP order differ"
    common = (frq["MAF"].values > 0.05) & (frq["MAF"].values < 0.95)
    os.makedirs(ld(a4, "annot_stats"), exist_ok=True)
    pd.DataFrame({"annotation": names, "n_snp_in": A.sum(axis=0).astype(int),
                  "n_common_in": A[common].sum(axis=0).astype(int),
                  "n_snp_total": len(bim), "n_common_total": int(common.sum())}).to_csv(
        ld(a4, "annot_stats", f"{c}.tsv"), sep="\t", index=False)
    print(f"chr{c}: {len(bim)} SNPs, {len(sets)} peaks, {int(cnt.sum())} SNP-peak pairs, {len(names)} annotations")


def read_m(path):
    with open(path) as fh:
        return [float(x) for x in fh.readline().split()]


def pack(a4, c, base_prefix):
    c = int(c)
    l2 = pd.read_csv(ld(a4, "l2", f"all.{c}.l2.ldscore.gz"), sep="\t")
    base = pd.read_csv(f"{base_prefix}{c}.l2.ldscore.gz", sep="\t", usecols=["CHR", "SNP", "BP"])
    if not np.array_equal(l2["SNP"].values, base["SNP"].values):
        raise SystemExit(f"chr{c}: LD score SNP column differs from baselineLD")
    names = [x[:-2] if x.endswith("L2") else x for x in l2.columns[3:]]
    if names != testable(a4):
        raise SystemExit(f"chr{c}: LD score columns differ from the testable annotations")
    os.makedirs(ld(a4, "packed"), exist_ok=True)
    np.save(ld(a4, "packed", f"l2.{c}.npy"), np.ascontiguousarray(l2.iloc[:, 3:].to_numpy(dtype=np.float32).T))
    l2[["CHR", "SNP", "BP"]].to_csv(ld(a4, "packed", f"snps.{c}.tsv.gz"), sep="\t", index=False)
    m = read_m(ld(a4, "l2", f"all.{c}.l2.M"))
    m550 = read_m(ld(a4, "l2", f"all.{c}.l2.M_5_50"))
    pd.DataFrame({"annotation": names, "M": m, "M_5_50": m550}).to_csv(ld(a4, "packed", f"M.{c}.tsv"), sep="\t", index=False)
    A = pd.read_csv(ld(a4, "annot_l2", f"all.{c}.annot.gz"), sep=r"\s+", dtype=np.float32)
    if list(A.columns) != names:
        raise SystemExit(f"chr{c}: annot columns differ from the LD score columns")
    idx = {n: i for i, n in enumerate(names)}
    for model, cols in MODELS.items():
        if not all(x in idx for x in cols):
            continue
        ii = [idx[x] for x in cols]
        md = ld(a4, "models", model)
        os.makedirs(md, exist_ok=True)
        pre = os.path.join(md, f"{model}.{c}")
        sub = l2.iloc[:, [0, 1, 2] + [3 + i for i in ii]].copy()
        sub.columns = ["CHR", "SNP", "BP"] + [x + "L2" for x in cols]
        with gzip.open(pre + ".l2.ldscore.gz", "wt", compresslevel=1) as fh:
            sub.to_csv(fh, sep="\t", index=False, float_format="%.3f")
        for suf, vec in ((".l2.M", m), (".l2.M_5_50", m550)):
            with open(pre + suf, "w") as fh:
                fh.write("\t".join(repr(vec[i]) for i in ii) + "\n")
        with gzip.open(pre + ".annot.gz", "wt", compresslevel=1) as fh:
            A[cols].to_csv(fh, sep="\t", index=False, float_format="%g")
    print(f"chr{c}: {len(l2)} SNPs x {len(names)} annotations packed")


if __name__ == "__main__":
    if sys.argv[1] == "annot":
        annot(*sys.argv[2:6])
    elif sys.argv[1] == "pack":
        pack(*sys.argv[2:5])
    else:
        raise SystemExit("mode must be annot or pack")
