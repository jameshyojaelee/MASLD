#!/usr/bin/env python3
"""A1 step 3: SNP-level S-LDSC annotations for one chromosome (adapted from H22).

Binary annotation: 1 if the SNP (1000G EUR bim, GRCh37) lies in any lifted peak
of the set (1-based closed interval). Continuous annotation (ACCESS, GC, RESP,
SIGN, NULLRESP): the standardized value of the tested peak containing the SNP,
mean over peaks if more than one, 0 outside tested peaks.

Writes, for testable annotations only,
  ldsc/annot_l2/all.<CHR>.annot.gz           thin annot, every testable column (--l2 input)
  ldsc/models/<M>/<M>.<CHR>.annot.gz         thin annot per non-null model (--overlap-annot)
  ldsc/annot_stats/<CHR>.tsv                 per column: n SNPs > 0, and over common
                                             SNPs (0.05 < MAF < 0.95): n, sum, sum of squares
Usage: 03_make_annot.py <OUT> <CHR> <BFILE_PREFIX> <FRQ_PREFIX>
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

NULL_KINDS = ("nullcont", "nulltert")


def load(out):
    pk = pd.read_csv(os.path.join(out, "peaks/peaks.tsv"), sep="\t")
    for c in ("lifted", "only", "prom", "tested", "gain", "lost", "gain10", "lost10"):
        pk[c] = pk[c].astype(str).str.upper().eq("TRUE")
    ann = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t", keep_default_na=False)
    ann["testable"] = ann["testable"].astype(str).str.upper().eq("TRUE")
    ann = ann[ann["testable"]].reset_index(drop=True)
    nulls = {}
    for L in ("hepatocyte", "stellate"):
        f = os.path.join(out, f"peaks/nulls_{L}.npz")
        if os.path.exists(f):
            nulls[L] = np.load(f)
    return pk, ann, nulls


def peak_values(pk, ann, nulls, sel):
    """Per-peak value matrix for the peaks in sel (len(sel) x n_annot)."""
    n = len(pk)
    V = np.zeros((len(sel), len(ann)), np.float32)
    lin = pk["lineage"].values
    for j, r in ann.iterrows():
        s, L = r["set"], r["lineage"]
        isL = lin == L
        if s == "UNION":
            v = np.ones(n, bool)
        elif s == "ALL":
            v = isL
        elif s == "ONLY":
            v = isL & pk["only"].values
        elif s == "PROM":
            v = isL & pk["only"].values & pk["prom"].values
        elif s == "DIST":
            v = isL & pk["only"].values & ~pk["prom"].values
        elif s == "TESTED":
            v = isL & pk["tested"].values
        elif s in ("GAIN", "LOST", "GAIN10", "LOST10"):
            v = isL & pk[s.lower()].values
        elif s in ("ACCESS", "GC", "RESP", "SIGN"):
            col = {"ACCESS": "access_std", "GC": "gc_std", "RESP": "resp_std", "SIGN": "sign_std"}[s]
            v = np.where(isL & pk["tested"].values, pk[col].values, 0.0)
        elif s in ("NULLRESP", "NULLGAIN", "NULLLOST"):
            z = nulls[L]
            k = int(r["draw"]) - 1
            v = np.zeros(n, np.float32)
            if s == "NULLRESP":
                v[z["idx"]] = z["resp"][k]
            else:
                v[z["idx"]] = (z["lab"][k] == (1 if s == "NULLGAIN" else -1)).astype(np.float32)
        else:
            raise ValueError(s)
        V[:, j] = np.nan_to_num(np.asarray(v, np.float32)[sel])
    return V


def fmt_column(col, cont):
    n = len(col)
    if not cont:
        return (col > 0).astype(np.uint8).reshape(n, 1) + ord("0")
    u, inv = np.unique(np.round(col.astype(np.float64), 5), return_inverse=True)
    strs = ["%.5f" % x for x in u]
    w = max(len(s) for s in strs)
    tab = np.frombuffer("".join(s.rjust(w) for s in strs).encode(), np.uint8).reshape(len(u), w)
    return tab[inv]


def write_thin(path, names, A, cont):
    n = A.shape[0]
    parts = []
    for j in range(A.shape[1]):
        parts.append(fmt_column(A[:, j], cont[j]))
        sep = np.full((n, 1), ord("\t") if j < A.shape[1] - 1 else ord("\n"), np.uint8)
        parts.append(sep)
    buf = np.hstack(parts)
    with gzip.open(path, "wb", compresslevel=1) as fh:
        fh.write(("\t".join(names) + "\n").encode())
        fh.write(buf.tobytes())


def main(out, c, bfile, frq_prefix):
    c = int(c)
    pk, ann, nulls = load(out)
    names = list(ann["annotation"])
    cont = list(ann["kind"] == "continuous")
    bim = pd.read_csv(f"{bfile}{c}.bim", sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    pos = bim["BP"].values
    assert np.all(np.diff(pos) >= 0), "bim not sorted by position"

    sel = np.where(pk["lifted"].values & (pk["chr19"].values == c))[0]
    V = peak_values(pk, ann, nulls, sel)
    st, en = pk["s19"].values[sel].astype(np.int64), pk["e19"].values[sel].astype(np.int64)
    lo = np.searchsorted(pos, st, side="left")
    hi = np.searchsorted(pos, en, side="right")
    cnt = hi - lo
    rows = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)]) if cnt.sum() else np.array([], int)
    cols = np.repeat(np.arange(len(sel)), cnt)
    S = sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(bim), len(sel)))

    A = np.zeros((len(bim), len(names)), np.float32)
    bj = [j for j in range(len(names)) if not cont[j]]
    cj = [j for j in range(len(names)) if cont[j]]
    if bj:
        A[:, bj] = (S @ (V[:, bj] > 0).astype(np.float32)) > 0
    if cj:
        # continuous: mean over the tested peaks of that lineage containing the SNP
        n_multi = 0
        for L in sorted(set(ann.loc[cj, "lineage"])):
            jj = [j for j in cj if ann.loc[j, "lineage"] == L]
            tmask = ((pk["lineage"].values[sel] == L) & pk["tested"].values[sel]).astype(np.float32)
            den = np.asarray(S @ tmask).ravel()
            num = S @ (V[:, jj] * tmask[:, None])
            A[:, jj] = np.where(den[:, None] > 0, num / np.maximum(den, 1)[:, None], 0.0)
            n_multi += int((den > 1).sum())
        print(f"chr{c}: SNPs in > 1 tested peak of a lineage: {n_multi}")

    d = os.path.join(out, "ldsc/annot_l2")
    os.makedirs(d, exist_ok=True)
    write_thin(os.path.join(d, f"all.{c}.annot.gz"), names, A, cont)

    frq = pd.read_csv(f"{frq_prefix}{c}.frq", sep=r"\s+")
    assert np.array_equal(frq["SNP"].values, bim["SNP"].values), "frq and bim SNP order differ"
    common = (frq["MAF"].values > 0.05) & (frq["MAF"].values < 0.95)
    Ac = A[common].astype(np.float64)
    d = os.path.join(out, "ldsc/annot_stats")
    os.makedirs(d, exist_ok=True)
    pd.DataFrame({"annotation": names, "kind": ann["kind"].values,
                  "n_snp_nonzero": (A != 0).sum(axis=0).astype(int),
                  "n_common_nonzero": (Ac != 0).sum(axis=0).astype(int),
                  "sum_common": Ac.sum(axis=0), "sumsq_common": (Ac ** 2).sum(axis=0),
                  "n_snp_total": len(bim), "n_common_total": int(common.sum())}).to_csv(
        os.path.join(d, f"{c}.tsv"), sep="\t", index=False)

    models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", dtype=str, keep_default_na=False)
    idx = {n: i for i, n in enumerate(names)}
    for _, r in models[~models["kind"].isin(NULL_KINDS)].iterrows():
        mcols = r["columns"].split(",")
        md = os.path.join(out, "ldsc/models", r["model"])
        os.makedirs(md, exist_ok=True)
        ii = [idx[x] for x in mcols]
        write_thin(os.path.join(md, f"{r['model']}.{c}.annot.gz"), mcols, A[:, ii], [cont[i] for i in ii])
    print(f"chr{c}: {len(bim)} SNPs, {len(sel)} lifted peaks, {int(cnt.sum())} SNP-peak pairs, {len(names)} columns")


if __name__ == "__main__":
    main(*sys.argv[1:5])
