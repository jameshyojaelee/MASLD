#!/usr/bin/env python3
"""H22 step 3: S-LDSC annotation files from the lifted peak sets.

Modes:
  models <OUT>
      Write <OUT>/ldsc/models.tsv: one row per regression model and the
      annotation columns it adds to baselineLD. Only annotations with
      testable == True (>= 1,000 lifted peaks) enter a model.
  l2annot <OUT> <CHR> <BFILE_PREFIX>
      ldsc/annot_l2/all.<CHR>.annot.gz: thin annot with testable columns only
      (the --l2 input).
  chrom <OUT> <CHR> <BFILE_PREFIX> <FRQ_PREFIX>
      For one chromosome: binary SNP-in-peak indicators for every annotation
      (1000G EUR bim, GRCh37, 1-based closed peak interval). Writes
        ldsc/annot_all/all.<CHR>.annot.gz             thin annot, all columns (for --l2)
        ldsc/models/<MODEL>/<MODEL>.<CHR>.annot.gz     thin annot per non-null model
        ldsc/annot_stats/<CHR>.tsv                     counts over common SNPs
      Common SNPs: 0.05 < FRQ < 0.95 in the 1000G frq file, the rule LDSC uses
      for --overlap-annot and M_5_50.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse


def load_meta(out):
    meta = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t")
    peaks = pd.read_csv(os.path.join(out, "peaks/peaks_hg19.tsv"), sep="\t")
    mem = np.fromfile(os.path.join(out, "peaks/membership.u8"), dtype=np.uint8)
    mem = mem.reshape(len(meta), len(peaks))
    return meta, peaks, mem


def write_models(out):
    meta = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t")
    ok = set(meta.loc[meta["testable"], "annotation"])
    rows = []
    lineages = sorted(meta["lineage"].dropna().unique())
    for L in lineages:
        acc = f"ACC_{L}"
        if f"SPEC_{L}" in ok and "UNION" in ok:
            rows.append(("CT_" + L, "ct", L, "", "UNION,SPEC_" + L))
        for tag, kind in (("STEA", "stea"), ("FIB", "fib"), ("FDR", "fdr")):
            cols = [c for c in (f"{tag}_GAIN_{L}", f"{tag}_LOST_{L}") if c in ok]
            if cols and acc in ok:
                rows.append((f"{tag}_{L}", kind, L, "", ",".join([acc] + cols)))
        nulls = meta[(meta["lineage"] == L) & meta["type"].isin(["null_gain", "null_lost"])]
        for k in sorted(nulls["draw"].dropna().astype(int).unique()):
            cols = [c for c in (f"NULL_GAIN_{L}_{k}", f"NULL_LOST_{L}_{k}") if c in ok]
            if cols and acc in ok:
                rows.append((f"NULL_{L}_{k}", "null", L, str(k), ",".join([acc] + cols)))
    df = pd.DataFrame(rows, columns=["model", "kind", "lineage", "draw", "columns"])
    os.makedirs(os.path.join(out, "ldsc"), exist_ok=True)
    df.to_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", index=False)
    print(df["kind"].value_counts().to_string())


def write_thin(path, names, arr):
    """Write a 0/1 uint8 matrix as a tab-separated thin annot file, quickly."""
    n, m = arr.shape
    buf = np.empty((n, 2 * m), dtype=np.uint8)
    buf[:, 0::2] = arr + ord("0")
    buf[:, 1::2] = ord("\t")
    buf[:, -1] = ord("\n")
    with gzip.open(path, "wb", compresslevel=1) as fh:
        fh.write(("\t".join(names) + "\n").encode())
        fh.write(buf.tobytes())


def snp_annot(out, c, bfile):
    meta, peaks, mem = load_meta(out)
    bim = pd.read_csv(f"{bfile}{c}.bim", sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    pos = bim["BP"].values
    assert np.all(np.diff(pos) >= 0), "bim not sorted by position"
    sel = np.where(peaks["chr"].values == c)[0]
    st, en = peaks["start"].values[sel], peaks["end"].values[sel]
    lo = np.searchsorted(pos, st, side="left")
    hi = np.searchsorted(pos, en, side="right")
    cnt = hi - lo
    rows = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)]) if cnt.sum() else np.array([], int)
    cols = np.repeat(np.arange(len(sel)), cnt)
    S = sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(bim), len(sel)))
    Mc = mem[:, sel].T.astype(np.float32)  # peaks_c x annot
    A = (S @ Mc) > 0
    return meta, bim, A.astype(np.uint8), len(sel), int(cnt.sum())


def l2annot(out, c, bfile):
    """Thin annot with testable columns only, the input to ldsc --l2.

    Untestable sets (< 1,000 lifted peaks, including their null draws) never
    enter a model; leaving them out of --l2 cuts memory and does not change the
    LD score of any other column (each column is summed independently).
    """
    c = int(c)
    meta, bim, A, _, _ = snp_annot(out, c, bfile)
    keep = np.where(meta["testable"].values)[0]
    d = os.path.join(out, "ldsc/annot_l2")
    os.makedirs(d, exist_ok=True)
    write_thin(os.path.join(d, f"all.{c}.annot.gz"), list(meta["annotation"].values[keep]), A[:, keep])
    print(f"chr{c}: {len(bim)} SNPs, {len(keep)} testable annotations")


def chrom(out, c, bfile, frq_prefix):
    c = int(c)
    meta, bim, A, n_peaks, n_pairs = snp_annot(out, c, bfile)
    names = list(meta["annotation"])

    d = os.path.join(out, "ldsc/annot_all")
    os.makedirs(d, exist_ok=True)
    write_thin(os.path.join(d, f"all.{c}.annot.gz"), names, A)

    frq = pd.read_csv(f"{frq_prefix}{c}.frq", sep=r"\s+")
    assert np.array_equal(frq["SNP"].values, bim["SNP"].values), "frq and bim SNP order differ"
    common = (frq["MAF"].values > 0.05) & (frq["MAF"].values < 0.95)
    st_df = pd.DataFrame({"annotation": names,
                          "n_snp_in": A.sum(axis=0).astype(int),
                          "n_common_in": A[common].sum(axis=0).astype(int),
                          "n_snp_total": len(bim), "n_common_total": int(common.sum())})
    d = os.path.join(out, "ldsc/annot_stats")
    os.makedirs(d, exist_ok=True)
    st_df.to_csv(os.path.join(d, f"{c}.tsv"), sep="\t", index=False)

    models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", keep_default_na=False)  # "null" is a kind, not NA
    idx = {n: i for i, n in enumerate(names)}
    for _, r in models[models["kind"] != "null"].iterrows():
        mcols = r["columns"].split(",")
        d = os.path.join(out, "ldsc/models", r["model"])
        os.makedirs(d, exist_ok=True)
        write_thin(os.path.join(d, f"{r['model']}.{c}.annot.gz"), mcols, A[:, [idx[x] for x in mcols]])
    print(f"chr{c}: {len(bim)} SNPs, {n_peaks} peaks, {n_pairs} SNP-peak pairs")


if __name__ == "__main__":
    if sys.argv[1] == "models":
        write_models(sys.argv[2])
    elif sys.argv[1] == "chrom":
        chrom(*sys.argv[2:6])
    elif sys.argv[1] == "l2annot":
        l2annot(*sys.argv[2:5])
    else:
        raise SystemExit("mode must be models, chrom or l2annot")
