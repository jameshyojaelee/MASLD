#!/usr/bin/env python3
"""JOIN7 step 1: per gene and trait, MAGMA gene p and GWAS variants within gene +/- 50 kb.

Traits: pdff (Pazoki 2022, GRCh37), cirr (Ghouse 2024, GRCh38), alt_ukb (Pan-UKBB EUR, GRCh37),
alt_mvp (MVP EUR, GRCh38). GRCh38 windows come from Ghouse positions of the Pazoki rsIDs inside the
GRCh37 window (PRESPEC_JOIN7.md section 4). Also writes rs738409 rows aligned to G.

Usage: 10_extract.py <project_root> <out_dir>
"""
import sys

import numpy as np
import pandas as pd

MR = "GWAS/MR_Data"
PDFF = f"{MR}/GCST90267352_PDFF_Pazoki2022.tsv.gz"
CIRR = f"{MR}/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz"
UKB = f"{MR}/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz"
MVP = f"{MR}/MVP/raw/labs/MVP_R4.1000G_AGR.ALT_Mean_INT.EUR.GIA.dbGaP.txt.gz"
GENES = "tools/pops/example/data/utils/gene_annot_jun10.txt"
MAGMA = "Analysis/Hypotheses/GSET8-20261110T235125Z/out_23188024/{}.genes.out"
FLANK = 50_000
GWS = 5e-8
CHUNK = 2_000_000

TEST = ("THRB GLP1R GIPR GCGR FGFR1 FGFR2 FGFR3 PPARA PPARD PPARG NR1H4 ACACA ACACB DGAT2 SCD FASN "
        "MAP3K5 CASP1 CASP2 CASP3 CASP4 CASP5 CASP6 CASP7 CASP8 CASP9 CASP10 CASP12 CASP14 CCR2 CCR5 "
        "LGALS3 LOXL2 SLC5A2 MPC1 MPC2 KHK FGF1 FGFR4 KLB ADORA3 AOC3 PNPLA3").split()
CONTEXT = "TM6SF2 GCKR HSD17B13 CIDEB MTARC1 INHBE GPAM MOGAT2 NR3C1 RORA CNR1".split()
RS = "rs738409"
RS37 = ("22", 44324727)


def windows(root):
    g = pd.read_csv(f"{root}/{GENES}", sep="\t", dtype={"CHR": str})
    g = g[g["NAME"].isin(TEST + CONTEXT)].drop_duplicates("NAME")
    g = g.assign(lo=(g["START"] - FLANK).clip(lower=0), hi=g["END"] + FLANK)
    return g.rename(columns={"NAME": "gene", "ENSGID": "ensgid", "CHR": "chr"})


def scan(path, usecols, chrcol, poscol, win, extra=None, **kw):
    """Rows of `path` inside any window (chr, lo, hi) of `win`; returns long frame with gene."""
    out, keep_extra = [], []
    hdr = pd.read_csv(path, sep="\t", nrows=0, **kw).columns
    names = [c.lstrip("#") for c in hdr]
    for ch in pd.read_csv(path, sep="\t", names=names, header=0, usecols=usecols,
                          dtype={chrcol: str}, chunksize=CHUNK, **kw):
        if extra is not None:
            keep_extra.append(extra(ch))
        for c, w in win.groupby("chr"):
            sub = ch[ch[chrcol] == c]
            if sub.empty:
                continue
            pos = sub[poscol].astype(np.int64).values
            for r in w.itertuples():
                m = (pos >= r.lo) & (pos <= r.hi)
                if m.any():
                    out.append(sub[m].assign(gene=r.gene))
    res = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=usecols + ["gene"])
    ex = pd.concat(keep_extra, ignore_index=True) if extra is not None else None
    return res, ex


def summarise(df, trait, pcol):
    rows = []
    for gene, d in df.groupby("gene"):
        p = d[pcol].astype(float)
        i = p.idxmin()
        rows.append(dict(gene=gene, trait=trait, n_var=len(d), n_gws=int((p < GWS).sum()),
                         min_p=float(p.min()), lead=d.loc[i, "lead_id"]))
    return pd.DataFrame(rows)


def align_g(ea, oa, beta):
    if ea == "G":
        return beta
    if oa == "G":
        return -beta
    return np.nan


def main(root, out):
    win = windows(root)
    win.to_csv(f"{out}/gene_windows_grch37.tsv", sep="\t", index=False)
    missing = sorted(set(TEST + CONTEXT) - set(win["gene"]))
    print("genes absent from gene file:", missing)
    ctrl = []

    # liver fat (GRCh37)
    pd_cols = ["variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele",
               "other_allele", "beta"]
    pdf, ex = scan(f"{root}/{PDFF}", pd_cols, "chromosome", "base_pair_location", win,
                   extra=lambda c: c[c["variant_id"] == RS])
    r = ex.iloc[0]
    ctrl.append(dict(file="pdff", ea=r["effect_allele"], oa=r["other_allele"], beta=float(r["beta"]),
                     beta_G=align_g(r["effect_allele"].upper(), r["other_allele"].upper(), float(r["beta"]))))
    pdf["lead_id"] = pdf["variant_id"]
    summ = [summarise(pdf, "pdff", "p_value")]
    pdf[["gene", "variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele",
         "beta", "p_value"]].to_csv(f"{out}/variants_pdff.tsv.gz", sep="\t", index=False)

    # Pan-UKBB ALT (GRCh37)
    u_cols = ["chr", "pos", "ref", "alt", "beta_EUR", "neglog10_pval_EUR", "low_confidence_EUR"]
    udf, ex = scan(f"{root}/{UKB}", u_cols, "chr", "pos", win, compression="gzip",
                   extra=lambda c: c[(c["chr"] == RS37[0]) & (c["pos"] == RS37[1])])
    udf = udf[udf["low_confidence_EUR"].astype(str).str.lower() == "false"].dropna(
        subset=["beta_EUR", "neglog10_pval_EUR"])
    udf = udf.assign(p=10 ** (-udf["neglog10_pval_EUR"].astype(float)),
                     lead_id=udf["chr"] + ":" + udf["pos"].astype(str) + ":" + udf["ref"] + ":" + udf["alt"])
    summ.append(summarise(udf, "alt_ukb", "p"))
    ex = ex[ex["low_confidence_EUR"].astype(str).str.lower() == "false"]
    for r in ex.itertuples():
        ctrl.append(dict(file="alt_ukb", ea=r.alt, oa=r.ref, beta=float(r.beta_EUR),
                         beta_G=align_g(r.alt.upper(), r.ref.upper(), float(r.beta_EUR))))

    # GRCh38 windows: Ghouse positions of Pazoki rsIDs inside each GRCh37 window
    rs2gene = pdf[["variant_id", "gene"]].drop_duplicates()
    rsset = set(rs2gene["variant_id"])
    g_cols = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "p_value", "rsid"]
    hits = []
    for ch in pd.read_csv(f"{root}/{CIRR}", sep="\t", usecols=g_cols, dtype={"chromosome": str, "rsid": str},
                          chunksize=CHUNK):
        hits.append(ch[ch["rsid"].isin(rsset | {RS})])
    hits = pd.concat(hits, ignore_index=True)
    m = hits.merge(rs2gene, left_on="rsid", right_on="variant_id")
    w37 = win.set_index("gene")["chr"]
    m = m[m["chromosome"] == m["gene"].map(w37)]
    win38 = (m.groupby("gene").agg(chr=("chromosome", "first"), lo=("base_pair_location", "min"),
                                   hi=("base_pair_location", "max"), n_map=("rsid", "size")).reset_index())
    win38.to_csv(f"{out}/gene_windows_grch38_mapped.tsv", sep="\t", index=False)
    r = hits[hits["rsid"] == RS].iloc[0]
    ctrl.append(dict(file="cirr", ea=r["effect_allele"], oa=r["other_allele"], beta=float(r["beta"]),
                     beta_G=align_g(r["effect_allele"].upper(), r["other_allele"].upper(), float(r["beta"]))))

    cdf, _ = scan(f"{root}/{CIRR}", g_cols, "chromosome", "base_pair_location", win38)
    cdf["lead_id"] = cdf["rsid"]
    summ.append(summarise(cdf, "cirr", "p_value"))
    cdf.to_csv(f"{out}/variants_cirr.tsv.gz", sep="\t", index=False)

    # MVP ALT (GRCh38)
    v_cols = ["SNP_ID", "chrom", "pos", "ref", "alt", "ea", "beta", "pval", "r2"]
    vdf, ex = scan(f"{root}/{MVP}", v_cols, "chrom", "pos", win38,
                   extra=lambda c: c[c["SNP_ID"] == RS])
    vdf = vdf[vdf["r2"].astype(float) >= 0.9].dropna(subset=["beta", "pval"])
    vdf = vdf[(vdf["ea"] == vdf["ref"]) | (vdf["ea"] == vdf["alt"])]
    vdf["lead_id"] = vdf["SNP_ID"]
    summ.append(summarise(vdf, "alt_mvp", "pval"))
    for r in ex.itertuples():
        oa = r.ref if r.ea == r.alt else r.alt
        ctrl.append(dict(file="alt_mvp", ea=r.ea, oa=oa, beta=float(r.beta),
                         beta_G=align_g(str(r.ea).upper(), str(oa).upper(), float(r.beta))))

    summ = pd.concat(summ, ignore_index=True)

    # MAGMA gene p (gene body)
    mags = []
    for t, f in (("pdff", "pdff"), ("cirr", "cirrhosis")):
        mg = pd.read_csv(f"{root}/{MAGMA.format(f)}", sep=r"\s+")
        mg = mg.merge(win[["gene", "ensgid"]], left_on="GENE", right_on="ensgid")
        mg = mg.rename(columns={"P": "magma_p", "ZSTAT": "magma_z", "NSNPS": "magma_nsnps"})
        mags.append(mg[["gene", "magma_p", "magma_z", "magma_nsnps"]].assign(trait=t))
    summ = summ.merge(pd.concat(mags, ignore_index=True), on=["gene", "trait"], how="outer")
    summ.to_csv(f"{out}/gene_trait_raw.tsv", sep="\t", index=False)

    ctrl = pd.DataFrame(ctrl)
    ctrl["pass"] = ctrl["beta_G"] > 0
    ctrl.to_csv(f"{out}/control_rs738409.tsv", sep="\t", index=False)
    print(ctrl.to_string(index=False))
    if not (ctrl["pass"].all() and set(ctrl["file"]) == {"pdff", "cirr", "alt_ukb", "alt_mvp"}):
        sys.exit("rs738409 alignment check failed")


if __name__ == "__main__":
    np.random.seed(20261111)
    main(sys.argv[1], sys.argv[2])
