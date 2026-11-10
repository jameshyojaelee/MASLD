#!/usr/bin/env python3
"""Follow-up step 6 (Q8): do the 436 COLOC r2 positive genes that rest on AST alone lean toward
endothelial or myeloid expression, compared with genes that rest on ALT alone?

Written before the comparison was run (2026-11-02). Gene support = Tier-1/2 studies with SuSiE
PP.H4 > 0.5 (the adopted COLOC r2 rule; reproduces 436 genes). Expression = liver snRNA/scRNA
pseudobulk counts per cell type, summed over samples, CPM per cell type; a gene's share in cell
type k = CPM_k / sum over the 11 cell types.
Groups: AST_only (support traits == {AST}) vs ALT_only (== {ALT}).
Primary family (BH, 2 tests, one-sided Mann-Whitney AST_only > ALT_only), unit = locus (genes
merged when TSSs are < 1 Mb apart; locus value = mean share): endothelial share; myeloid share
(Macrophages + Mono_mono_derived_cells). Positive control (one-sided ALT_only > AST_only):
hepatocyte share. Gene-level tests are descriptive.
Usage: 06_coloc_gene_celltype.py <COLOC_DIR> <TIER_TSV> <PSEUDOBULK_DIR> <GTF> <OUT_DIR>
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

DIRECT = {"NAFLD", "NASH", "PDFF"}


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def gtf_genes(path):
    rows = []
    pat = re.compile(r'gene_id "([^".]+)[^"]*";.*gene_name "([^"]+)"')
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene" or not f[0].startswith("chr") or "_" in f[0]:
                continue
            m = pat.search(f[8])
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            rows.append((m.group(1), m.group(2), f[0], tss))
    return pd.DataFrame(rows, columns=["ensembl", "symbol", "chrom", "tss"]).drop_duplicates("ensembl")


def loci(df):
    df = df.sort_values(["chrom", "tss"]).copy()
    lid, prev, k = [], (None, -10**10), 0
    for c, t in zip(df["chrom"], df["tss"]):
        if c != prev[0] or t - prev[1] >= 1_000_000:
            k += 1
        lid.append(k)
        prev = (c, t)
    df["locus"] = lid
    return df


def main(coloc, tier_tsv, pb, gtf, out):
    os.makedirs(out, exist_ok=True)
    tier = pd.read_csv(tier_tsv, sep="\t")
    tier = tier[tier["tier"].isin([1, 2])].set_index("study_name")["trait"]
    gl = pd.read_csv(os.path.join(coloc, "gene_level_coloc_tier12.csv"))
    pos = set(gl.loc[gl["susie_state_t12"] == "positive", "ensembl"])
    al = pd.read_csv(os.path.join(coloc, "susie_coloc_all_gwas.csv"), usecols=["gwas_name", "ensembl", "PP.H4.susie"])
    al = al[al["gwas_name"].isin(tier.index) & (al["PP.H4.susie"] > 0.5)]
    al["trait"] = al["gwas_name"].map(tier)
    sup = al.groupby("ensembl")["trait"].agg(lambda s: "+".join(sorted(set(s))))
    sup = sup[sup.index.isin(pos)]
    assert len(sup) == len(pos) == 436, (len(sup), len(pos))

    g = gtf_genes(gtf).set_index("ensembl")
    d = pd.DataFrame({"support": sup}).join(g, how="left")
    d["group"] = np.where(d["support"] == "AST", "AST_only", np.where(d["support"] == "ALT", "ALT_only", "other"))

    shares = {}
    for f in sorted(os.listdir(pb)):
        if f.startswith("pseudobulk_"):
            ct = f[len("pseudobulk_"):-len(".tsv.gz")]
            x = pd.read_csv(os.path.join(pb, f), sep="\t", index_col=0).sum(axis=1)
            shares[ct] = x / x.sum() * 1e6
    cpm = pd.DataFrame(shares).fillna(0)
    cpm = cpm[cpm.sum(axis=1) > 0]
    share = cpm.div(cpm.sum(axis=1), axis=0)
    share["myeloid"] = share["Macrophages"] + share["Mono_mono_derived_cells"]
    d = d.join(share[["Endothelial_cells", "myeloid", "Hepatocytes"]], on="symbol")
    d["expressed"] = d["Endothelial_cells"].notna()
    d["top_ct"] = d["symbol"].map(cpm.idxmax(axis=1))
    d.to_csv(os.path.join(out, "Q8_gene_table.tsv"), sep="\t")

    e = d[d["expressed"] & d["group"].isin(["AST_only", "ALT_only"])]
    e = loci(e)
    loc = e.groupby(["group", "locus"])[["Endothelial_cells", "myeloid", "Hepatocytes"]].mean().reset_index()
    rows = []
    for var, alt in (("Endothelial_cells", "greater"), ("myeloid", "greater"), ("Hepatocytes", "less")):
        for unit, tab in (("locus", loc), ("gene", e)):
            a = tab.loc[tab["group"] == "AST_only", var]
            b = tab.loc[tab["group"] == "ALT_only", var]
            p = mannwhitneyu(a, b, alternative=alt).pvalue
            rows.append(dict(variable=var, unit=unit, alternative="AST_only " + alt, n_ast=len(a), n_alt=len(b),
                             median_ast=a.median(), median_alt=b.median(), p=p))
    r = pd.DataFrame(rows)
    prim = (r["unit"] == "locus") & r["variable"].isin(["Endothelial_cells", "myeloid"])
    r["q_BH_primary"] = np.nan
    r.loc[prim, "q_BH_primary"] = bh(r.loc[prim, "p"])
    r.to_csv(os.path.join(out, "Q8_tests.tsv"), sep="\t", index=False)

    comp = d.groupby("support").size().sort_values(ascending=False)
    comp.to_csv(os.path.join(out, "Q8_support_combinations.tsv"), sep="\t", header=["n_genes"])
    top = pd.crosstab(d.loc[d["expressed"], "group"], d.loc[d["expressed"], "top_ct"])
    top.to_csv(os.path.join(out, "Q8_top_celltype.tsv"), sep="\t")
    print(comp.head(8).to_string())
    print(f"expressed: {int(d['expressed'].sum())} of 436; any direct-trait support: "
          f"{int(d['support'].str.split('+').map(lambda s: bool(DIRECT & set(s))).sum())}")
    print(r.to_string(index=False))
    print(top.to_string())


if __name__ == "__main__":
    main(*sys.argv[1:6])
