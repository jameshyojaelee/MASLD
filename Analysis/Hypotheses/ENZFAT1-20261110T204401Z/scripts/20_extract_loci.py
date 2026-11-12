#!/usr/bin/env python3
"""ENZFAT1 S2 step 1: raw effect rows at the 169 ALT lead variants (plus control variants) in each GWAS.

GRCh37 files (liver fat, Pan-UKBB) are matched by chr:pos and allele pair. The rsID comes from the
liver-fat file; GRCh38 files (MVP, Ghouse, FinnGen) are matched by rsID and allele pair.
MVP rows need imputation r2 >= 0.9; Pan-UKBB rows need low_confidence_EUR false.
Output: one row per variant and trait with effect allele, other allele, beta, se.

Usage: 20_extract_loci.py <project_root> <loci.tsv> <out.tsv>
"""
import sys

import numpy as np
import pandas as pd

MR = "GWAS/MR_Data"
LAB = f"{MR}/MVP/raw/labs/MVP_R4.1000G_AGR.{{}}_Mean_INT.EUR.GIA.dbGaP.txt.gz"
PHE = f"{MR}/MVP/raw/phecodes/MVP_R4.1000G_AGR.Phe_{{}}.EUR.GIA.dbGaP.txt.gz"
# control variants outside the locus table (GRCh37 positions from PRESPEC_H10.md line 120).
# PNPLA3 and GOT1 have no lead variant among the 169 loci, so they are added here as control rows
# (outside the ashr fit). GOT1 lead: smallest Pan-UKBB EUR AST p within GOT1 +/- 100 kb (gene from
# tools/magma/ref/NCBI37.3.gene.loc), MAF >= 1%, present in the liver-fat file with the same allele
# pair; this follows the H10 rule for the GPT negative control.
EXTRA = [("ctrl_HSD17B13_rs72613567", "4", 88231393), ("ctrl_HSD17B13_rs6834314", "4", 88213808),
         ("ctrl_PNPLA3_rs738409", "22", 44324727)]
GOT1 = ("10", 101156627, 101190530)
PDFF = f"{MR}/GCST90267352_PDFF_Pazoki2022.tsv.gz"


def pairs_equal(a1, a2, b1, b2):
    return ((a1 == b1) & (a2 == b2)) | ((a1 == b2) & (a2 == b1))


def scan(path, usecols, keycol_fn, wanted, **kw):
    hdr = pd.read_csv(path, sep="\t", nrows=0, **kw).columns
    names = [c.lstrip("#") for c in hdr]
    out = []
    for ch in pd.read_csv(path, sep="\t", names=names, header=0, usecols=usecols, dtype=str,
                          chunksize=2_000_000, **kw):
        ch = ch.assign(key=keycol_fn(ch))
        out.append(ch[ch["key"].isin(wanted)])
    return pd.concat(out, ignore_index=True)


def got1_lead(root):
    c, start, end = GOT1
    win = {f"{c}:{p}" for p in range(start - 100_000, end + 100_001)}
    u = scan(f"{root}/{MR}/PanUKBB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz",
             ["chr", "pos", "ref", "alt", "af_EUR", "neglog10_pval_EUR", "low_confidence_EUR"],
             lambda d: d["chr"] + ":" + d["pos"], win, compression="gzip")
    u = u[u["low_confidence_EUR"].str.lower() == "false"].dropna(subset=["af_EUR", "neglog10_pval_EUR"])
    af = u["af_EUR"].astype(float)
    u = u[(af >= 0.01) & (af <= 0.99)]
    f = scan(f"{root}/{PDFF}", ["chromosome", "base_pair_location", "effect_allele", "other_allele"],
             lambda d: d["chromosome"] + ":" + d["base_pair_location"], win)
    m = u.merge(f, on="key")
    m = m[pairs_equal(m["alt"], m["ref"], m["effect_allele"], m["other_allele"])]
    m = m[~m["key"].duplicated(keep=False)]
    lead = m.loc[m["neglog10_pval_EUR"].astype(float).idxmax()]
    print(f"GOT1 control lead: {lead['key']} {lead['ref']}/{lead['alt']} AST -log10 p {lead['neglog10_pval_EUR']} "
          f"af_EUR {lead['af_EUR']} ({len(m)} candidate variants)")
    return ("ctrl_GOT1_lead", c, int(lead["pos"]))


def main(root, loci_path, out_path):
    loci = pd.read_csv(loci_path, sep="\t", dtype=str)[["locus", "chr", "pos", "a1", "a2"]]
    extra = pd.DataFrame([dict(locus=n, chr=c, pos=str(p), a1=None, a2=None)
                          for n, c, p in EXTRA + [got1_lead(root)]])
    loci = pd.concat([loci, extra], ignore_index=True)
    loci["key"] = loci["chr"] + ":" + loci["pos"]
    wanted = set(loci["key"])
    rows = []

    def add(df, trait):
        df = df.assign(trait=trait)
        rows.append(df[["locus", "trait", "rsid", "ea", "oa", "beta", "se"]])

    def match_alleles(df):
        # control variants carry no allele pair in the locus table: accept the single row at that position
        known = df["a1"].notna()
        ok = ~known | pairs_equal(df["ea"], df["oa"], df["a1"], df["a2"])
        return df[ok]

    # liver fat (GRCh37, has rsID)
    f = scan(f"{root}/{PDFF}",
             ["variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
              "standard_error"], lambda c: c["chromosome"] + ":" + c["base_pair_location"], wanted)
    f = f.rename(columns={"variant_id": "rsid", "effect_allele": "ea", "other_allele": "oa",
                          "standard_error": "se"}).merge(loci, on="key")
    f = match_alleles(f)
    f = f[~f["locus"].duplicated(keep=False) | f["a1"].notna()]
    add(f, "pdff")
    rsmap = f.drop_duplicates("locus").set_index("rsid")["locus"]
    loci = loci.merge(f[["locus", "rsid"]].drop_duplicates("locus"), on="locus", how="left")
    pair = loci.set_index("locus")[["a1", "a2"]]
    # fill the allele pair of control variants from the liver-fat row
    for _, r in f[f["a1"].isna()].iterrows():
        pair.loc[r["locus"]] = [r["ea"], r["oa"]]
    print(f"variants with an rsID from the liver-fat file: {loci['rsid'].notna().sum()} of {len(loci)}")

    # Pan-UKBB (GRCh37, no rsID)
    for trait, code in [("alt_ukb", "30620-ALT"), ("ast_ukb", "30650-AST"), ("ggt_ukb", "30730-GGT")]:
        u = scan(f"{root}/{MR}/PanUKBB/biomarkers-{code}-both_sexes-irnt.tsv.bgz",
                 ["chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "low_confidence_EUR"],
                 lambda c: c["chr"] + ":" + c["pos"], wanted, compression="gzip")
        u = u[u["low_confidence_EUR"].str.lower() == "false"]
        u = u.rename(columns={"alt": "ea", "ref": "oa", "beta_EUR": "beta", "se_EUR": "se"}).merge(
            loci[["locus", "key", "rsid"]], on="key")
        u = u[pairs_equal(u["ea"], u["oa"], pair.loc[u["locus"], "a1"].values, pair.loc[u["locus"], "a2"].values)]
        add(u, trait)

    rs_wanted = set(rsmap.index)

    def by_rsid(df, trait):
        df = df.assign(locus=df["rsid"].map(rsmap))
        df = df[pairs_equal(df["ea"], df["oa"], pair.loc[df["locus"], "a1"].values,
                            pair.loc[df["locus"], "a2"].values)]
        add(df, trait)

    # MVP (GRCh38)
    mvp = [("alt_mvp", LAB.format("ALT"), False), ("albumin_mvp", LAB.format("Albumin"), False),
           ("platelet_mvp", LAB.format("Platelet"), False), ("nafld_mvp", PHE.format("571_5"), True),
           ("cirr_mvp", PHE.format("571_51"), True)]
    for trait, path, binary in mvp:
        cols = ["SNP_ID", "ref", "alt", "ea", "r2"] + (["or", "ci"] if binary else ["beta", "sebeta"])
        m = scan(f"{root}/{path}", cols, lambda c: c["SNP_ID"], rs_wanted)
        m = m[m["r2"].astype(float) >= 0.9]
        m = m.assign(rsid=m["SNP_ID"], oa=np.where(m["ea"] == m["alt"], m["ref"], m["alt"]))
        if binary:
            ci = m["ci"].str.split(",", expand=True).astype(float)
            m = m.assign(beta=np.log(m["or"].astype(float)),
                         se=(np.log(ci[1]) - np.log(ci[0])) / (2 * 1.959964))
        else:
            m = m.rename(columns={"sebeta": "se"})
        by_rsid(m, trait)

    # Ghouse cirrhosis and FinnGen NAFLD (GRCh38, rsID)
    g = scan(f"{root}/{MR}/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz",
             ["rsid", "effect_allele", "other_allele", "beta", "standard_error"], lambda c: c["rsid"], rs_wanted)
    by_rsid(g.rename(columns={"effect_allele": "ea", "other_allele": "oa", "standard_error": "se"}), "cirr_ghouse")
    fg = scan(f"{root}/{MR}/FinnGen/finngen_R12_NAFLD.gz", ["rsids", "ref", "alt", "beta", "sebeta"],
              lambda c: c["rsids"], rs_wanted)
    by_rsid(fg.rename(columns={"rsids": "rsid", "alt": "ea", "ref": "oa", "sebeta": "se"}), "nafld_finngen")

    out = pd.concat(rows, ignore_index=True)
    out["ea"] = out["ea"].str.upper()
    out["oa"] = out["oa"].str.upper()
    out = out[~out.duplicated(["locus", "trait"], keep=False)]
    out.to_csv(out_path, sep="\t", index=False)
    print(out.groupby("trait")["locus"].nunique().to_string())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
