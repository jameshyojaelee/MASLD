#!/usr/bin/env python3
"""ENZFAT1 S0: one raw GWAS file -> LDSC munge input (SNP A1 A2 BETA P N), HapMap3 SNPs only.

A1 is the effect allele. BETA is the signed effect (ln OR for the MVP phecodes).
Also writes the raw rows of the orientation-check SNPs (all SNPs, not only HapMap3).
Filters: MVP imputation r2 >= 0.9; Pan-UKBB low_confidence_EUR false.

Usage: 00_build_input.py <trait> <project_root> <out_dir>
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy import stats

MR = "GWAS/MR_Data"
TRAITS = {
    "pdff": dict(kind="generic", file=f"{MR}/GCST90267352_PDFF_Pazoki2022.tsv.gz", snp="variant_id",
                 ea="effect_allele", oa="other_allele", beta="beta", se="standard_error", p="p_value", n="n"),
    "cirr_ghouse": dict(kind="generic", file=f"{MR}/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz",
                        snp="rsid", ea="effect_allele", oa="other_allele", beta="beta", se="standard_error",
                        p="p_value", n="samplesize"),
    "nafld_finngen": dict(kind="generic", file=f"{MR}/FinnGen/finngen_R12_NAFLD.gz", snp="rsids", ea="alt",
                          oa="ref", beta="beta", se="sebeta", p="pval", n=438857),
    "alt_mvp": dict(kind="mvp", file=f"{MR}/MVP/raw/labs/MVP_R4.1000G_AGR.ALT_Mean_INT.EUR.GIA.dbGaP.txt.gz",
                    binary=False),
    "cirr_mvp": dict(kind="mvp", file=f"{MR}/MVP/raw/phecodes/MVP_R4.1000G_AGR.Phe_571_51.EUR.GIA.dbGaP.txt.gz",
                     binary=True),
    "nafld_mvp": dict(kind="mvp", file=f"{MR}/MVP/raw/phecodes/MVP_R4.1000G_AGR.Phe_571_5.EUR.GIA.dbGaP.txt.gz",
                      binary=True),
    "alt_ukb": dict(kind="panukbb", file=f"{MR}/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz", n=400822),
    "ast_ukb": dict(kind="panukbb", file=f"{MR}/PanUKBB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz", n=399482),
    "ggt_ukb": dict(kind="panukbb", file=f"{MR}/PanUKBB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz", n=400751),
}
# rsID -> GRCh37 chr:pos (positions read from the liver-fat file, which is GRCh37)
CHECK = {"rs738409": "22:44324727", "rs58542926": "19:19379549", "rs1260326": "2:27730940",
         "rs641738": "19:54676763"}
CHECK_POS = {v: k for k, v in CHECK.items()}
OUT_COLS = ["SNP", "A1", "A2", "BETA", "P", "N"]
CHK_COLS = ["SNP", "A1", "A2", "BETA", "Z", "P", "N"]


def z_from_p(beta, p):
    p = np.clip(p.astype(float), 1e-300, 1.0)
    return np.sign(beta) * stats.norm.isf(p / 2.0)


def hm3_reference(root):
    ma = pd.read_csv(f"{root}/GWAS/ldsc_ref/hm3_merge_alleles.txt", sep="\t", dtype=str)
    ma = ma[~ma["SNP"].duplicated(keep=False)]
    mp = pd.read_csv(f"{root}/GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv", sep="\t", dtype=str)
    mp.columns = ["key", "SNP"]
    mp = mp.drop_duplicates()
    mp = mp[~mp["key"].duplicated(keep=False) & ~mp["SNP"].duplicated(keep=False)]
    return set(ma["SNP"]), mp.merge(ma, on="SNP").set_index("key")["SNP"]


def chunks_generic(cfg, path):
    use = [cfg[k] for k in ("snp", "ea", "oa", "beta", "se", "p")]
    if isinstance(cfg["n"], str):
        use.append(cfg["n"])
    hdr = pd.read_csv(path, sep="\t", nrows=0).columns
    names = [c.lstrip("#") for c in hdr]
    for ch in pd.read_csv(path, sep="\t", names=names, header=0, usecols=use, dtype={cfg["snp"]: str},
                          chunksize=2_000_000):
        ch = ch.dropna(subset=[cfg["snp"], cfg["beta"], cfg["se"], cfg["p"]])
        ch = ch[ch[cfg["se"]].astype(float) > 0]
        n = ch[cfg["n"]].values if isinstance(cfg["n"], str) else float(cfg["n"])
        beta = ch[cfg["beta"]].astype(float)
        yield len(ch), 0, pd.DataFrame({
            "SNP": ch[cfg["snp"]].values, "A1": ch[cfg["ea"]].str.upper().values,
            "A2": ch[cfg["oa"]].str.upper().values, "BETA": beta.values,
            "Z": (beta / ch[cfg["se"]].astype(float)).values,
            "P": np.clip(ch[cfg["p"]].astype(float).values, 1e-300, 1.0), "N": n})


def chunks_mvp(cfg, path):
    eff = "or" if cfg["binary"] else "beta"
    use = ["SNP_ID", "ref", "alt", "ea", "num_samples", eff, "pval", "r2"]
    if not cfg["binary"]:
        use.append("sebeta")
    for ch in pd.read_csv(path, sep="\t", usecols=use, dtype={"SNP_ID": str}, chunksize=2_000_000):
        n_in = len(ch)
        ch = ch.dropna(subset=["SNP_ID", eff, "pval", "r2"])
        low = int((ch["r2"].astype(float) < 0.9).sum())
        ch = ch[ch["r2"].astype(float) >= 0.9]
        ea = ch["ea"].str.upper()
        ref, alt = ch["ref"].str.upper(), ch["alt"].str.upper()
        ch = ch[(ea == ref) | (ea == alt)]
        ea, ref, alt = ea.loc[ch.index], ref.loc[ch.index], alt.loc[ch.index]
        oa = np.where(ea == alt, ref, alt)
        if cfg["binary"]:
            ch = ch.assign(_b=np.log(ch["or"].astype(float)))
            ok = np.isfinite(ch["_b"].values)
            ch, ea, oa = ch[ok], ea[ok], oa[ok]
            beta = ch["_b"]
            z = z_from_p(beta.values, ch["pval"].values)
        else:
            beta = ch["beta"].astype(float)
            z = (beta / ch["sebeta"].astype(float)).values
        # MVP prints p = 0 below 1e-321; rebuild P from Z where available
        p = ch["pval"].astype(float).values
        p = np.where(p > 0, p, 2.0 * stats.norm.sf(np.abs(z)))
        yield n_in, low, pd.DataFrame({"SNP": ch["SNP_ID"].values, "A1": ea.values, "A2": oa,
                                       "BETA": beta.values, "Z": z, "P": np.clip(p, 1e-300, 1.0),
                                       "N": ch["num_samples"].values})


def chunks_panukbb(cfg, path, key2snp):
    use = ["chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"]
    wanted = set(key2snp.index) | set(CHECK_POS)
    for ch in pd.read_csv(path, sep="\t", usecols=use, dtype={"chr": str}, compression="gzip",
                          chunksize=2_000_000):
        n_in = len(ch)
        ch = ch.assign(key=ch["chr"] + ":" + ch["pos"].astype(str))
        ch = ch[ch["key"].isin(wanted)]
        low = int((ch["low_confidence_EUR"].astype(str).str.lower() != "false").sum())
        ch = ch[ch["low_confidence_EUR"].astype(str).str.lower() == "false"]
        ch = ch.dropna(subset=["beta_EUR", "se_EUR", "neglog10_pval_EUR"])
        snp = ch["key"].map(key2snp)
        snp = snp.where(snp.notna(), ch["key"].map(CHECK_POS))
        beta = ch["beta_EUR"].astype(float)
        p = np.power(10.0, -ch["neglog10_pval_EUR"].astype(float).values)
        yield n_in, low, pd.DataFrame({"SNP": snp.values, "A1": ch["alt"].values, "A2": ch["ref"].values,
                                       "BETA": beta.values, "Z": (beta / ch["se_EUR"].astype(float)).values,
                                       "P": np.clip(p, 1e-300, 1.0), "N": float(cfg["n"])})


def main(trait, root, out_dir):
    cfg = TRAITS[trait]
    path = f"{root}/{cfg['file']}"
    hm3, key2snp = hm3_reference(root)
    if cfg["kind"] == "generic":
        gen = chunks_generic(cfg, path)
    elif cfg["kind"] == "mvp":
        gen = chunks_mvp(cfg, path)
    else:
        gen = chunks_panukbb(cfg, path, key2snp)
    kept, checks, read, dropped = [], [], 0, 0
    for n_in, low, df in gen:
        read += n_in
        dropped += low
        checks.append(df[df["SNP"].isin(CHECK)][CHK_COLS])
        kept.append(df[df["SNP"].isin(hm3)][OUT_COLS])
    df = pd.concat(kept, ignore_index=True)
    n_dup = int(df["SNP"].duplicated(keep=False).sum())
    df = df.drop_duplicates(subset="SNP", keep=False)
    with gzip.open(f"{out_dir}/inputs/{trait}.tsv.gz", "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)
    chk = pd.concat(checks, ignore_index=True)
    chk.insert(0, "trait", trait)
    chk.to_csv(f"{out_dir}/checks/{trait}.tsv", sep="\t", index=False)
    print(f"{trait}: read {read} rows; dropped by quality filter {dropped}; wrote {len(df)} HapMap3 SNPs; "
          f"dropped {n_dup} rows with duplicated rsID; N median {df['N'].median():.0f}; "
          f"check SNPs found {sorted(chk['SNP'])}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
