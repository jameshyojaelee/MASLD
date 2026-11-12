#!/usr/bin/env python3
"""ENZFAT1 F2 step 1 (Addendum 2): one F2 GWAS file -> LDSC munge input (SNP A1 A2 BETA P N), HapMap3 only.

Also writes orientation-check rows (all variants, not only HapMap3) and a build check.
A1 is the effect allele. Binary trait (type 2 diabetes): BETA = ln OR.
Usage: 40_build_f2.py <trait> <project_root> <run_dir> <out_dir>
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy import stats

trait, root, run, out = sys.argv[1:5]
DL = f"{run}/downloads_untrusted"
T2D = f"{root}/Analysis/Hypotheses/DIR6-20261111T024808Z/downloads_untrusted/GCST90475667/GCST90475667.h.tsv.gz"
CFG = {
    "bmi": dict(file=f"{root}/GWAS/MR_Data/Yengo_2018_BMI_GWAS.txt.gz", snp="SNP", ea="Tested_Allele",
                oa="Other_Allele", beta="BETA", se="SE", p="P", n="N", chr="CHR", pos="POS"),
    "tg": dict(file=f"{DL}/tg_GCST90239664/GCST90239664_buildGRCh37.tsv", snp="variant_id", ea="effect_allele",
               oa="other_allele", beta="beta", se="standard_error", p="p_value", n="n", chr="chromosome",
               pos="base_pair_location"),
    "fi": dict(file=f"{DL}/fi_GCST90002238/GCST90002238_buildGRCh37.tsv.gz", snp=None, ea="effect_allele",
               oa="other_allele", beta="beta", se="standard_error", p="p_value", n="sample_size", chr="chromosome",
               pos="base_pair_location"),
    "dpw": dict(file=f"{DL}/dpw_GCST007461/DrinksPerWeek.txt.gz", snp="RSID", ea="ALT", oa="REF", beta="BETA",
                se="SE", p="PVALUE", n="N", chr="CHROM", pos="POS"),
    "height": dict(file=f"{DL}/height_GCST90245992/GCST90245992_buildGRCh37.tsv", snp="variant_id",
                   ea="effect_allele", oa="other_allele", beta="beta", se="standard_error", p="p_value", n="n",
                   chr="chromosome", pos="base_pair_location"),
    "t2d": dict(file=T2D, snp="rsid", ea="effect_allele", oa="other_allele", beta="beta", se="standard_error",
                p="p_value", n=432648, chr=None, pos=None, binary=True),
}
CHECK = {"rs1558902": "16:53803574", "rs964184": "11:116648917", "rs1260326": "2:27730940",
         "rs7903146": "10:114758349", "rs1229984": "4:100239319", "rs2943641": "2:227093745",
         "rs1042725": "12:66358347"}
CHECK_POS = {v: k for k, v in CHECK.items()}
cfg = CFG[trait]
path = cfg["file"]
hdr = list(pd.read_csv(path, sep="\t", nrows=0).columns)
ma = pd.read_csv(f"{root}/GWAS/ldsc_ref/hm3_merge_alleles.txt", sep="\t", dtype=str)
hm3 = set(ma["SNP"])
mp = pd.read_csv(f"{root}/GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv", sep="\t", dtype=str)
mp.columns = ["key", "SNP"]
mp = mp[~mp["key"].duplicated(keep=False) & ~mp["SNP"].duplicated(keep=False)]
key2snp = mp.set_index("key")["SNP"]
ld = pd.concat([pd.read_csv(f"{root}/GWAS/ldsc_ref/eur_w_ld_chr/{c}.l2.ldscore.gz", sep="\t",
                            usecols=["SNP", "CHR", "BP"], dtype=str) for c in range(1, 23)])
ld_pos = (ld["CHR"] + ":" + ld["BP"]).set_axis(ld["SNP"])

if cfg.get("binary"):
    snpcol = cfg["snp"] if cfg["snp"] in hdr else "hm_rsid"
    effcol = "beta" if "beta" in hdr else "odds_ratio"
    ncol = "n" if "n" in hdr else None
    cfg.update(snp=snpcol, beta=effcol)
    if ncol:
        cfg["n"] = ncol
    print(f"t2d columns used: snp={snpcol} effect={effcol} n={'column' if ncol else cfg['n']}")
use = [c for c in (cfg["snp"], cfg["ea"], cfg["oa"], cfg["beta"], cfg["se"], cfg["p"], cfg["chr"], cfg["pos"])
       if c]
if isinstance(cfg["n"], str):
    use.append(cfg["n"])
kept, checks, pos_rows, read = [], [], [], 0
qc = dict(se_na=0, se_finite=0, p_zero=0, p_lt_1e300=0, p_na=0)
for ch in pd.read_csv(path, sep="\t", usecols=use, dtype=str, chunksize=2_000_000):
    read += len(ch)
    _se = pd.to_numeric(ch[cfg["se"]], errors="coerce"); _p = pd.to_numeric(ch[cfg["p"]], errors="coerce")
    qc["se_na"] += int(_se.isna().sum()); qc["se_finite"] += int(np.isfinite(_se).sum())
    qc["p_na"] += int(_p.isna().sum()); qc["p_zero"] += int((_p == 0).sum()); qc["p_lt_1e300"] += int((_p < 1e-300).sum())
    # Binary trait (MVP T2D): standard_error is NA in the file, so SE is not required (Z comes from P below).
    req = [cfg["beta"], cfg["p"]] + ([] if cfg.get("binary") else [cfg["se"]])
    ch = ch.dropna(subset=req)
    for c in (cfg["chr"], cfg["pos"]):
        if c:
            ch[c] = ch[c].str.strip()  # GLGC TG positions are space-padded
    key = (ch[cfg["chr"]].str.replace("chr", "") + ":" + ch[cfg["pos"]]) if cfg["chr"] else None
    if cfg["snp"] is None:
        snp = key.map(key2snp)
        snp = snp.where(snp.notna(), key.map(CHECK_POS))
    else:
        snp = ch[cfg["snp"]]
    beta = ch[cfg["beta"]].astype(float)
    if cfg.get("binary") and cfg["beta"] == "odds_ratio":
        beta = np.log(beta)
    se = ch[cfg["se"]].astype(float)
    pval = np.clip(ch[cfg["p"]].astype(float).values, 1e-300, 1.0)
    z = (beta / se).values
    se_ok = np.isfinite(se.values) & (se.values > 0)
    if cfg.get("binary"):
        z = np.where(se_ok, z, np.sign(beta.values) * stats.norm.isf(pval / 2))
        se_ok = np.ones(len(z), dtype=bool)
    n = ch[cfg["n"]].astype(float).values if isinstance(cfg["n"], str) else float(cfg["n"])
    df = pd.DataFrame({"SNP": snp.values, "A1": ch[cfg["ea"]].str.upper().values,
                       "A2": ch[cfg["oa"]].str.upper().values, "BETA": beta.values, "Z": z,
                       "P": pval, "N": n, "KEY": key.values if key is not None else None})
    df = df[np.isfinite(df["BETA"]) & np.isfinite(df["Z"]) & se_ok]
    checks.append(df[df["SNP"].isin(CHECK)].drop(columns="KEY"))
    sel = df["SNP"].isin(hm3)
    kept.append(df[sel][["SNP", "A1", "A2", "BETA", "P", "N"]])
    if cfg["chr"] and cfg["snp"] is not None:
        pos_rows.append(pd.DataFrame({"SNP": df.loc[sel, "SNP"].values, "key": df.loc[sel, "KEY"].values}))
d = pd.concat(kept, ignore_index=True)
n_dup = int(d["SNP"].duplicated(keep=False).sum())
d = d.drop_duplicates("SNP", keep=False)
with gzip.open(f"{out}/inputs/{trait}.tsv.gz", "wt") as fh:
    d.to_csv(fh, sep="\t", index=False)
chk = pd.concat(checks, ignore_index=True)
chk.insert(0, "trait", trait)
chk.to_csv(f"{out}/checks/{trait}.tsv", sep="\t", index=False)
if pos_rows:
    pr = pd.concat(pos_rows).drop_duplicates("SNP", keep=False)
    pr = pr[pr["SNP"].isin(ld_pos.index)]
    frac = float((pr["key"].values == ld_pos.loc[pr["SNP"]].values).mean())
    build = f"GRCh37 position match {frac:.4f} of {len(pr)} HapMap3 rsIDs"
else:
    frac = np.nan if cfg["snp"] is not None else len(d)
    build = ("matched by GRCh37 chr:pos to " + str(len(d)) + " HapMap3 SNPs") if cfg["snp"] is None \
        else "no position column used (rsID only)"
pd.DataFrame([dict(trait=trait, build_check=build, value=frac)]).to_csv(f"{out}/checks/{trait}_build.tsv",
                                                                         sep="\t", index=False)
print(f"{trait} raw-row QC: {qc}")
print(f"{trait}: read {read} rows; wrote {len(d)} HapMap3 SNPs; dropped {n_dup} duplicated-rsID rows; "
      f"N median {d['N'].median():.0f}; {build}; check SNPs {sorted(chk['SNP'].unique())}")
