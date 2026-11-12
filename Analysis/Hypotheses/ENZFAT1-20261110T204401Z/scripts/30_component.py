#!/usr/bin/env python3
"""ENZFAT1 S3 step 1: per-SNP statistics for the ALT component not shared with liver fat (closed form).

Model (GWAS-by-subtraction form, Demange 2021): fat = lf*F1; ALT = la1*F1 + la2*F2; F1, F2 unit variance,
uncorrelated. Genetic parameters from the S1 Genomic SEM ldsc() matrix (alt_ukb, pdff):
  k = S_af / S_ff,  la2 = sqrt(S_aa - S_af^2 / S_ff).
Per SNP, on the standardized scale b = z / sqrt(N):
  gamma2 = (b_a - k * b_f) / la2
  var    = (c_aa / N_a + k^2 c_ff / N_f - 2 k I_af / sqrt(N_a N_f)) / la2^2
c_xx = max(1, LDSC intercept of trait x) (userGWAS GC = "standard"); I_af = cross-trait intercept.
Variants: GRCh37 chr:pos present in both the liver-fat file and Pan-UKBB EUR ALT (low_confidence_EUR false),
same allele pair. Effect allele = Pan-UKBB alt. Effective N = 1 / var.
Also writes chromosome 22 HapMap3 input files for the userGWAS check.
Usage: 30_component.py <project_root> <S1 out dir> <S3 out dir>
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy import stats

root, s1, out = sys.argv[1:4]
MR = f"{root}/GWAS/MR_Data"
S = pd.read_csv(f"{s1}/S_observed.tsv", sep="\t", index_col=0)
I = pd.read_csv(f"{s1}/ldsc_intercepts.tsv", sep="\t", index_col=0)
I.columns = I.index
S_aa, S_ff, S_af = S.loc["alt_ukb", "alt_ukb"], S.loc["pdff", "pdff"], S.loc["alt_ukb", "pdff"]
I_aa, I_ff, I_af = I.loc["alt_ukb", "alt_ukb"], I.loc["pdff", "pdff"], I.loc["alt_ukb", "pdff"]
k = S_af / S_ff
la2 = np.sqrt(S_aa - S_af ** 2 / S_ff)
c_aa, c_ff = max(1.0, I_aa), max(1.0, I_ff)
pars = pd.DataFrame({"param": ["S_aa", "S_ff", "S_af", "I_aa", "I_ff", "I_af", "k", "la2", "lf", "la1",
                               "share_u_fat_check"],
                     "value": [S_aa, S_ff, S_af, I_aa, I_ff, I_af, k, la2, np.sqrt(S_ff), S_af / np.sqrt(S_ff),
                               la2 ** 2 / S_aa]})
pars.to_csv(f"{out}/closed_form_parameters.tsv", sep="\t", index=False)
print(pars.to_string(index=False))

f = pd.read_csv(f"{MR}/GCST90267352_PDFF_Pazoki2022.tsv.gz", sep="\t",
                usecols=["variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele",
                         "beta", "standard_error", "n"],
                dtype={"chromosome": str, "variant_id": str})
f = f.dropna(subset=["beta", "standard_error"])
f = f[f["standard_error"] > 0]
f["key"] = f["chromosome"] + ":" + f["base_pair_location"].astype(str)
f = f[~f["key"].duplicated(keep=False)]
print(f"liver-fat rows with unique position: {len(f)}")
keys = set(f["key"])
parts = []
for ch in pd.read_csv(f"{MR}/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz", sep="\t", compression="gzip",
                      usecols=["chr", "pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR", "low_confidence_EUR"],
                      dtype={"chr": str}, chunksize=2_000_000):
    ch = ch.assign(key=ch["chr"] + ":" + ch["pos"].astype(str))
    ch = ch[ch["key"].isin(keys) & (ch["low_confidence_EUR"].astype(str).str.lower() == "false")]
    parts.append(ch.dropna(subset=["beta_EUR", "se_EUR"]))
u = pd.concat(parts, ignore_index=True)
u = u[~u["key"].duplicated(keep=False)]
m = u.merge(f, on="key")
fa, fo = m["effect_allele"].str.upper(), m["other_allele"].str.upper()
same = (fa == m["alt"]) & (fo == m["ref"])
flip = (fa == m["ref"]) & (fo == m["alt"])
m = m[same | flip].copy()
sgn = np.where(same[same | flip], 1.0, -1.0)
z_a = (m["beta_EUR"] / m["se_EUR"]).values
z_f = sgn * (m["beta"] / m["standard_error"]).values
N_a, N_f = 400822.0, m["n"].astype(float).values
b_a, b_f = z_a / np.sqrt(N_a), z_f / np.sqrt(N_f)
gam = (b_a - k * b_f) / la2
var = (c_aa / N_a + k ** 2 * c_ff / N_f - 2 * k * I_af / np.sqrt(N_a * N_f)) / la2 ** 2
se = np.sqrt(var)
z = gam / se
res = pd.DataFrame({"SNP": m["variant_id"].values, "CHR": m["chr"].values, "BP": m["pos"].values,
                    "A1": m["alt"].values, "A2": m["ref"].values, "AF_EUR_A1": m["af_EUR"].values,
                    "BETA": gam, "SE": se, "Z": z, "P": 2 * stats.norm.sf(np.abs(z)), "N": 1.0 / var,
                    "Z_ALT": z_a, "Z_FAT": z_f})
with gzip.open(f"{out}/alt_not_fat_component_GRCh37.tsv.gz", "wt") as fh:
    res.to_csv(fh, sep="\t", index=False, float_format="%.6g")
print(f"variants written: {len(res)}; median N_eff {np.median(1 / var):.0f}; "
      f"genome-wide significant (p < 5e-8): {(res['P'] < 5e-8).sum()}")
# chromosome 22 HapMap3 inputs for the userGWAS check (rsID from the liver-fat file)
hm3 = pd.read_csv(f"{root}/GWAS/ldsc_ref/hm3_merge_alleles.txt", sep="\t")
c22 = m[(m["chr"] == "22") & m["variant_id"].isin(set(hm3["SNP"]))]
pd.DataFrame({"SNP": c22["variant_id"], "A1": c22["alt"], "A2": c22["ref"], "BETA": c22["beta_EUR"],
              "SE": c22["se_EUR"], "P": 2 * stats.norm.sf(np.abs(c22["beta_EUR"] / c22["se_EUR"])), "N": N_a}
             ).to_csv(f"{out}/chr22_alt_ukb.tsv", sep="\t", index=False)
pd.DataFrame({"SNP": c22["variant_id"], "A1": c22["effect_allele"].str.upper(),
              "A2": c22["other_allele"].str.upper(), "BETA": c22["beta"], "SE": c22["standard_error"],
              "P": 2 * stats.norm.sf(np.abs(c22["beta"] / c22["standard_error"])), "N": c22["n"]}
             ).to_csv(f"{out}/chr22_pdff.tsv", sep="\t", index=False)
ld = pd.read_csv(f"{root}/GWAS/ldsc_ref/eur_w_ld_chr/22.l2.ldscore.gz", sep="\t")
ref = ld[["SNP", "CHR", "BP", "MAF"]].merge(hm3, on="SNP")
ref[["SNP", "CHR", "BP", "MAF", "A1", "A2"]].to_csv(f"{out}/chr22_ref_eur_w_ld_hm3.tsv", sep="\t", index=False)
print(f"chr22 HapMap3 variants for userGWAS: {len(c22)}; reference rows {len(ref)}")
