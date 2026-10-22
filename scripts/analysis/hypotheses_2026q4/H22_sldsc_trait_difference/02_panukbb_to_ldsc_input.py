#!/usr/bin/env python3
"""H22 step 2: Pan-UKBB per-phenotype file -> LDSC munge input, EUR columns only.

Pan-UKBB flat files are GRCh37 and carry no rsID. chr:pos is mapped to rsID with
the project's HapMap3 GRCh37 map; a row is kept only if its {ref, alt} alleles
equal the HapMap3 merge-allele pair for that rsID (removes multi-allelic
duplicates before munging). beta_EUR is the effect of the alt allele.
Rows with low_confidence_EUR true or missing EUR statistics are dropped.
P = 10^-neglog10_pval_EUR, floored at 1e-300 (such SNPs exceed the S-LDSC
chi-square cap anyway).

Usage: 02_panukbb_to_ldsc_input.py <panukbb.tsv.bgz> <hm3_chrbp_to_rsid.tsv>
       <hm3_merge_alleles.txt> <N_EUR> <out.tsv.gz>
"""
import gzip
import sys

import numpy as np
import pandas as pd


def main(src, map_fh, ma_fh, n_eur, out):
    mp = pd.read_csv(map_fh, sep="\t", dtype=str)
    mp.columns = ["key", "SNP"]
    mp = mp.drop_duplicates()  # the map file repeats 78,310 identical rows
    mp = mp[~mp["key"].duplicated(keep=False) & ~mp["SNP"].duplicated(keep=False)]
    ma = pd.read_csv(ma_fh, sep="\t", dtype=str)
    ma = ma[~ma["SNP"].duplicated(keep=False)]  # 12,123 rsIDs appear with two allele pairs
    ref = mp.merge(ma, on="SNP", how="inner")
    ref = ref.set_index("key")
    assert ref.index.is_unique
    cols = ["chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "low_confidence_EUR"]
    kept, read = [], 0
    for ch in pd.read_csv(src, sep="\t", usecols=cols, dtype={"chr": str}, compression="gzip",
                          chunksize=2_000_000):
        read += len(ch)
        ch = ch[ch["low_confidence_EUR"].astype(str).str.lower() == "false"]
        ch = ch.dropna(subset=["beta_EUR", "se_EUR", "neglog10_pval_EUR"])
        ch = ch.assign(key=ch["chr"] + ":" + ch["pos"].astype(str))
        ch = ch[ch["key"].isin(ref.index)]
        r = ref.loc[ch["key"]]
        ok = ((ch["alt"].values == r["A1"].values) & (ch["ref"].values == r["A2"].values)) | \
             ((ch["alt"].values == r["A2"].values) & (ch["ref"].values == r["A1"].values))
        ch = ch[ok]
        r = r[ok]
        p = np.power(10.0, -ch["neglog10_pval_EUR"].astype(float).values)
        kept.append(pd.DataFrame({"SNP": r["SNP"].values, "A1": ch["alt"].values, "A2": ch["ref"].values,
                                  "BETA": ch["beta_EUR"].values, "P": np.maximum(p, 1e-300),
                                  "N": float(n_eur)}))
    df = pd.concat(kept, ignore_index=True)
    n_dup = int(df["SNP"].duplicated(keep=False).sum())
    df = df.drop_duplicates(subset="SNP", keep=False)
    with gzip.open(out, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)
    print(f"read {read} rows; wrote {len(df)} SNPs; dropped {n_dup} rows with duplicated rsID")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
