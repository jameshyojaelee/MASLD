#!/usr/bin/env python3
"""ENZFAT1 F2: collect LDSC rg results; family F2 BH over 5 traits; height negative control.
Usage: 42_f2_summary.py <out_dir>"""
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

out = sys.argv[1]
NAMES = {"bmi": "BMI", "tg": "triglycerides", "fi": "fasting insulin", "dpw": "alcohol intake (drinks/week)",
         "t2d": "type 2 diabetes", "height": "height (negative control)"}


def parse(log):
    lines = open(log).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("Summary of Genetic Correlation Results")][0]
    hdr = lines[i + 1].split()
    rows = []
    for l in lines[i + 2:]:
        if not l.strip() or l.startswith("Analysis finished"):
            break
        rows.append(l.split())
    d = pd.DataFrame(rows, columns=hdr)
    d["p2"] = d["p2"].str.split("/").str[-1].str.replace(".sumstats.gz", "", regex=False)
    for c in hdr[2:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


res = []
for src in ("component", "alt_ukb", "pdff"):
    d = parse(f"{out}/ldsc/rg_{src}.log")
    d.insert(0, "source", src)
    res.append(d)
r = pd.concat(res, ignore_index=True)
r["trait"] = r["p2"]
r["lo95"], r["hi95"] = r["rg"] - 1.959964 * r["se"], r["rg"] + 1.959964 * r["se"]
# Addendum 2 section 2: trait single-trait h2 z < 4 stays in the family but its rg is flagged unreliable
r["trait_h2_z"] = r["h2_obs"] / r["h2_obs_se"]
r["rg_unreliable_trait_h2z_lt_4"] = r["trait_h2_z"] < 4
fam = r[(r["source"] == "component") & r["trait"].isin(["bmi", "tg", "fi", "dpw", "t2d"])].copy()
fam["q_BH"] = multipletests(fam["p"], method="fdr_bh")[1]
fam["F2_q_lt_0.05"] = fam["q_BH"] < 0.05
r = r.merge(fam[["source", "trait", "q_BH", "F2_q_lt_0.05"]], on=["source", "trait"], how="left")
h = r[(r["source"] == "component") & (r["trait"] == "height")].iloc[0]
neg = dict(rg=h["rg"], lo95=h["lo95"], hi95=h["hi95"],
           pass_=bool((h["lo95"] <= 0 <= h["hi95"]) or abs(h["rg"]) < 0.1))
r["trait_name"] = r["trait"].map(NAMES)
r.to_csv(f"{out}/f2_rg.tsv", sep="\t", index=False)
pd.DataFrame([neg]).rename(columns={"pass_": "pass"}).to_csv(f"{out}/f2_negative_control.tsv", sep="\t",
                                                             index=False)
print(r[["source", "trait", "rg", "se", "lo95", "hi95", "p", "q_BH", "gcov_int", "h2_obs", "h2_obs_se", "trait_h2_z"]]
      .round(4).to_string(index=False))
print("height negative control:", neg)
