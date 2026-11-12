#!/usr/bin/env python
"""Collate GSET8 follow-up MAGMA set results into one table plus MVP calibration and positive-control genes."""
import argparse
import glob
import os

import pandas as pd

POS = {"PNPLA3": "ENSG00000100344", "TM6SF2": "ENSG00000213996", "HSD17B13": "ENSG00000170509",
       "APOE": "ENSG00000130203", "GCKR": "ENSG00000084734", "MBOAT7": "ENSG00000125505"}


def read_gsa(path):
    d = pd.read_csv(path, sep=r"\s+", comment="#")
    if "FULL_NAME" in d.columns:
        d["VARIABLE"] = d["FULL_NAME"]
    return d[d["TYPE"] == "SET"] if "TYPE" in d.columns else d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = ap.parse_args().out
    rows = []
    for f in sorted(glob.glob(f"{out}/*.gsa.out")):
        trait, tag = os.path.basename(f).split(".")[:2]
        d = read_gsa(f)
        d = d[~d.VARIABLE.str.startswith("rnd_")]
        for _, r in d.iterrows():
            rows.append(dict(trait=trait, analysis=tag, set=r.VARIABLE, NGENES=r.NGENES, BETA=r.BETA,
                             BETA_STD=r.BETA_STD, SE=r.SE, P=r.P))
    res = pd.DataFrame(rows)
    res.to_csv(f"{out}/followup_results.tsv", sep="\t", index=False)
    d = read_gsa(f"{out}/mvp_cirrhosis.primary.gsa.out")
    cal = []
    for s in ["hep_up", "hep_down", "complexI"]:
        obs = d[d.VARIABLE == s].P.iloc[0]
        rnd = d[d.VARIABLE.str.startswith(f"rnd_{s}_")].P
        cal.append(dict(trait="mvp_cirrhosis", set=s, n_random=len(rnd), frac_p_lt_005=(rnd < 0.05).mean(),
                        frac_p_lt_001=(rnd < 0.01).mean(), observed_p=obs,
                        empirical_p=(1 + (rnd <= obs).sum()) / (1 + len(rnd))))
    cal = pd.DataFrame(cal)
    cal.to_csv(f"{out}/mvp_calibration_random_sets.tsv", sep="\t", index=False)
    g = pd.read_csv(f"{out}/mvp_cirrhosis.genes.out", sep=r"\s+")
    g["rank"] = g.P.rank(method="min")
    pc = g[g.GENE.isin(POS.values())].assign(gene=lambda x: x.GENE.map({v: k for k, v in POS.items()}))
    pc[["gene", "GENE", "NSNPS", "ZSTAT", "P", "rank"]].to_csv(f"{out}/mvp_positive_control_genes.tsv",
                                                              sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(res.to_string(index=False)); print(cal.to_string(index=False)); print(pc.to_string(index=False))


if __name__ == "__main__":
    main()
