#!/usr/bin/env python
"""GSET8: collate primary tests (Holm over 6), AveExpr sensitivity, random-set calibration,
positive-control gene z."""
import argparse

import numpy as np
import pandas as pd

TRAITS = ["pdff", "cirrhosis"]
SETS = ["hep_up", "hep_down", "complexI"]
POS = {"PNPLA3": "ENSG00000100344", "TM6SF2": "ENSG00000213996", "GCKR": "ENSG00000084734",
       "APOE": "ENSG00000130203", "MBOAT7": "ENSG00000125505", "GPAM": "ENSG00000119927",
       "MTARC1 (MARC1)": "ENSG00000186205", "TRIB1": "ENSG00000173334", "HSD17B13": "ENSG00000170509"}


def read_gsa(path):
    d = pd.read_csv(path, sep=r"\s+", comment="#")
    if "FULL_NAME" in d.columns:
        d["VARIABLE"] = d["FULL_NAME"]
    return d


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, (m - k) * p[i])
        adj[i] = min(1.0, run)
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = ap.parse_args().out
    cols = ["trait", "set", "NGENES", "BETA", "BETA_STD", "SE", "P"]

    prim, cal, sens = [], [], []
    for t in TRAITS:
        d = read_gsa(f"{out}/{t}.primary.gsa.out")
        d = d[d["TYPE"] == "SET"] if "TYPE" in d.columns else d
        for s in SETS:
            r = d[d["VARIABLE"] == s].iloc[0]
            rnd = d[d["VARIABLE"].str.startswith(f"rnd_{s}_")]
            assert len(rnd) == 1000, (t, s, len(rnd))
            prim.append(dict(trait=t, set=s, **{c: r[c] for c in cols[2:]}))
            cal.append(dict(
                trait=t, set=s, n_random=len(rnd), random_ngenes_median=int(rnd["NGENES"].median()),
                frac_p_lt_005=(rnd["P"] < 0.05).mean(), frac_p_lt_001=(rnd["P"] < 0.01).mean(),
                random_p_median=rnd["P"].median(),
                observed_p=r["P"], empirical_p=(1 + (rnd["P"] <= r["P"]).sum()) / (1 + len(rnd))))
        e = read_gsa(f"{out}/{t}.sens_aveexpr.gsa.out")
        e = e[e["TYPE"] == "SET"] if "TYPE" in e.columns else e
        for s in SETS:
            r = e[e["VARIABLE"] == s].iloc[0]
            sens.append(dict(trait=t, set=s, **{c: r[c] for c in cols[2:]}))

    prim = pd.DataFrame(prim)
    prim["P_holm6"] = holm(prim["P"])
    prim.to_csv(f"{out}/results_primary.tsv", sep="\t", index=False)
    sens = pd.DataFrame(sens)
    sens.to_csv(f"{out}/results_sensitivity_aveexpr.tsv", sep="\t", index=False)
    cal = pd.DataFrame(cal)
    cal.to_csv(f"{out}/calibration_random_sets.tsv", sep="\t", index=False)

    pc = []
    for t in TRAITS:
        g = pd.read_csv(f"{out}/{t}.genes.out", sep=r"\s+")
        g["rank"] = g["P"].rank(method="min")
        for name, eid in POS.items():
            r = g[g["GENE"] == eid]
            if len(r):
                r = r.iloc[0]
                pc.append(dict(trait=t, gene=name, ensg=eid, NSNPS=r["NSNPS"], ZSTAT=r["ZSTAT"], P=r["P"],
                               rank=int(r["rank"]), n_genes=len(g)))
            else:
                pc.append(dict(trait=t, gene=name, ensg=eid, NSNPS=np.nan, ZSTAT=np.nan, P=np.nan,
                               rank=np.nan, n_genes=len(g)))
    pc = pd.DataFrame(pc)
    pc.to_csv(f"{out}/positive_control_genes.tsv", sep="\t", index=False)

    pd.set_option("display.width", 200)
    for name, df in [("PRIMARY", prim), ("SENSITIVITY AveExpr", sens), ("CALIBRATION", cal),
                     ("POSITIVE CONTROL", pc)]:
        print(f"== {name}\n{df.to_string(index=False)}\n")


if __name__ == "__main__":
    main()
