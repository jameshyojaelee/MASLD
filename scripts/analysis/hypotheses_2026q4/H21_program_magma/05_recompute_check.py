#!/usr/bin/env python
"""H21 step 5: recompute reported numbers from the raw MAGMA tables.

Independent of 04: own parser, statsmodels BH, p recomputed from BETA/SE,
rho recomputed from the gene-level files. Exits non-zero on disagreement.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

FAMILY = ["alt", "pdff", "nafld", "cirrhosis"]


def parse(fp: Path) -> pd.DataFrame:
    rows, header = [], None
    for line in open(fp):
        if line.startswith("#"):
            continue
        f = line.split()
        if header is None:
            header = f
            continue
        rows.append(dict(zip(header, f)))
    d = pd.DataFrame(rows)
    d = d[d.TYPE == "SET"]
    d["name"] = d["FULL_NAME"] if "FULL_NAME" in d else d["VARIABLE"]
    for c in ("NGENES", "BETA", "SE", "P"):
        d[c] = d[c].astype(float)
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)
    summ = json.loads((out / "results/SUMMARY.json").read_text())
    ok = True
    tabs = {t: parse(out / "geneset" / f"{t}.gsa.out") for t in FAMILY}
    fam = []
    for t, d in tabs.items():
        d = d[~d.name.str.startswith("RND__") & (d.NGENES >= 10)]
        fam.append(d.assign(trait=t))
        ndf = None
    fam = pd.concat(fam, ignore_index=True)
    fam["q"] = multipletests(fam.P, method="fdr_bh")[1]
    print(f"family tests: {len(fam)} (reported {summ['n_family_tests']})")
    ok &= len(fam) == summ["n_family_tests"]
    for t in FAMILY:
        n = int(((fam.trait == t) & (fam.q < 0.05)).sum())
        rep = summ["n_BH_q_lt_0.05_per_trait_family"][t]
        nominal = int(((fam.trait == t) & (fam.P < 0.05)).sum())
        print(f"{t}: BH q<0.05 recomputed {n}, reported {rep}; nominal p<0.05 (not significant) {nominal}; "
              f"min p {fam[fam.trait == t].P.min():.3g}; min q {fam[fam.trait == t].q.min():.3g}")
        ok &= n == rep
    # p recomputed from beta/se, one-sided normal approximation to MAGMA's t
    z = fam.BETA / fam.SE
    p_norm = stats.norm.sf(z)
    dev = np.max(np.abs(np.log10(p_norm) - np.log10(fam.P)))
    print(f"max |log10 p(normal from BETA/SE) - log10 p(MAGMA)| = {dev:.3f}")
    ok &= dev < 0.15

    # difference test, recomputed
    c = tabs["cirrhosis"].set_index("name")
    p = tabs["pdff"].set_index("name")
    rho = summ["difference"]["ghouse"]["rho_residual_geneZ"]
    names = [n for n in c.index if not n.startswith("RND__") and c.loc[n, "NGENES"] >= 10]
    d = c.loc[names, "BETA"] - p.loc[names, "BETA"]
    se = np.sqrt(c.loc[names, "SE"] ** 2 + p.loc[names, "SE"] ** 2
                 - 2 * rho * c.loc[names, "SE"] * p.loc[names, "SE"])
    pz = 2 * stats.norm.sf(np.abs(d / se))
    rep = pd.read_csv(out / "results/difference_cirrhosis_ghouse_minus_pdff.tsv", sep="\t").set_index("set")
    if summ["difference"]["ghouse"]["normal_p_calibrated"]:
        q = multipletests(pz, method="fdr_bh")[1]
        n_pos = int(((q < 0.05) & (d.values > 0)).sum())
        print(f"difference q<0.05 and d>0: recomputed {n_pos}, reported "
              f"{summ['difference']['ghouse']['n_difference_q_lt_0.05_positive']}")
        ok &= n_pos == summ["difference"]["ghouse"]["n_difference_q_lt_0.05_positive"]
    md = float(np.max(np.abs(d.values - rep.loc[names, "d"].values)))
    print(f"max |d recomputed - d reported| = {md:.2e}")
    ok &= md < 1e-9
    qfam = fam[fam.trait == "cirrhosis"].set_index("name")["q"]
    sup = [n for n in names if qfam[n] < 0.05 and d[n] > 0 and rep.loc[n, "q_d"] < 0.05 and p.loc[n, "P"] >= 0.05]
    rep_sup = sorted(x["set"] for x in summ["difference"]["ghouse"]["supported"])
    print(f"supported programs recomputed: {sorted(sup)}; reported: {rep_sup}")
    ok &= sorted(sup) == rep_sup

    # rho recomputed from MAGMA per-gene files when present
    try:
        gc = pd.read_csv(out / "geneset/cirrhosis.gsa.genes.out", sep=r"\s+", comment="#", dtype={"GENE": str})
        gp = pd.read_csv(out / "geneset/pdff.gsa.genes.out", sep=r"\s+", comment="#", dtype={"GENE": str})
        col = [x for x in gc.columns if "RESID" in x.upper()] or ["ZSTAT"]
        m = gc[["GENE", col[0]]].merge(gp[["GENE", col[0]]], on="GENE")
        r = np.corrcoef(m.iloc[:, 1], m.iloc[:, 2])[0, 1]
        print(f"gene-level correlation of {col[0]} (cirrhosis, PDFF), {len(m)} genes, before expression "
              f"covariates: {r:.4f}; rho used: {rho:.4f}")
    except Exception as e:  # file layout differs
        print("rho cross-check skipped:", e)
    print("RECOMPUTE CHECK:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
