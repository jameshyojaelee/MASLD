#!/usr/bin/env python
"""H21 step 4: family BH, cirrhosis-minus-PDFF difference test, controls.

Reads MAGMA .gsa.out tables written by 03_geneset.sbatch. Rules are those in
PRESPEC_H21.md; nothing here is tuned on the outcome.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

FAMILY = ["alt", "pdff", "nafld", "cirrhosis"]
OUTSIDE = ["cirrhosis_mvp", "bmi"]
MIN_GENES = 10
POS_PRIMARY = "hotspot_hepatocytes_e359eda2adde8948"
POS_PANEL_MODULES = ["5", "10", "16", "19", "28"]


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def read_gsa(fp: Path) -> pd.DataFrame:
    d = pd.read_csv(fp, sep=r"\s+", comment="#")
    d = d[d["TYPE"] == "SET"].copy()
    name = d["FULL_NAME"] if "FULL_NAME" in d.columns else d["VARIABLE"]
    d["set"] = name.values
    return d[["set", "NGENES", "BETA", "BETA_STD", "SE", "P"]]


def gene_residuals(out: Path, trait: str, covar: pd.DataFrame, include: set[str]) -> pd.Series:
    """Gene-Z residuals after technical covariates and the two expression
    covariates. Uses MAGMA's own base residual when the per-gene file has it."""
    gi = out / "geneset" / f"{trait}.gsa.genes.out"
    src = None
    if gi.exists():
        g = pd.read_csv(gi, sep=r"\s+", comment="#", dtype={"GENE": str})
        rc = [c for c in g.columns if "RESID" in c.upper()]
        if rc:
            g["y"] = g[rc[0]]
            src = f"{gi.name}:{rc[0]}"
            tech = []
    if src is None:
        g = pd.read_csv(out / "gene_analysis" / f"{trait}.genes.out", sep=r"\s+",
                        dtype={"GENE": str})
        z = g["ZSTAT"].to_numpy(float)
        z = np.clip(z, -3 * 1.0, z.mean() + 6 * z.std())  # MAGMA default outlier bounds
        g["y"] = z
        g["log_nsnps"] = np.log(g["NSNPS"])
        g["dens"] = g["NPARAM"] / g["NSNPS"]
        g["log_dens"] = np.log(g["dens"])
        g["inv_n"] = 1.0 / g["N"]
        tech = ["NSNPS", "log_nsnps", "dens", "log_dens", "inv_n"]
        src = f"{trait}.genes.out:ZSTAT with own technical covariates"
    g = g[g["GENE"].isin(include)].merge(covar, left_on="GENE", right_on="GENE", how="inner")
    X = g[tech + ["hep_spec", "mean_expr"]].to_numpy(float)
    X = X[:, X.std(axis=0) > 0]
    X = np.column_stack([np.ones(len(g)), X])
    y = g["y"].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = pd.Series(y - X @ beta, index=g["GENE"].values, name=trait)
    r.attrs["source"] = src
    return r


def load_sets(fp: Path) -> dict[str, list[str]]:
    sets = {}
    with open(fp) as fh:
        for line in fh:
            parts = line.split()
            sets[parts[0]] = parts[1:]
    return sets


def difference(tab_c: pd.DataFrame, tab_p: pd.DataFrame, rho: float, sd_c: float, sd_p: float):
    m = tab_c.merge(tab_p, on="set", suffixes=("_c", "_p"))
    m["d"] = m["BETA_c"] - m["BETA_p"]
    m["se_d"] = np.sqrt(m["SE_c"] ** 2 + m["SE_p"] ** 2 - 2 * rho * m["SE_c"] * m["SE_p"])
    m["z_d"] = m["d"] / m["se_d"]
    m["p_d_normal"] = 2 * stats.norm.sf(np.abs(m["z_d"]))
    bc, bp = m["BETA_c"] / sd_c, m["BETA_p"] / sd_p
    sc, sp = m["SE_c"] / sd_c, m["SE_p"] / sd_p
    m["d_scaled"] = bc - bp
    m["z_d_scaled"] = m["d_scaled"] / np.sqrt(sc ** 2 + sp ** 2 - 2 * rho * sc * sp)
    m["p_d_scaled_normal"] = 2 * stats.norm.sf(np.abs(m["z_d_scaled"]))
    m["pdff_upper95"] = m["BETA_p"] + 1.959964 * m["SE_p"]
    return m


def ols_paired(resid_c: pd.Series, resid_p: pd.Series, sets: dict, names) -> pd.DataFrame:
    """Cross-check that ignores gene-gene LD: regress the per-gene residual
    difference on set membership."""
    common = resid_c.index.intersection(resid_p.index)
    dd = (resid_c.loc[common] - resid_p.loc[common])
    rows = []
    for s in names:
        ind = dd.index.isin(sets[s])
        n1 = int(ind.sum())
        if n1 < 2:
            rows.append({"set": s, "d_ols": np.nan, "z_d_ols": np.nan})
            continue
        a, b = dd[ind], dd[~ind]
        sp2 = ((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(dd) - 2)
        se = np.sqrt(sp2 * (1 / len(a) + 1 / len(b)))
        rows.append({"set": s, "d_ols": a.mean() - b.mean(), "z_d_ols": (a.mean() - b.mean()) / se})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    res = out / "results"
    res.mkdir(exist_ok=True)

    man = pd.read_csv(out / "inputs/program_set_manifest.tsv", sep="\t", dtype=str)
    covar = pd.read_csv(out / "inputs/gene_covar.tsv", sep="\t", dtype={"GENE": str})
    include = set(pd.read_csv(out / "inputs/gene_include.txt", header=None, dtype=str)[0])
    sets = load_sets(out / "inputs/all_sets.txt")
    gsa = {t: read_gsa(out / "geneset" / f"{t}.gsa.out") for t in FAMILY + OUTSIDE}

    # ---------- per-trait program table and the 468-test family ----------
    long = []
    for t, d in gsa.items():
        d = d.assign(trait=t, is_random=d["set"].str.startswith("RND__"))
        long.append(d)
    long = pd.concat(long, ignore_index=True)
    prog = long[~long["is_random"]].merge(
        man[["program_uid", "cell_type", "module", "module_name"]],
        left_on="set", right_on="program_uid", how="left")
    prog["testable"] = prog["NGENES"] >= MIN_GENES
    fam = prog["trait"].isin(FAMILY) & prog["testable"]
    prog["q_family"] = np.nan
    prog.loc[fam, "q_family"] = bh(prog.loc[fam, "P"].to_numpy())
    for t in OUTSIDE:
        k = (prog["trait"] == t) & prog["testable"]
        prog.loc[k, "q_within_trait_outside_family"] = bh(prog.loc[k, "P"].to_numpy())
    prog.to_csv(res / "program_by_trait.tsv", sep="\t", index=False)

    n_reported = prog[prog.trait == "pdff"]["set"].nunique()
    n_testable = int(prog[(prog.trait == "pdff") & prog.testable].shape[0])
    per_trait = {t: int(((prog.trait == t) & (prog.q_family < 0.05)).sum()) for t in FAMILY}
    per_trait_out = {t: int(((prog.trait == t) & (prog.q_within_trait_outside_family < 0.05)).sum())
                     for t in OUTSIDE}

    # ---------- residual correlation and difference tests ----------
    resid = {t: gene_residuals(out, t, covar, include) for t in ["cirrhosis", "pdff", "cirrhosis_mvp"]}
    summary_diff = {}
    diff_tables = {}
    for label, tc in [("ghouse", "cirrhosis"), ("mvp", "cirrhosis_mvp")]:
        common = resid[tc].index.intersection(resid["pdff"].index)
        rho = float(np.corrcoef(resid[tc].loc[common], resid["pdff"].loc[common])[0, 1])
        sd_c, sd_p = float(resid[tc].loc[common].std()), float(resid["pdff"].loc[common].std())
        m = difference(gsa[tc], gsa["pdff"], rho, sd_c, sd_p)
        m = m.merge(ols_paired(resid[tc], resid["pdff"], sets, m["set"].tolist()), on="set")
        m["is_random"] = m["set"].str.startswith("RND__")
        rnd = m[m.is_random & (m.NGENES_c >= MIN_GENES)]
        n_rnd = len(rnd)
        k_rnd = int((rnd.p_d_normal < 0.05).sum())
        ci = stats.binomtest(k_rnd, n_rnd).proportion_ci(0.95, method="exact")
        calibrated = bool(ci.low <= 0.05 <= ci.high)
        zr = np.abs(rnd["z_d"].to_numpy())
        m["p_d_empirical"] = [(1 + (zr >= abs(z)).sum()) / (1 + n_rnd) for z in m["z_d"]]
        m["p_d_primary"] = m["p_d_normal"] if calibrated else m["p_d_empirical"]
        pg = m[~m.is_random & (m.NGENES_c >= MIN_GENES)].copy()
        pg["q_d"] = bh(pg["p_d_primary"].to_numpy())
        pg["q_d_scaled"] = bh(pg["p_d_scaled_normal"].to_numpy())
        pg = pg.merge(man[["program_uid", "cell_type", "module", "module_name"]],
                      left_on="set", right_on="program_uid", how="left")
        if label == "ghouse":
            qc = prog[(prog.trait == "cirrhosis")].set_index("set")["q_family"]
        else:
            qc = prog[(prog.trait == "cirrhosis_mvp")].set_index("set")["q_within_trait_outside_family"]
        pg["q_cirrhosis"] = pg["set"].map(qc)
        pg["supported"] = (pg.q_cirrhosis < 0.05) & (pg.d > 0) & (pg.q_d < 0.05) & (pg.P_p >= 0.05)
        pg.sort_values("p_d_primary").to_csv(
            res / f"difference_cirrhosis_{label}_minus_pdff.tsv", sep="\t", index=False)
        rnd.to_csv(res / f"difference_random_sets_{label}.tsv", sep="\t", index=False)
        diff_tables[label] = pg
        summary_diff[label] = {
            "rho_residual_geneZ": rho, "sd_resid_cirrhosis": sd_c, "sd_resid_pdff": sd_p,
            "n_genes_for_rho": int(len(common)),
            "residual_source_cirrhosis": resid[tc].attrs["source"],
            "residual_source_pdff": resid["pdff"].attrs["source"],
            "random_sets_n": n_rnd, "random_sets_p_lt_0.05": k_rnd,
            "random_rate": k_rnd / n_rnd, "random_rate_ci95": [ci.low, ci.high],
            "normal_p_calibrated": calibrated,
            "primary_p": "normal" if calibrated else "empirical from random sets",
            "n_programs_tested": int(len(pg)),
            "n_difference_q_lt_0.05": int((pg.q_d < 0.05).sum()),
            "n_difference_q_lt_0.05_positive": int(((pg.q_d < 0.05) & (pg.d > 0)).sum()),
            "n_difference_scaled_q_lt_0.05_positive": int(((pg.q_d_scaled < 0.05) & (pg.d_scaled > 0)).sum()),
            "min_q_d": float(pg.q_d.min()),
            "n_supported": int(pg.supported.sum()),
            "supported": pg.loc[pg.supported, ["set", "cell_type", "module", "module_name"]].to_dict("records"),
            "cor_z_d_magma_vs_ols": float(pg[["z_d", "z_d_ols"]].corr().iloc[0, 1]),
        }

    # ---------- random-set calibration per trait ----------
    cal = []
    rnd_long = long[long.is_random & (long.NGENES >= MIN_GENES)]
    for t, d in rnd_long.groupby("trait"):
        k, n = int((d.P < 0.05).sum()), len(d)
        ci = stats.binomtest(k, n).proportion_ci(0.95, method="exact")
        cal.append({"trait": t, "n_random_sets": n, "n_p_lt_0.05": k, "rate": k / n,
                    "ci95_low": ci.low, "ci95_high": ci.high, "min_p": d.P.min(),
                    "mean_beta": d.BETA.mean()})
    cal = pd.DataFrame(cal)
    cal.to_csv(res / "random_set_calibration.tsv", sep="\t", index=False)

    # ---------- positive control ----------
    pc = prog[(prog.trait == "pdff") & (prog.cell_type == "hepatocytes")
              & prog.module.isin(POS_PANEL_MODULES)].copy()
    pc["primary"] = pc["set"] == POS_PRIMARY
    pc["passes_0.01"] = pc["P"] < 0.01
    pc.to_csv(res / "positive_control_pdff.tsv", sep="\t", index=False)

    summ = {
        "family_traits": FAMILY, "n_programs_reported_by_magma": int(n_reported),
        "n_programs_testable": n_testable, "n_family_tests": int(fam.sum()),
        "n_BH_q_lt_0.05_per_trait_family": per_trait,
        "n_BH_q_lt_0.05_within_trait_outside_family": per_trait_out,
        "min_family_q_per_trait": {t: float(prog[(prog.trait == t)].q_family.min()) for t in FAMILY},
        "difference": summary_diff,
        "random_set_calibration": cal.to_dict("records"),
        "positive_control_pdff": {
            "primary_program": POS_PRIMARY,
            "primary_p": float(pc.loc[pc.primary, "P"].iloc[0]) if pc.primary.any() else None,
            "panel_n": int(len(pc)), "panel_min_p": float(pc.P.min()) if len(pc) else None,
            "passes": bool(pc["passes_0.01"].any()) if len(pc) else False},
    }
    (res / "SUMMARY.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
