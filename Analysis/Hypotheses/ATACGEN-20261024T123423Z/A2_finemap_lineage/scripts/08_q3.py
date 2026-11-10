#!/usr/bin/env python3
"""A2 Q3: posterior in MASH-responsive vs stable peaks (hepatocyte, stellate), conditional on
falling in jointly testable DA peaks. Family Q3: 2 lineages x 2 universes, BH."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from a2common import RES, SEED, W, bh, boot_p, ci

MIN_POS = 5  # Amendment A5: units with posterior in responsive and in stable peaks

DL = ["hepatocyte", "stellate"]
CLS = ["responsive", "stable", "resp_up", "resp_down"]


def flags(df, lin):
    d = df[f"da_{lin}"].astype(str)
    dr = df[f"dadir_{lin}"].astype(str)
    return {"tested": d.isin(["responsive", "middle", "stable"]).to_numpy(),
            "responsive": (d == "responsive").to_numpy(), "stable": (d == "stable").to_numpy(),
            "resp_up": ((d == "responsive") & (dr == "up")).to_numpy(),
            "resp_down": ((d == "responsive") & (dr == "down")).to_numpy()}


def per_window(df, key, wcol, lin, q_by_bin):
    """O and E (TSS-matched pooled q, and within-window) per window for one lineage."""
    f = flags(df, lin)
    w = df[wcol].to_numpy(float)
    tb = df.tss_bin.to_numpy()
    out = pd.DataFrame({key: df[key].to_numpy()})
    out["M_tested"] = w * f["tested"]
    out["n_tested"] = f["tested"].astype(float)
    for c in CLS:
        out[f"O_{c}"] = w * f[c]
        out[f"E_{c}"] = w * f["tested"] * q_by_bin[c][tb]
        out[f"n_{c}"] = f[c].astype(float)
    g = out.groupby(key).sum()
    for c in CLS:
        with np.errstate(divide="ignore", invalid="ignore"):
            g[f"Ew_{c}"] = np.where(g.n_tested > 0, g.M_tested * g[f"n_{c}"] / g.n_tested, 0.0)
    return g


def pooled_q(df, lin):
    f = flags(df, lin)
    t = df.tss_bin.to_numpy()
    q = {}
    for c in CLS:
        num = np.bincount(t[f["tested"]], weights=f[c][f["tested"]].astype(float), minlength=5)
        den = np.bincount(t[f["tested"]], minlength=5).astype(float)
        q[c] = np.where(den > 0, num / np.maximum(den, 1), 0.0)
    return q, np.bincount(t[f["tested"]], minlength=5)


def ratio_boot(units, rng, nboot=2000, ecol="E"):
    O_R, E_R = units["O_responsive"].to_numpy(), units[f"{ecol}_responsive"].to_numpy()
    O_S, E_S = units["O_stable"].to_numpy(), units[f"{ecol}_stable"].to_numpy()
    n = len(units)
    C = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(nboot)])
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = (C @ O_R) / (C @ E_R)
        rs = (C @ O_S) / (C @ E_S)
    return (O_R.sum() / E_R.sum(), O_S.sum() / E_S.sum(),
            (O_R.sum() / E_R.sum()) / (O_S.sum() / E_S.sum()), rr, rs, rr / rs)


rows = []
# primary: COLOC r2 primary pairs, gene units
cv = pd.read_parquet(W / "q2_coloc_annot.parquet")
# secondary: Q1 liver-class credible sets, LD-block units
mem = pd.read_parquet(W / "q1_members_annot.parquet")
win = pd.read_parquet(W / "q1_windows_annot.parquet")
mem = mem[mem.trait_class.isin(["direct", "liver_enzyme", "sequela"])]
win = win[win.wid.isin(mem.wid.unique())]
for lin in DL:
    rng = np.random.default_rng(SEED)
    # pool for q: each pair's variant set once
    q, nb = pooled_q(cv, lin)
    g = per_window(cv, "pid", "w", lin, q)
    g = g.join(cv.drop_duplicates("pid").set_index("pid")[["gene"]])
    units = g.groupby("gene").mean()
    for ecol, tag in (("E", "tss_matched_pooled"), ("Ew", "within_window")):
        r_r, r_s, R, rr, rs, RR = ratio_boot(units, rng, ecol=ecol)
        rows.append({"universe": "coloc_r2_primary_pairs", "unit": "gene", "lineage": lin,
                     "expectation": tag, "n_units": len(units),
                     "n_units_with_tested_mass": int((units.M_tested > 0).sum()),
                     "mass_in_tested_peaks": units.M_tested.mean(),
                     "obs_responsive": units.O_responsive.mean(),
                     "exp_responsive": units[f"{ecol}_responsive"].mean(),
                     "obs_stable": units.O_stable.mean(), "exp_stable": units[f"{ecol}_stable"].mean(),
                     "oe_responsive": r_r, "oe_stable": r_s, "R": R,
                     "R_lo": ci(RR)[0], "R_hi": ci(RR)[1],
                     "n_units_obs_responsive": int((units.O_responsive > 0).sum()),
                     "n_units_obs_stable": int((units.O_stable > 0).sum()),
                     "p_boot": boot_p(np.log(RR)) if min((units.O_responsive > 0).sum(),
                                                       (units.O_stable > 0).sum()) >= MIN_POS
                     else np.nan,
                     "oe_up": units.O_resp_up.sum() / units[f"{ecol}_resp_up"].sum(),
                     "oe_down": units.O_resp_down.sum() / units[f"{ecol}_resp_down"].sum(),
                     "tested_variants_by_tss_bin": ",".join(map(str, nb))})
    # secondary
    q2, nb2 = pooled_q(win, lin)
    gw = per_window(mem, "cs_uid", "w", lin, q2)
    # within-window expectation for CSs uses their own window variants
    fw = flags(win, lin)
    wn = pd.DataFrame({"wid": win.wid.to_numpy(), "n_tested": fw["tested"].astype(float),
                       **{f"n_{c}": fw[c].astype(float) for c in CLS}}).groupby("wid").sum()
    cmap = mem.drop_duplicates("cs_uid").set_index("cs_uid")
    gw = gw.join(cmap[["wid", "block"]])
    for c in CLS:
        nt = gw.wid.map(wn.n_tested).to_numpy()
        nc = gw.wid.map(wn[f"n_{c}"]).to_numpy()
        gw[f"Ew_{c}"] = np.where(nt > 0, gw.M_tested * nc / np.maximum(nt, 1), 0.0)
    units = gw.drop(columns=["wid"]).groupby("block").mean()
    for ecol, tag in (("E", "tss_matched_pooled"), ("Ew", "within_window")):
        r_r, r_s, R, rr, rs, RR = ratio_boot(units, rng, ecol=ecol)
        rows.append({"universe": "q1_liver_credible_sets", "unit": "LD block", "lineage": lin,
                     "expectation": tag, "n_units": len(units),
                     "n_units_with_tested_mass": int((units.M_tested > 0).sum()),
                     "mass_in_tested_peaks": units.M_tested.mean(),
                     "obs_responsive": units.O_responsive.mean(),
                     "exp_responsive": units[f"{ecol}_responsive"].mean(),
                     "obs_stable": units.O_stable.mean(), "exp_stable": units[f"{ecol}_stable"].mean(),
                     "oe_responsive": r_r, "oe_stable": r_s, "R": R,
                     "R_lo": ci(RR)[0], "R_hi": ci(RR)[1],
                     "n_units_obs_responsive": int((units.O_responsive > 0).sum()),
                     "n_units_obs_stable": int((units.O_stable > 0).sum()),
                     "p_boot": boot_p(np.log(RR)) if min((units.O_responsive > 0).sum(),
                                                       (units.O_stable > 0).sum()) >= MIN_POS
                     else np.nan,
                     "oe_up": units.O_resp_up.sum() / units[f"{ecol}_resp_up"].sum(),
                     "oe_down": units.O_resp_down.sum() / units[f"{ecol}_resp_down"].sum(),
                     "tested_variants_by_tss_bin": ",".join(map(str, nb2))})
res = pd.DataFrame(rows)
fam = res.expectation == "tss_matched_pooled"
res["family"] = np.where(fam, "Q3 (BH over 4)", "sensitivity")
res["q_bh"] = np.nan
res.loc[fam, "q_bh"] = bh(res.loc[fam, "p_boot"].fillna(1.0))
res.loc[res.p_boot.isna(), "q_bh"] = np.nan
res.to_csv(RES / "q3_responsive_vs_stable.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(res.drop(columns=["tested_variants_by_tss_bin"]).round(4).to_string())

# variant-level listing of colocalized posterior in responsive peaks (context)
lst = []
da = pd.read_csv(W / "da_peaks.tsv", sep="\t")
for lin in DL:
    s = cv[(cv[f"da_{lin}"].astype(str) == "responsive") & (cv.w >= 0.05)]
    for r in s.itertuples():
        lst.append({"lineage": lin, "gene": r.gene, "gwas_name": r.gwas_name, "pp_h4": r.pp_h4,
                    "variant_hg19": f"{r.chromosome}:{r.position}:{r.allele1}:{r.allele2}",
                    "hg38": f"chr{r.chromosome}:{r.pos38}", "snp_pp_h4_norm": r.w,
                    "direction": getattr(r, f"dadir_{lin}")})
lst = pd.DataFrame(lst)
if not lst.empty:
    lst = lst.sort_values(["lineage", "snp_pp_h4_norm"], ascending=[True, False])
lst.to_csv(RES / "q3_coloc_variants_in_responsive_peaks.tsv", sep="\t", index=False)
print(lst.head(40).to_string())
