"""Independent recomputation of Figure 1B values from the adopted bulk fit and COLOC r2.

Inputs (read-only):
  bulk fit:  RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv
  COLOC r2:  GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv
Definitions follow docs/manuscript/05_figure_legends.md (Fig 1B) and PRESPEC_H57.md (universe, bins).
Nothing here reads the panel's own source table or plotting script.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = ROOT / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1"
DEG = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
DEG_REPRO = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
SEED = 20261017


def load(deg_path):
    deg = pd.read_csv(deg_path)
    deg["ens"] = deg["gene"].str.replace(r"\.\d+$", "", regex=True)
    assert deg["ens"].is_unique
    col = pd.read_csv(COLOC)
    assert col["ensembl"].is_unique
    deg["is_deg"] = (deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.5)
    full = deg.merge(col, left_on="ens", right_on="ensembl", how="inner")
    return deg, col, full


def bins_for(u, nb):
    u = u.sort_values("ens").reset_index(drop=True)
    r = u["AveExpr"].rank(method="first")
    u["bin"] = pd.qcut(r, nb, labels=False)
    return u


def strat_stats(u, gcol, dcol, nb=10):
    u = bins_for(u, nb)
    rows = []
    for k, b in u.groupby("bin"):
        n = len(b); G = int(b[gcol].sum()); D = int(b[dcol].sum()); O = int((b[gcol] & b[dcol]).sum())
        rows.append(dict(bin=k, n=n, G=G, D=D, O=O, E=G * D / n,
                         amin=b["AveExpr"].min(), amax=b["AveExpr"].max()))
    bt = pd.DataFrame(rows)
    E = bt["E"].sum(); O = bt["O"].sum()
    tabs = [np.array([[r.O, r.G - r.O], [r.D - r.O, r.n - r.G - r.D + r.O]]) for r in bt.itertuples()]
    st = StratifiedTable(tabs)
    orr = st.oddsratio_pooled; lo, hi = st.oddsratio_pooled_confint()
    # exact convolution null of the total overlap
    dist = np.array([1.0])
    for r in bt.itertuples():
        ks = np.arange(0, min(r.G, r.D) + 1)
        pk = hypergeom.pmf(ks, r.n, r.D, r.G)
        dist = np.convolve(dist, pk)
    xs = np.arange(len(dist))
    p = dist[np.abs(xs - E) >= abs(O - E) - 1e-9].sum()
    cdf = np.cumsum(dist)
    q025 = int(xs[np.searchsorted(cdf, 0.025)]); q975 = int(xs[np.searchsorted(cdf, 0.975)])
    return dict(O=int(O), E=E, MH_OR=orr, lo=lo, hi=hi, p_exact=p, null_q025=q025, null_q975=q975), bt


def main():
    res = {}
    deg, col, full = load(DEG)
    res["n_bulk_genes"] = len(deg)
    res["n_bulk_deg_total"] = int(deg["is_deg"].sum())
    res["n_deg_up"] = int((deg["is_deg"] & (deg["logFC"] > 0)).sum())
    res["n_deg_down"] = int((deg["is_deg"] & (deg["logFC"] < 0)).sum())
    res["n_coloc_t12_genes"] = len(col)
    res["U_full"] = len(full)
    res["U_full_state_by_deg"] = pd.crosstab(full["susie_state_t12"], full["is_deg"]).to_dict()
    res["D_full"] = int(full["is_deg"].sum())
    ev = full[full["susie_state_t12"].isin(["positive", "tested_le_0.5"])].copy()
    ev["is_gen"] = ev["susie_state_t12"].eq("positive")
    # consistency of state vs PP4
    res["gen_by_pp4_gt_0.5"] = int((ev["coloc_best_susie_pp4"] > 0.5).sum())
    res["N"] = len(ev); res["G"] = int(ev["is_gen"].sum()); res["D"] = int(ev["is_deg"].sum())
    res["O"] = int((ev["is_gen"] & ev["is_deg"]).sum())
    res["gen_only"] = res["G"] - res["O"]; res["deg_only"] = res["D"] - res["O"]
    res["frac_O_of_G"] = res["O"] / res["G"]; res["frac_O_of_D"] = res["O"] / res["D"]
    res["frac_genonly_of_G"] = res["gen_only"] / res["G"]; res["frac_degonly_of_D"] = res["deg_only"] / res["D"]
    res["deg_no_susie_result_in_U_full"] = int((full["is_deg"] & full["susie_state_t12"].eq("no_susie_result")).sum())
    res["genes_no_susie_result_in_U_full"] = int(full["susie_state_t12"].eq("no_susie_result").sum())
    res["untestable_in_U_full"] = int(full["susie_state_t12"].eq("untestable").sum())
    res["E_unstrat"] = res["G"] * res["D"] / res["N"]
    t = np.array([[res["O"], res["G"] - res["O"]], [res["D"] - res["O"], res["N"] - res["G"] - res["D"] + res["O"]]])
    fo, fp = fisher_exact(t)
    res["fisher_OR_unstrat"] = fo; res["fisher_p_unstrat"] = fp
    for nb in (5, 10, 20, 50):
        s, bt = strat_stats(ev, "is_gen", "is_deg", nb)
        res[f"strat_{nb}bins"] = s
        if nb == 10:
            bt.to_csv(OUT / "recompute_1b_bins10.tsv", sep="\t", index=False)
    s10 = res["strat_10bins"]
    res["expected_frac_of_G_10bins"] = s10["E"] / res["G"]
    res["expected_frac_of_D_10bins"] = s10["E"] / res["D"]
    res["null95_frac_of_G"] = [s10["null_q025"] / res["G"], s10["null_q975"] / res["G"]]
    res["null95_frac_of_D"] = [s10["null_q025"] / res["D"], s10["null_q975"] / res["D"]]

    # NC-1 style: expression-matched lists of non-genetic genes, one-to-one nearest AveExpr, seeded
    rng = np.random.default_rng(SEED)
    gen = ev[ev["is_gen"]]; non = ev[~ev["is_gen"]].reset_index(drop=True)
    nonA = non["AveExpr"].to_numpy(); nonD = non["is_deg"].to_numpy()
    Os = []; OEs = []
    for i in range(200):
        order = rng.permutation(len(gen))
        used = np.zeros(len(non), bool); pick = []
        for a in gen["AveExpr"].to_numpy()[order]:
            d = np.abs(nonA - a); d[used] = np.inf
            j = int(np.argmin(d)); used[j] = True; pick.append(j)
        lab = np.zeros(len(ev), bool)
        evc = ev.copy(); evc["is_list"] = False
        ids = set(non.loc[pick, "ens"])
        evc["is_list"] = evc["ens"].isin(ids)
        s, _ = strat_stats(evc, "is_list", "is_deg", 10)
        Os.append(s["O"]); OEs.append(s["O"] / s["E"])
    res["nc1_n_lists"] = 200; res["nc1_O_mean"] = float(np.mean(Os)); res["nc1_O_min"] = int(min(Os)); res["nc1_O_max"] = int(max(Os))
    res["nc1_strat_OE_mean"] = float(np.mean(OEs))
    res["nc1_E_mean_implied"] = float(np.mean(Os) / np.mean(OEs))

    # display alternatives
    alt = {}
    for thr in (0.8,):
        g2 = ev["coloc_best_susie_pp4"] > thr
        alt[f"G_pp4_gt_{thr}"] = int(g2.sum()); alt[f"O_pp4_gt_{thr}"] = int((g2 & ev["is_deg"]).sum())
        e2 = ev.copy(); e2["g2"] = g2
        alt[f"strat10_pp4_gt_{thr}"] = strat_stats(e2, "g2", "is_deg", 10)[0]
    d2 = ev["padj"] < 0.05
    e2 = ev.copy(); e2["d2"] = d2
    alt["D_padj_only"] = int(d2.sum()); alt["O_padj_only"] = int((d2 & ev["is_gen"]).sum())
    alt["strat10_padj_only"] = strat_stats(e2, "is_gen", "d2", 10)[0]
    # whole interface (old universe): genes without SuSiE result counted as non-genetic
    f2 = full.copy(); f2["is_gen"] = f2["susie_state_t12"].eq("positive")
    alt["U_full_O"] = int((f2["is_gen"] & f2["is_deg"]).sum())
    alt["U_full_E_unstrat"] = f2["is_gen"].sum() * f2["is_deg"].sum() / len(f2)
    alt["U_full_strat10"] = strat_stats(f2, "is_gen", "is_deg", 10)[0]
    res["alternatives"] = alt

    # positive controls
    pcs = {}
    for sym in ["PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MBOAT7", "APOE", "MARC1", "MTARC1", "GPAM", "SERPINA1", "TRIB1"]:
        r = deg[deg["symbol"] == sym]
        if r.empty:
            pcs[sym] = "not in bulk table"; continue
        r = r.iloc[0]
        c = col[col["ensembl"] == r["ens"]]
        pcs[sym] = dict(ens=r["ens"], logFC=float(r["logFC"]), padj=float(r["padj"]), AveExpr=float(r["AveExpr"]),
                        is_deg=bool(r["is_deg"]),
                        susie_state=None if c.empty else c.iloc[0]["susie_state_t12"],
                        susie_pp4=None if c.empty else (None if pd.isna(c.iloc[0]["coloc_best_susie_pp4"]) else float(c.iloc[0]["coloc_best_susie_pp4"])),
                        abf_pp4=None if c.empty else (None if pd.isna(c.iloc[0]["coloc_best_abf_pp4"]) else float(c.iloc[0]["coloc_best_abf_pp4"])),
                        susie_driving_gwas=None if c.empty else c.iloc[0]["susie_driving_gwas"],
                        in_U_eval=bool(r["ens"] in set(ev["ens"])))
    res["positive_controls"] = pcs

    # same computation on the reproduction copy of the fit
    if DEG_REPRO.exists():
        d_r, _, full_r = load(DEG_REPRO)
        ev_r = full_r[full_r["susie_state_t12"].isin(["positive", "tested_le_0.5"])].copy()
        ev_r["is_gen"] = ev_r["susie_state_t12"].eq("positive")
        m = deg[["ens", "logFC", "padj", "AveExpr", "is_deg"]].merge(d_r[["ens", "logFC", "padj", "AveExpr", "is_deg"]], on="ens", suffixes=("", "_r"))
        res["repro_copy"] = dict(N=len(ev_r), D=int(ev_r["is_deg"].sum()), O=int((ev_r["is_gen"] & ev_r["is_deg"]).sum()),
                                 E10=strat_stats(ev_r, "is_gen", "is_deg", 10)[0]["E"],
                                 deg_calls_differ=int((m["is_deg"] != m["is_deg_r"]).sum()),
                                 max_abs_logFC_diff=float((m["logFC"] - m["logFC_r"]).abs().max()),
                                 max_abs_AveExpr_diff=float((m["AveExpr"] - m["AveExpr_r"]).abs().max()))

    ev[["ens", "symbol", "AveExpr", "logFC", "padj", "is_deg", "susie_state_t12", "coloc_best_susie_pp4", "is_gen"]].to_csv(
        OUT / "recompute_1b_U_eval_genes.tsv.gz", sep="\t", index=False)
    with open(OUT / "recompute_1b.json", "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps({k: v for k, v in res.items() if k not in ("positive_controls",)}, indent=1, default=str))
    print(json.dumps(res["positive_controls"], indent=1, default=str))


if __name__ == "__main__":
    main()
