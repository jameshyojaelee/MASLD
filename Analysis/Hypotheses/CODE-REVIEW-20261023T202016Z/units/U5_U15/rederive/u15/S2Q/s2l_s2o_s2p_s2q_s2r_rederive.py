"""Independent retabulation of Figure S2L, S2O, S2P, S2Q and S2R.

COLOC: adopted r2 master restricted to the 35 primary strata (common/coloc35.csv.gz, written by
common/extract_coloc35.py). Direct = tier 1 (15 strata), enzyme = tier 2 (20 strata).
Per gene and stratum family: SuSiE max, ABF max (finite values only), best = max of the two.
Single-signal-only star: ABF max > 0.5 and SuSiE max missing or <= 0.5.
Bulk: 774-unit fit deg_results.csv; canonical call = BH padj < 0.05 and |log2FC| > 0.5,
matched by symbol (every Ensembl record for the symbol is reported). No randomness.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
U = Path(__file__).resolve().parents[1]
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)

c = pd.read_csv(U / "common/coloc35.csv.gz")
print("coloc rows (35 strata):", len(c), "strata:", c["gwas_name"].nunique(),
      "tier counts (strata):", c.drop_duplicates("gwas_name")["tier"].value_counts().to_dict())
deg = pd.read_csv(R / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv", float_precision="round_trip")
print("bulk genes:", len(deg), "unique ENSG:", deg["gene"].nunique())
deg["canon"] = (deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.5)
print("canonical calls:", int(deg["canon"].sum()))


def fam(genes, tiers):
    sub = c[c["gene"].isin(genes) & c["tier"].isin(tiers)]
    g = sub.groupby("gene").agg(susie=("PP.H4.susie", "max"), abf=("PP.H4.abf", "max"),
                                n_strata=("gwas_name", "nunique"),
                                n_susie=("PP.H4.susie", lambda s: s.notna().sum()))
    g = g.reindex(genes)
    g["best"] = g[["susie", "abf"]].max(axis=1)
    g["state"] = np.where(g["susie"] > 0.5, "multi",
                 np.where(g["abf"] > 0.5, "single-only",
                 np.where(g["n_strata"].fillna(0) > 0, "evaluated<=0.5", "not evaluable")))
    return g


def fmt(v):
    if pd.isna(v):
        return "NA"
    if v >= 0.995:
        return "1.00"
    if 0.495 <= v < 0.505:
        return f"{v:.3f}"
    return f"{v:#.2g}"


def bulk(genes):
    rows = []
    for g in genes:
        d = deg[deg["symbol"] == g]
        if len(d) == 0:
            rows.append((g, "NA", np.nan, np.nan, "not in bulk universe (symbol)"))
        for _, r in d.iterrows():
            call = ("DEG " + ("up" if r["logFC"] > 0 else "down") + f" ({r['padj']:.0e})") if r["canon"] else "no canonical call"
            rows.append((g, r["gene"], r["logFC"], r["padj"], call))
    return pd.DataFrame(rows, columns=["gene", "ensembl", "logFC", "padj", "bulk"])


# ---------------- S2L ----------------
print("\n=== S2L F2RL1")
f = c[(c["gene"] == "F2RL1") & c["gwas_name"].isin(["UKBB_GGT", "PanUKBB_CSA_GGT", "BBJ_GGT", "PanUKBB_AFR_GGT"])]
print(f[["gwas_name", "ancestry", "PP.H4.abf", "PP.H4.susie", "method", "n_cs_pairs"]].to_string(index=False))
d = deg[deg["symbol"] == "F2RL1"].iloc[0]
z = stats.norm.ppf(0.975)
print(f"logFC {d['logFC']:.4f} SE {d['SE']:.4f} t {d['t']:.3f} (logFC/SE {d['logFC']/d['SE']:.3f}) padj {d['padj']:.3g}")
print(f"95% CI (z): {d['logFC']-z*d['SE']:.4f} to {d['logFC']+z*d['SE']:.4f}")
out = {"S2L": f}

# ---------------- S2O / S2P ----------------
print("\n=== S2O / S2P")
ann = pd.read_csv(R / "RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv", sep="\t")
six = ["MTARC1", "HSD17B13", "ERLIN1", "MTTP", "GCKR", "LGALS1"]
o = fam(six, [1, 2])
o = o.join(ann.set_index("gene")[["gobeil_pqtl_coloc_pph4", "layer_relationship", "our_max_susie_pp4", "our_max_abf_pp4"]])
o["best_fmt"] = o["best"].map(fmt)
print(o.to_string())
o.to_csv(U / "S2O/S2O_S2P_rederived.csv")

# ---------------- S2Q ----------------
print("\n=== S2Q")
el = pd.read_csv(R / "data/external/elison2025_sc_chromatin/elison2025_sc_chromatin_calls.csv")
print("Elison rows:", len(el))
q_genes = ["RORA", "PPP1R3B", "KRT8", "EFHD1", "CUX2", "RELB", "CREB5", "THRB", "HNF4A",
           "TRPS1", "FCGR2B", "HLA-DQA1", "FLACC1", "PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7"]
print("Elison genes not drawn:", sorted(set(el["gene"]) - set(q_genes)))


def table(genes):
    dd = fam(genes, [1]).add_prefix("dir_")
    ee = fam(genes, [2]).add_prefix("enz_")
    t = dd.join(ee)
    t["dir_lab"] = [fmt(v) + ("*" if s == "single-only" else "") for v, s in zip(t["dir_best"], t["dir_state"])]
    t["enz_lab"] = [fmt(v) + ("*" if s == "single-only" else "") for v, s in zip(t["enz_best"], t["enz_state"])]
    b = bulk(genes)
    return t, b


tq, bq = table(q_genes)
print(tq[["dir_susie", "dir_abf", "dir_n_strata", "dir_lab", "enz_susie", "enz_abf", "enz_n_strata", "enz_lab"]].to_string())
print(bq.to_string(index=False))
tq.to_csv(U / "S2Q/S2Q_rederived.csv")
bq.to_csv(U / "S2Q/S2Q_bulk_rederived.csv", index=False)

# ---------------- S2R ----------------
print("\n=== S2R")
zh = pd.read_csv(R / "data/external/zhu2026_masld_mpra/zhu2026_masld_mpra_calls.csv")
r_genes = list(zh["gene"])
print("Zhu genes:", len(r_genes))
tr, br = table(r_genes)
tr = tr.join(zh.set_index("gene")[["zhu_celltype", "zhu_mpra_log2fc"]])
tr["mpra_1dp"] = tr["zhu_mpra_log2fc"].round(1)
sup = []
for g, r in tr.iterrows():
    s = []
    if r["dir_best"] > 0.5:
        s.append("direct" + ("*" if r["dir_state"] == "single-only" else ""))
    if r["enz_best"] > 0.5:
        s.append("enzyme" + ("*" if r["enz_state"] == "single-only" else ""))
    bb = br[br["gene"] == g]
    if bb["bulk"].str.startswith("DEG").any():
        s.append("DEG")
    sup.append(" + ".join(s) if s else "none")
tr["support"] = sup
print(tr[["dir_susie", "dir_abf", "dir_lab", "enz_susie", "enz_abf", "enz_lab", "mpra_1dp", "zhu_celltype", "support"]].to_string())
print(br.to_string(index=False))
tr.to_csv(U / "S2R/S2R_rederived.csv")
br.to_csv(U / "S2R/S2R_bulk_rederived.csv", index=False)

# NCAN / DLG5: look for their Ensembl IDs in the bulk table under any symbol
print("\n=== NCAN / DLG5 Ensembl IDs from COLOC master and presence in bulk table")
for g in ["NCAN", "DLG5"]:
    ids = sorted(set(c.loc[c["gene"] == g, "ensembl"].dropna()))
    for i in ids:
        hit = deg[deg["gene"].str.split(".").str[0] == i.split(".")[0]]
        print(g, i, "rows in bulk table:", len(hit), hit[["gene", "symbol", "logFC", "padj"]].to_dict("records"))
