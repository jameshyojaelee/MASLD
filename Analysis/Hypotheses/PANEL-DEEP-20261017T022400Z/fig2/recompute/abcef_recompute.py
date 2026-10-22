"""Independent recomputation of Figure 2 A (COLOC columns), B, C (lead selection), E, F
from the COLOC r2 master (35 tier-1/2 studies), the fine-mapping tables and GWAS sumstats.

Inputs: upstream/ written by extract_upstream.py (job 9820).
Gene key: stable Ensembl ID (anonymous rows without an Ensembl ID are excluded).
"""
import os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
U = f"{HERE}/upstream"
OUT = f"{HERE}/out_ABCEF"
os.makedirs(OUT, exist_ok=True)
reg = pd.read_csv(f"{U}/registry_t12.tsv", sep="\t")
m = pd.read_csv(f"{U}/r2_master_t12.tsv.gz", sep="\t", dtype={"chr": str, "top_snp": str})
m["s4"] = pd.to_numeric(m["PP.H4.susie"], errors="coerce")
m["a4"] = pd.to_numeric(m["PP.H4.abf"], errors="coerce")
m = m.merge(reg[["study_name", "trait", "tier", "ancestry"]].rename(columns={"study_name": "gwas_name", "ancestry": "gwas_anc"}),
            on="gwas_name", how="left")
m["anc5"] = m.gwas_anc.replace({"CSA": "SAS"})
print("ancestry labels in registry:", m.gwas_anc.value_counts().to_dict())
g = m[m.ensembl.notna() & (m.ensembl.astype(str).str.startswith("ENSG"))].copy()
print("rows with Ensembl", len(g), "of", len(m))
summary = {}

# ---------- gene-level states ----------
gs = g.groupby("ensembl").agg(max_s4=("s4", "max"), max_a4=("a4", "max"),
                              untestable=("method", lambda s: (s == "susie_untestable_insufficient_shared_posterior").any()),
                              symbol=("gene", lambda s: s.dropna().iloc[0] if s.notna().any() else ""))
gs["susie_pos"] = gs.max_s4 > 0.5
gs["abf_pos"] = gs.max_a4 > 0.5
gs["union"] = gs.susie_pos | gs.abf_pos
gs["abf_only"] = gs.abf_pos & ~gs.susie_pos
summary["n_genes_finite_susie"] = int(gs.max_s4.notna().sum())
summary["multi_signal"] = int(gs.susie_pos.sum())
summary["abf_pos"] = int(gs.abf_pos.sum())
summary["single_signal_only"] = int(gs.abf_only.sum())
summary["union"] = int(gs.union.sum())
summary["union_unique_symbols"] = int(gs[gs.union].symbol.nunique())
b = gs[gs.susie_pos].max_s4
summary["multi_bands_gt0.9_0.8-0.9_0.5-0.8"] = [int((b > 0.9).sum()), int(((b > 0.8) & (b <= 0.9)).sum()), int(((b > 0.5) & (b <= 0.8)).sum())]
a_ = gs[gs.abf_only].max_a4
summary["abf_only_bands_gt0.9_0.8-0.9_0.5-0.8"] = [int((a_ > 0.9).sum()), int(((a_ > 0.8) & (a_ <= 0.9)).sum()), int(((a_ > 0.5) & (a_ <= 0.8)).sum())]
ao = gs[gs.abf_only]
summary["abf_only_split"] = {"no_finite_susie_and_not_untestable": int((ao.max_s4.isna() & ~ao.untestable).sum()),
                             "finite_susie_le_0.5": int(ao.max_s4.notna().sum()),
                             "untestable_no_finite": int((ao.max_s4.isna() & ao.untestable).sum()),
                             "untestable_any": int(ao.untestable.sum())}
summary["untestable_genes_t12_not_positive"] = sorted(gs[gs.untestable & ~gs.susie_pos].symbol.tolist())
gs.to_csv(f"{OUT}/gene_states.tsv", sep="\t")

# ---------- Panel A: COLOC genes by driving trait (for ribbon check) ----------
def driver(rowset, col):
    return rowset.loc[rowset[col].idxmax()]
pos_rows = g[g.ensembl.isin(gs.index[gs.susie_pos])]
drv_s = pos_rows.loc[pos_rows.groupby("ensembl").s4.idxmax()]
summary["multi_by_driving_trait"] = drv_s.trait.value_counts().to_dict()
summary["multi_by_driving_trait_ancestry"] = {f"{k[0]}|{k[1]}": int(v) for k, v in drv_s.groupby(["trait", "anc5"]).size().items()}
abf_rows = g[g.ensembl.isin(gs.index[gs.abf_only])]
drv_a = abf_rows.loc[abf_rows.groupby("ensembl").a4.idxmax()]
summary["abf_only_by_driving_trait"] = drv_a.trait.value_counts().to_dict()
summary["abf_only_by_driving_trait_ancestry"] = {f"{k[0]}|{k[1]}": int(v) for k, v in drv_a.groupby(["trait", "anc5"]).size().items()}
summary["multi_driven_by_MVP_EUR_strata"] = int(drv_s.gwas_name.isin(["MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR"]).sum())
# genes supported ONLY by MVP EUR strata (SuSiE > 0.5 only there)
sp = g[g.s4 > 0.5]
only_mvp = sp.groupby("ensembl").gwas_name.agg(lambda s: set(s) <= {"MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR"})
summary["multi_supported_only_by_MVP_EUR"] = int(only_mvp.sum())
drv_s[["ensembl", "gene", "gwas_name", "trait", "anc5", "s4", "top_snp"]].to_csv(f"{OUT}/A_multi_drivers.tsv", sep="\t", index=False)
drv_a[["ensembl", "gene", "gwas_name", "trait", "anc5", "a4", "top_snp"]].to_csv(f"{OUT}/A_abfonly_drivers.tsv", sep="\t", index=False)

# ---------- Panel B ----------
fs = g[g.s4.notna()]
bl = fs.loc[fs.groupby("ensembl").s4.idxmax()].copy()
# ties at the max
mx = fs.groupby("ensembl").s4.transform("max")
ties = fs[fs.s4 == mx].groupby("ensembl").size()
summary["B_genes_with_tied_max_rows"] = ties[ties > 1].index.map(lambda e: gs.loc[e, "symbol"]).tolist()
bl[["c", "p"]] = bl.top_snp.str.split(":", expand=True).iloc[:, :2]
bl["p"] = bl.p.astype(int)
fm = pd.read_csv(f"{U}/combined_finemapping_at_leads.tsv.gz", sep="\t", dtype={"chromosome": str})
fm["p"] = fm.position.astype(int)
agg = fm.groupby(["study", "chromosome", "p"]).agg(n_rec=("susie_pip", "size"), susie_pip_max=("susie_pip", "max"),
                                                   susie_pip_min=("susie_pip", "min"),
                                                   susie_pip_clean_max=("susie_pip_clean", "max")).reset_index()
bl = bl.merge(agg, left_on=["gwas_name", "c", "p"], right_on=["study", "chromosome", "p"], how="left")
summary["B_genes_with_finite_susie"] = len(bl)
summary["B_with_fm_record"] = int(bl.n_rec.notna().sum())
summary["B_with_susie_pip"] = int(bl.susie_pip_max.notna().sum())
summary["B_max_from_MVP_EUR"] = int(bl.gwas_name.isin(["MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR"]).sum())
sh = bl[bl.susie_pip_max.notna() & ((bl.s4 > 0.05) | (bl.susie_pip_max >= 0.99))].copy()
summary["B_shown"] = len(sh)
summary["B_corner_pip>=0.9_pp4<0.5"] = int(((sh.susie_pip_max >= 0.9) & (sh.s4 < 0.5)).sum())
summary["B_topright_pip>=0.9_pp4>0.5"] = int(((sh.susie_pip_max >= 0.9) & (sh.s4 > 0.5)).sum())
sh2 = bl[bl.susie_pip_min.notna() & ((bl.s4 > 0.05) | (bl.susie_pip_min >= 0.99))]
summary["B_shown_if_min_pip"] = len(sh2)
summary["B_corner_if_min_pip"] = int(((sh2.susie_pip_min >= 0.9) & (sh2.s4 < 0.5)).sum())
# alternative: all genes with a SuSiE PIP (no display filter)
al = bl[bl.susie_pip_max.notna()]
summary["B_all_with_pip_corner"] = int(((al.susie_pip_max >= 0.9) & (al.s4 < 0.5)).sum())
summary["B_all_with_pip_topright"] = int(((al.susie_pip_max >= 0.9) & (al.s4 > 0.5)).sum())
for gene in ["RORA", "HSD17B13", "PNPLA3", "TM6SF2", "MTARC1", "GCKR", "THBS3", "POLL"]:
    r = bl[bl.gene == gene]
    if len(r):
        r = r.iloc[0]
        summary[f"B_{gene}"] = {"gwas": r.gwas_name, "lead": r.top_snp, "pp4_susie": float(r.s4),
                                "pip_max": None if pd.isna(r.susie_pip_max) else float(r.susie_pip_max),
                                "pip_min": None if pd.isna(r.susie_pip_min) else float(r.susie_pip_min)}
bl.to_csv(f"{OUT}/B_gene_leads.tsv", sep="\t", index=False)
sh.to_csv(f"{OUT}/B_shown.tsv", sep="\t", index=False)

# ---------- Panel E ----------
anc = sp.groupby("ensembl").anc5.agg(lambda s: set(s))
eur_only = anc.map(lambda s: s == {"EUR"})
shared = anc.map(lambda s: "EUR" in s and len(s) > 1)
non = anc.map(lambda s: "EUR" not in s)
summary["E_eur_only"] = int(eur_only.sum()); summary["E_shared"] = int(shared.sum()); summary["E_noneur_only"] = int(non.sum())
# EUR-only genes: tested in any non-EUR study (finite SuSiE)?
fin = g[g.s4.notna()].groupby("ensembl").anc5.agg(lambda s: set(s))
eo = anc[eur_only].index
summary["E_eur_only_with_finite_nonEUR_susie"] = int(sum(len(fin.get(e, set()) - {"EUR"}) > 0 for e in eo))
summary["E_eur_only_without_any_nonEUR_susie"] = int(sum(len(fin.get(e, set()) - {"EUR"}) == 0 for e in eo))
abfnon = g[(g.a4 > 0.5) & (g.anc5 != "EUR")].ensembl.unique()
summary["E_eur_only_with_nonEUR_ABF_gt0.5"] = int(sum(e in set(abfnon) for e in eo))
# non-EUR-only: GWS of lead variant
gp = pd.read_csv(f"{U}/gwas_p_at_leads.tsv.gz", sep="\t", dtype={"chromosome": str})
gp["key"] = gp.study + "|" + gp.chromosome + ":" + gp.position.astype(str)
pmin = gp.groupby("key").pval.min()
nr = sp[sp.ensembl.isin(non[non].index)].copy()
nr["key"] = nr.gwas_name + "|" + nr.top_snp.str.split(":").str[:2].str.join(":")
nr["gwas_p"] = nr.key.map(pmin)
lead_max = nr.loc[nr.groupby("ensembl").s4.idxmax()]
lead_minp = nr.loc[nr.groupby("ensembl").gwas_p.idxmin()]
summary["E_noneur_GWS_at_maxPP4_lead"] = int((lead_max.gwas_p < 5e-8).sum())
summary["E_noneur_GWS_any_positive_lead"] = int((lead_minp.gwas_p < 5e-8).sum())
summary["E_noneur_p_missing"] = int(lead_max.gwas_p.isna().sum())
gws = lead_minp[lead_minp.gwas_p < 5e-8]
summary["E_GWS_genes"] = sorted(gws.gene.tolist())
summary["E_GWS_ancestry_of_min_p_lead"] = gws.anc5.value_counts().to_dict()
gws2 = lead_max[lead_max.gwas_p < 5e-8]
summary["E_GWS_genes_maxPP4_lead"] = sorted(gws2.gene.tolist())
summary["E_GWS_ancestry_maxPP4_lead"] = gws2.anc5.value_counts().to_dict()
lead_max.to_csv(f"{OUT}/E_noneur_only_leads_maxPP4.tsv", sep="\t", index=False)
lead_minp.to_csv(f"{OUT}/E_noneur_only_leads_minp.tsv", sep="\t", index=False)
# nonEUR-only genes: neighbour within 1 Mb of a EUR positive lead
eurpos = sp[sp.anc5 == "EUR"].copy()
eurpos[["c", "p"]] = eurpos.top_snp.str.split(":", expand=True).iloc[:, :2]
eurpos["p"] = eurpos.p.astype(int)
def has_eur_nb(row):
    c, p = row.top_snp.split(":")[:2]
    e = eurpos[eurpos.c == c]
    return bool(((e.p - int(p)).abs() <= 1_000_000).any())
lead_max = lead_max.assign(eur_nb=lead_max.apply(has_eur_nb, axis=1))
summary["E_noneur_only_with_EUR_neighbour_1Mb"] = int(lead_max.eur_nb.sum())
summary["E_GWS_genes_without_EUR_neighbour"] = sorted(lead_max[(lead_max.gwas_p < 5e-8) & ~lead_max.eur_nb].gene.tolist())

# ---------- Panel F ----------
cell = g.groupby(["ensembl", "anc5"]).agg(s4=("s4", "max"), a4=("a4", "max")).reset_index()
cell["supp"] = (cell.s4 > 0.5) | (cell.a4 > 0.5)
br = cell.groupby("ensembl").supp.sum()
sel = br[br >= 4].index
summary["F_n_genes_ge4"] = len(sel)
summary["F_n_genes_all5"] = int((br == 5).sum())
fc = cell[cell.ensembl.isin(sel)].copy()
fc["symbol"] = fc.ensembl.map(gs.symbol)
fc["state"] = np.where(fc.s4 > 0.5, "multi", np.where(fc.a4 > 0.5, "single_only", "le_0.5"))
fc["shown"] = np.where(fc.state == "multi", fc.s4, np.where(fc.state == "single_only", fc.a4, fc[["s4", "a4"]].max(axis=1)))
grid = pd.MultiIndex.from_product([sel, ["EUR", "AFR", "AMR", "EAS", "SAS"]], names=["ensembl", "anc5"])
fc = fc.set_index(["ensembl", "anc5"]).reindex(grid).reset_index()
fc["symbol"] = fc.ensembl.map(gs.symbol)
summary["F_cells"] = len(fc); summary["F_cells_missing"] = int(fc.shown.isna().sum())
summary["F_state_counts"] = fc.state.value_counts(dropna=False).to_dict()
# number of studies per ancestry
summary["F_studies_per_ancestry"] = reg.assign(a=reg.ancestry.replace({"CSA": "SAS"})).groupby("a").size().to_dict()
order = fc.groupby("symbol").agg(breadth=("state", lambda s: (s == "multi").sum()), mean_s4=("s4", "mean")).reset_index()
order = order.sort_values(["breadth", "mean_s4", "symbol"], ascending=[False, False, True])
summary["F_row_order_breadth_meanS4"] = order.symbol.tolist()
fc.to_csv(f"{OUT}/F_cells.tsv", sep="\t", index=False)

# ---------- Panel C lead selection (annotation done in R) ----------
un = gs[gs.union].index
rows = g[g.ensembl.isin(un)]
cl = []
for e in un:
    r = rows[rows.ensembl == e]
    if gs.loc[e, "susie_pos"]:
        rr = r.loc[r.s4.idxmax()]; src = "susie_max"
    else:
        rr = r.loc[r.a4.idxmax()]; src = "abf_max"
    cl.append({"ensembl": e, "symbol": gs.loc[e, "symbol"], "gwas": rr.gwas_name, "top_snp": rr.top_snp, "src": src})
cl = pd.DataFrame(cl)
cl.to_csv(f"{OUT}/C_union_leads.tsv", sep="\t", index=False)
summary["C_union_genes"] = len(cl); summary["C_distinct_leads"] = int(cl.top_snp.nunique())

with open(f"{OUT}/summary.json", "w") as f:
    json.dump(summary, f, indent=1, default=str)
print(json.dumps(summary, indent=1, default=str))
