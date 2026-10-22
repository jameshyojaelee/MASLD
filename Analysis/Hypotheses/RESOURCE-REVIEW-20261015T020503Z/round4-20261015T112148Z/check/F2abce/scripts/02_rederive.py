"""Blind check F2abce, step 2: re-derive Fig 2B same-row values, Fig 2C membership
and leads, Fig 2E locus context, from COLOC r2 + fine-mapping release tables."""
import os, sys, json
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CK = sys.argv[1]
OUT = os.path.join(CK, "out")
res = {}

d = pd.read_csv(os.path.join(OUT, "coloc_main.tsv.gz"), sep="\t",
                dtype={"gene": str, "top_snp": str, "gwas_name": str})
d = d.sort_values("row").reset_index(drop=True)
d["pp4s"] = pd.to_numeric(d["PP.H4.susie"], errors="coerce")
d["pp4a"] = pd.to_numeric(d["PP.H4.abf"], errors="coerce")
has_snp = d.top_snp.notna() & (d.top_snp != "")
d["snp_chr"] = np.nan; d["snp_pos"] = np.nan
sp = d.loc[has_snp, "top_snp"].str.split(":", expand=True)
d.loc[has_snp, "snp_chr"] = sp[0].astype(str)
d.loc[has_snp, "snp_pos"] = sp[1].astype(str)
d["key"] = d.gwas_name + "," + d.snp_chr.astype(str) + "," + d.snp_pos.astype(str)
d.loc[~has_snp, "key"] = np.nan
res["main_rows"] = int(len(d)); res["main_rows_without_top_snp"] = int((~has_snp).sum())

# ---- R16-08 cross-check from pandas
s = d[(d.gwas_name == "2023_36280732_NAFLD_UKBB_EUR") & (d.chr.astype(str) == "1")]
res["R16_08_chr1_rows"] = int(len(s)); res["R16_08_chr1_susie_nonNA"] = int(s.pp4s.notna().sum())
res["R16_08_chr1_susie_gt0.5"] = int((s.pp4s > 0.5).sum())

# ---- fine-mapping rows at COLOC keys
fc = pd.read_csv(os.path.join(OUT, "fm_canon_filtered.csv"), dtype={"study": str, "susie_reliable": str},
                 usecols=["chromosome", "position", "allele1", "allele2", "study", "ancestry", "susie_pip",
                          "carma_pip", "susie_reliable", "recommended_pip", "susiex_pip", "mesusie_pip"])
fc["src"] = "canonical"
ft = pd.read_csv(os.path.join(OUT, "fm_targ_filtered.csv"), dtype={"study": str, "susie_reliable": str},
                 usecols=["chromosome", "position", "allele1", "allele2", "study", "ancestry", "susie_pip",
                          "carma_pip", "susie_reliable", "recommended_pip"])
ft["src"] = "targeted"
fm = pd.concat([fc, ft], ignore_index=True)
for c in ["susie_pip", "carma_pip", "recommended_pip", "susiex_pip", "mesusie_pip"]:
    fm[c] = pd.to_numeric(fm[c], errors="coerce")
fm["key"] = fm.study + "," + fm.chromosome.astype(str) + "," + fm.position.astype(str)
# producer's x: pmax(recommended_pip, susie_pip, na.rm)
fm["causal"] = fm[["recommended_pip", "susie_pip"]].max(axis=1, skipna=True)
# within-ancestry rule (aggregator lines 175-184, before joint PIPs are merged in)
rel = fm.susie_reliable.astype(str).str.upper().eq("TRUE")
eur = fm.ancestry.astype(str).eq("EUR")
fm["wa"] = np.where(eur & rel, fm[["susie_pip", "carma_pip"]].max(axis=1, skipna=True), fm.carma_pip)
res["fm_rows_canon"] = int((fm.src == "canonical").sum()); res["fm_rows_targ"] = int((fm.src == "targeted").sum())
g = fm.groupby("key")
kx = pd.DataFrame({
    "x_susie_max": g.susie_pip.max(),
    "x_susie_min": g.susie_pip.min(),
    "n_rows": g.size(),
    "n_alleles": g.apply(lambda t: t[["allele1", "allele2"]].drop_duplicates().shape[0]),
    "x_wa_max": g.wa.max(),
    "x_causal_max": g.causal.max(),
})
cs = fm[fm.src == "canonical"].groupby("key").susie_pip.max()
kx["x_susie_canon_first"] = cs.reindex(kx.index)
ts = fm[fm.src == "targeted"].groupby("key").susie_pip.max()
kx["x_susie_canon_first"] = kx.x_susie_canon_first.fillna(ts.reindex(kx.index))
kx.to_csv(os.path.join(OUT, "key_x.tsv.gz"), sep="\t")

roster = ["HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA"]

def gate(df, xcol):
    pp = df.pp4.clip(lower=0)
    return df[xcol].notna() & ((pp > 0.05) | (df[xcol] >= 0.99) | df.gene.isin(roster))

# ---- installed replication (producer rule, gene-level max x over any matched lead)
dy = d[has_snp & d.pp4s.notna()]
idx = dy.groupby("gene", sort=False).pp4s.idxmax()
y = dy.loc[idx, ["gene", "gwas_name", "top_snp", "key", "pp4s"]].rename(columns={"pp4s": "pp4"}).set_index("gene")
res["genes_finite_pp4_susie"] = int(len(y))
res["genes_finite_pp4_susie_incl_rows_without_snp"] = int(d[d.pp4s.notna()].gene.nunique())
# ties at the max with a different lead
mx = dy.groupby("gene").pp4s.transform("max")
tie = dy[dy.pp4s == mx].groupby("gene").apply(lambda t: t[["gwas_name", "top_snp"]].drop_duplicates().shape[0])
y["n_max_leads"] = tie.reindex(y.index)
res["genes_with_tied_max_pp4_leads"] = int((y.n_max_leads > 1).sum())

m = d[has_snp].merge(kx[["x_causal_max"]], left_on="key", right_index=True, how="inner")
gx = m.groupby("gene").x_causal_max.max()
inst = y.copy(); inst["x"] = gx.reindex(inst.index); inst = inst.reset_index()
inst = inst[inst.x.notna()]
pl = inst[gate(inst, "x")]
res["installed_rule_rows"] = int(len(inst)); res["installed_rule_plotted"] = int(len(pl))
res["installed_rule_corner"] = int(((pl.x >= 0.9) & (pl.pp4 < 0.5)).sum())
res["installed_rule_topright"] = int(((pl.x >= 0.9) & (pl.pp4 > 0.5)).sum())

# ---- same-row x
sr = y.reset_index().merge(kx, left_on="key", right_index=True, how="left")
sr["fm_record"] = sr.n_rows.notna()
res["same_row_genes_with_fm_record"] = int(sr.fm_record.sum())
for xcol in ["x_susie_max", "x_susie_canon_first", "x_susie_min", "x_wa_max", "x_causal_max"]:
    t = sr.copy(); t["x"] = t[xcol]
    p = t[gate(t, "x")]
    res[f"same_row[{xcol}]"] = {
        "genes_with_finite_x": int(t.x.notna().sum()),
        "plotted": int(len(p)),
        "corner_pip_ge_0.9_pp4_lt_0.5": int(((p.x >= 0.9) & (p.pp4 < 0.5)).sum()),
        "topright_pip_ge_0.9_pp4_gt_0.5": int(((p.x >= 0.9) & (p.pp4 > 0.5)).sum()),
        "roster_plotted": sorted(set(p.gene) & set(roster)),
    }
res["keys_with_multiple_distinct_susie_pip"] = int(((kx.x_susie_max - kx.x_susie_min).abs() > 1e-12).sum())
sr["is_plotted_susie"] = gate(sr, "x_susie_max").values
sr["is_plotted_wa"] = gate(sr, "x_wa_max").values
res["plotted_susie_genes_whose_lead_has_gt1_allele_pair"] = int(sr[sr.is_plotted_susie].n_alleles.gt(1).sum())
for gname in roster:
    r = sr[sr.gene == gname]
    res[f"roster_{gname}"] = r[["gwas_name", "top_snp", "pp4", "x_susie_max", "x_wa_max", "x_causal_max",
                                 "n_max_leads"]].to_dict("records")
sr.to_csv(os.path.join(OUT, "fig2b_same_row_rederived.tsv"), sep="\t", index=False)

# ---- Fig 2C membership (max of SuSiE and ABF > 0.5; lead = row with max, ties by gwas_name)
c = d.copy()
c["pp4_best"] = c[["pp4s", "pp4a"]].max(axis=1, skipna=True)
res["fig2c_rows_both_pp4_missing"] = int(c.pp4_best.isna().sum())
c2 = c[c.pp4_best.notna()].sort_values(["gene", "pp4_best", "gwas_name"], ascending=[True, False, True])
best = c2.groupby("gene", sort=False).head(1)
best = best[best.pp4_best > 0.5]
res["fig2c_union_genes"] = int(len(best))
# NA-first ordering check (data.table setorder default na.last=FALSE)
na_first_genes = set(c[c.pp4_best.isna()].gene)
res["fig2c_genes_with_any_NA_pp4_best_row"] = int(len(na_first_genes & set(best.gene)))
best[["gene", "gwas_name", "top_snp", "pp4_best"]].to_csv(os.path.join(OUT, "fig2c_leads.tsv"), sep="\t", index=False)

# ---- variants to annotate: Fig 2C leads + Fig 2B plotted same-row leads
va = pd.concat([best[["gene", "top_snp"]].assign(panel="2C"),
                sr[sr.is_plotted_susie][["gene", "top_snp"]].assign(panel="2B_susie"),
                sr[sr.is_plotted_wa][["gene", "top_snp"]].assign(panel="2B_wa")])
va.to_csv(os.path.join(OUT, "variants_to_annotate.tsv"), sep="\t", index=False)

# ---- Fig 2E
reg = pd.read_csv(f"{ROOT}/GWAS/finemapping/config/gwas_registry.tsv", sep="\t")
anc = dict(zip(reg.study_name, reg.ancestry))
e = d.copy(); e["anc"] = e.gwas_name.map(anc)
sig = e[(e.pp4s > 0.5) & e.anc.notna()]
ancset = sig.groupby("gene").anc.apply(set)
part = ancset.apply(lambda a: "EUR only" if a == {"EUR"} else ("shared" if "EUR" in a else "nonEUR"))
res["fig2e_partition"] = part.value_counts().to_dict()
memb = pd.read_csv(f"{ROOT}/figures/main/fig2_genetics/panels/Fig2E_ancestry_unique_coloc_GWS_membership.tsv", sep="\t")
nonE = memb[memb.display_category.str.startswith("non_EUR")].copy()
res["fig2e_nonEUR_set_matches_membership"] = set(nonE.gene) == set(part[part == "nonEUR"].index)
eursig = sig[(sig.anc == "EUR") & sig.snp_pos.notna()][["gene", "snp_chr", "snp_pos"]].copy()
eursig["snp_pos"] = eursig.snp_pos.astype(int)

def neighbours(chrom, pos, W=1_000_000):
    t = eursig[(eursig.snp_chr == str(chrom)) & ((eursig.snp_pos - int(pos)).abs() <= W)]
    return sorted(set(t.gene))

rows = []
for _, r in nonE.iterrows():
    ch, ps = r.best_snp.split(":")[:2]
    nb_minp = neighbours(ch, ps)
    gs = sig[(sig.gene == r.gene) & (sig.anc != "EUR")]
    top = gs.loc[gs.pp4s.idxmax()]
    nb_maxpp = neighbours(top.snp_chr, top.snp_pos)
    nb_any = sorted(set().union(*[neighbours(a, b) for a, b in zip(gs.snp_chr, gs.snp_pos)]))
    rows.append(dict(gene=r.gene, cat=r.display_category, best_gwas=r.best_gwas, minp_lead=r.best_snp,
                     maxpp4_gwas=top.gwas_name, maxpp4_lead=top.top_snp,
                     nb_minp_lead=",".join(nb_minp), n_minp=len(nb_minp),
                     nb_maxpp4_lead=",".join(nb_maxpp), n_maxpp4=len(nb_maxpp),
                     nb_any_lead=",".join(nb_any), n_any=len(nb_any)))
E = pd.DataFrame(rows)
E.to_csv(os.path.join(OUT, "fig2e_neighbours.tsv"), sep="\t", index=False)
res["fig2e_n_nonEUR"] = int(len(E))
for col in ["n_minp", "n_maxpp4", "n_any"]:
    res[f"fig2e_genes_with_EUR_neighbour[{col}]"] = int((E[col] > 0).sum())
G = E[E.cat == "non_EUR_GWS_unique"].copy()
G["c"] = G.minp_lead.str.split(":").str[0]; G["p"] = G.minp_lead.str.split(":").str[1].astype(int)
G = G.sort_values(["c", "p"])
loc = []; lid = 0; prev = None
for _, r in G.iterrows():
    if prev is None or r.c != prev[0] or r.p - prev[1] > 1_000_000:
        lid += 1
    loc.append(lid); prev = (r.c, r.p)
G["locus"] = loc
res["fig2e_gws_loci_1Mb_single_linkage"] = int(G.locus.nunique())
nn = G[G.n_minp == 0]
res["fig2e_gws_without_EUR_neighbour"] = {"genes": sorted(nn.gene), "n_genes": int(len(nn)), "n_loci": int(nn.locus.nunique())}
res["fig2e_gws_table"] = G[["gene", "minp_lead", "maxpp4_lead", "locus", "nb_minp_lead", "nb_maxpp4_lead"]].to_dict("records")

with open(os.path.join(OUT, "rederive_summary.json"), "w") as f:
    json.dump(res, f, indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
