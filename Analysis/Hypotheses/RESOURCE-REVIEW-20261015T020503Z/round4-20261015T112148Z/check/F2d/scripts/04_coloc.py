"""Independent COLOC counts: driving tier, genes per tier, candidate definitions
of eligible signal pairs, and Fig 2A fine-mapped-window support."""
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out"
tier = pd.read_csv(ROOT + "/GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
main = tier[tier.placement == "main"]
print("MAIN strata", len(main))
cols = ["gwas_name", "gene", "ensembl", "chr", "PP.H4.susie", "PP.H4.abf", "n_cs_pairs", "method",
        "top_snp", "ancestry", "n_cs_pairs_eligible", "eligibility_source"]
c = pd.read_csv(ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
                usecols=cols, low_memory=False)
print("rows", len(c), "strata", c.gwas_name.nunique())
c = c[c.gwas_name.isin(main.study_name)].copy()
print("MAIN rows", len(c), "MAIN strata in coloc", c.gwas_name.nunique())
c = c.merge(main[["study_name", "trait", "tier"]], left_on="gwas_name", right_on="study_name")
c = c[c.gene.notna() & (c.gene.astype(str) != "")]
print("genes tested", c.gene.nunique())
print("method counts\n", c.method.value_counts(dropna=False).to_string())
print("eligibility_source\n", c.eligibility_source.value_counts(dropna=False).to_string())
pp = c["PP.H4.susie"]
pos = c[pp > 0.5]
genes = pos.gene.unique()
print("multi-signal genes (max PP.H4.susie>0.5)", len(genes), "positive gene-study rows", len(pos))
best = c[c.gene.isin(genes)].sort_values("PP.H4.susie", ascending=False).drop_duplicates("gene")
print("driving tier", best.tier.value_counts().to_dict(), "driving trait", best.trait.value_counts().to_dict())
g1 = set(pos[pos.tier == 1].gene); g2 = set(pos[pos.tier == 2].gene)
print("genes tier1", len(g1), "tier2", len(g2), "both", len(g1 & g2), "tier1 only", len(g1 - g2))

# candidate definitions of eligible pairs per tier
for t in (1, 2):
    ct = c[c.tier == t]
    has = ct["PP.H4.susie"].notna()
    print(f"TIER {t}: rows {len(ct)}; rows with susie PP.H4 {int(has.sum())}; hits(rows) {int((ct['PP.H4.susie']>0.5).sum())};"
          f" studies with susie rows {ct[has].gwas_name.nunique()};"
          f" sum n_cs_pairs {ct.n_cs_pairs.sum()}; sum n_cs_pairs_eligible {ct.n_cs_pairs_eligible.sum()};"
          f" rows n_cs_pairs_eligible>0 {int((ct.n_cs_pairs_eligible>0).sum())}")
    print("   per-study susie rows / hits:", ct[has].groupby("gwas_name").apply(lambda z: (len(z), int((z['PP.H4.susie']>0.5).sum()))).to_dict())

# Fig 2A support
cs = pd.read_csv(ROOT + "/GWAS/finemapping/results/credible_sets.csv",
                 usecols=["chromosome", "position", "locus", "study"], dtype={"locus": str})
cs = cs[cs.study.isin(main.study_name)]
print("MAIN strata with credible sets", cs.study.nunique(), "missing:", sorted(set(main.study_name) - set(cs.study)))
win = cs.groupby(["study", "chromosome", "locus"]).position.agg(["min", "max"]).reset_index()
print("windows", len(win), "max span", (win["max"] - win["min"]).max())
print("UKBB windows", win[win.study.str.startswith("UKBB_")].groupby("study").size().to_dict())
pos = pos.copy()
pos["tchr"] = pos.top_snp.astype(str).str.split(":").str[0]
pos["tpos"] = pd.to_numeric(pos.top_snp.astype(str).str.split(":").str[1], errors="coerce")
print("positive rows with unparsable top_snp", int(pos.tpos.isna().sum()))
wb = {k: v for k, v in win.groupby(["study", "chromosome"])}
wchr = {k: v for k, v in win.groupby("chromosome")}

def inwin(r, pad, same):
    key = (r.gwas_name, int(r.tchr)) if str(r.tchr).isdigit() else None
    if same:
        w = wb.get(key)
    else:
        w = wchr.get(int(r.tchr)) if str(r.tchr).isdigit() else None
    if w is None:
        return False
    return bool(((w["min"] - pad <= r.tpos) & (r.tpos <= w["max"] + pad)).any())

for pad in (0, 250_000):
    pos[f"same{pad}"] = [inwin(r, pad, True) for r in pos.itertuples()]
    pos[f"any{pad}"] = [inwin(r, pad, False) for r in pos.itertuples()]
    print(f"pad {pad}: genes same-study window {pos[pos[f'same{pad}']].gene.nunique()}; any window {pos[pos[f'any{pad}']].gene.nunique()}")
nocs = set(main.study_name) - set(cs.study)
gs = pos.groupby("gene").gwas_name.agg(lambda x: set(x) <= nocs)
print("genes whose positive pairs all in strata without CS", int(gs.sum()))
# per cell
anc = c.drop_duplicates("gwas_name").set_index("gwas_name").ancestry
cellrows = []
for (tr, an), sub in main.assign(anc=main.study_name.map(anc)).groupby(["trait", "anc"]):
    st = set(sub.study_name)
    p_ = pos[pos.gwas_name.isin(st)]
    cellrows.append(dict(cell=f"{tr}|{an}", finemapped_windows=int(win.study.isin(st).sum()), strata=len(st),
                         genes_positive=p_.gene.nunique(), genes_same250=p_[p_["same250000"]].gene.nunique()))
print(pd.DataFrame(cellrows).to_string())
