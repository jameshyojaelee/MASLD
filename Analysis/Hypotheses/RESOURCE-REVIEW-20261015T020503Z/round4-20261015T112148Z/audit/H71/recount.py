# audit-H71: independent recount of every n stated in PRESPEC_H71.md.
# Margins only. Never prints a per-bin overlap, a stratified expectation, a p value
# or any control result. G2 is reported as a count of qualifying bins only.
import pandas as pd, numpy as np
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
G4 = R + "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG10/"
CT = R + "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B = R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
P = R + "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/"
E = R + "data/external/elison2025_sc_chromatin/elison2025_sc_chromatin_calls.csv"

gi = pd.read_parquet(P + "passport_gene_index.parquet")
ev = pd.read_parquet(P + "passport_evidence_long.parquet")
print("catalog genes", len(gi), "unique ensembl", gi.ensembl_id.nunique(),
      "dup symbols", int(gi.symbol.duplicated(keep=False).sum()),
      "symbol_collision", int(gi.symbol_collision.fillna(False).sum()),
      "ensembl has version", bool(gi.ensembl_id.str.contains(r"\.").any()))
print("classes", gi.primary_evidence_class.value_counts().to_dict())
tr = ev[ev.evidence_domain == "transcriptomics"].set_index("ensembl_id")
ge = ev[ev.evidence_domain == "genetics"].set_index("ensembl_id")
cls = gi.set_index("ensembl_id").primary_evidence_class
x = pd.DataFrame({"cls": cls.astype(str), "tr": tr.call_state.astype(str), "trt": tr.testability_state.astype(str),
                  "ge": ge.call_state.astype(str), "get": ge.testability_state.astype(str)})
print("class x (tr call, ge call)")
print(x.groupby(["cls", "tr", "ge"]).size().to_string())
print("untested & tr supported & ge untestable", int(((x.cls == "untested") & (x.tr == "supported") & (x.ge == "untestable")).sum()))
print("untested & tr supported (any ge)", int(((x.cls == "untested") & (x.tr == "supported")).sum()))
print("concordant & tr supported", int(((x.cls == "concordant") & (x.tr == "supported")).sum()))

# bulk table and canonical rule
d = pd.read_csv(B)
d["g"] = d.gene.astype(str).str.split(".").str[0]
print("bulk rows", len(d), "unique g", d.g.nunique())
d = d.drop_duplicates("g")
d["deg"] = (d.padj < 0.05) & (d.logFC.abs() > 0.5)
print("canonical DEGs (padj<0.05 & |logFC|>0.5)", int(d.deg.sum()))
# transcriptomics supported vs canonical rule
trs = set(tr.index[tr.call_state == "supported"])
canon = set(d.g[d.deg])
print("tr supported", len(trs), "canonical", len(canon), "both", len(trs & canon),
      "tr-only", len(trs - canon), "canon-only", len(canon - trs))
# compare estimate with logFC
t2 = tr[["estimate", "q_value"]].join(d.set_index("g")[["logFC", "padj"]], how="inner")
print("tr rows joined to bulk", len(t2), "max |estimate-logFC|", float((t2.estimate - t2.logFC).abs().max()),
      "max |q-padj|", float((t2.q_value - t2.padj).abs().max()))

# COLOC tier 1/2
c = pd.read_csv(CT)
c["g"] = c.ensembl.astype(str).str.split(".").str[0]
print("coloc rows", len(c), "unique g", c.g.nunique())
c = c.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
c["susie"] = pd.to_numeric(c.coloc_best_susie_pp4, errors="coerce")
print("coloc states", c.susie_state_t12.value_counts().to_dict())
print("positive vs susie>=0.5 crosstab")
print(pd.crosstab(c.susie_state_t12, c.susie >= 0.5).to_string())
print("susie_driving_gwas non-null among susie>=0.5", int(c.loc[c.susie >= 0.5, "susie_driving_gwas"].notna().sum()))
jt = d.merge(c[["g", "susie", "susie_state_t12", "susie_driving_gwas"]], on="g", how="inner")
ue = jt[jt.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
print("U_full", len(jt), "G", int((jt.susie >= 0.5).sum()), "states", jt.susie_state_t12.value_counts().to_dict())
print("U_eval", len(ue), "G", int((ue.susie >= 0.5).sum()), "positive", int((ue.susie_state_t12 == "positive").sum()))
print("DEGs in U_full", int(jt.deg.sum()), "DEGs in U_eval", int(ue.deg.sum()))
print("catalog genes in bulk table", int(gi.ensembl_id.isin(d.g).sum()))
print("catalog genetics supported", int((ev[ev.evidence_domain == 'genetics'].call_state == 'supported').sum()),
      "of which in U_full", int(ge.index[ge.call_state == 'supported'].isin(jt.g).sum()))

# Kumak lists
lists = {}
for s in ["NAS", "Fibrosis"]:
    k = pd.read_csv(G4 + f"kumak_{s}.tsv", sep="\t")
    k["dir"] = np.where(k.total_upregulated > k.total_downregulated, "up",
                        np.where(k.total_upregulated < k.total_downregulated, "down", "tie"))
    m = k.merge(gi[["symbol", "ensembl_id", "primary_evidence_class"]], left_on="Gene", right_on="symbol", how="left")
    print(f"\n{s} rows {len(k)} unique {k.Gene.nunique()} dir {k.dir.value_counts().to_dict()}")
    print(s, "GTEX", k.GTEX.value_counts().to_dict(), "pubmed", k.pubmed_mentioned.value_counts().to_dict())
    print(s, "total min (0.95 quantile floor)", int(k.total.min()), "max", int(k.total.max()))
    print(s, "matched", int(m.ensembl_id.notna().sum()), "unmatched", sorted(m.loc[m.ensembl_id.isna(), "Gene"]))
    print(s, "classes", m.primary_evidence_class.value_counts().to_dict(),
          "anchored+concordant", int(m.primary_evidence_class.isin(["genetically_anchored", "concordant"]).sum()))
    mm = m.dropna(subset=["ensembl_id"])
    L = set(mm.ensembl_id)
    lists[s] = L
    mb = mm.merge(d[["g", "logFC", "padj", "deg"]], left_on="ensembl_id", right_on="g", how="inner")
    ag = ((mb.dir == "up") & (mb.logFC > 0)) | ((mb.dir == "down") & (mb.logFC < 0))
    sg = mb.padj < 0.05
    print(s, "in bulk table", len(mb), "sign agree", int(ag.sum()), "| q<0.05", int(sg.sum()), "agree", int((ag & sg).sum()))
    print(s, "bulk-matched logFC>0", int((mb.logFC > 0).sum()), "<0", int((mb.logFC < 0).sum()), "==0", int((mb.logFC == 0).sum()))
    inf = jt[jt.g.isin(L)]
    print(s, "in U_full", len(inf), "states", inf.susie_state_t12.value_counts().to_dict(),
          "in U_eval", int(ue.g.isin(L).sum()))
    print(s, "bulk-supported in matched (tr supported)", int(mm.ensembl_id.isin(trs).sum()))
    for gene in ["COL1A1", "COL3A1", "THBS2", "FBLN5", "PDGFA", "LAMA2", "LAMA3", "GNMT", "MAT1A", "CYP2C19"]:
        r = m[m.Gene == gene]
        if len(r):
            print("  ", gene, r.dir.iloc[0], r.primary_evidence_class.iloc[0])
print("\nNAS & Fib overlap", len(lists["NAS"] & lists["Fibrosis"]), "union", len(lists["NAS"] | lists["Fibrosis"]))
print("U_eval overlap of the two lists", int(ue.g.isin(lists["NAS"] & lists["Fibrosis"]).sum()))

# G2 feasibility: count of qualifying bins only (no per-bin vectors)
ue = ue.sort_values("g").reset_index(drop=True)
r = ue.AveExpr.rank(method="first")
for nb in (10, 5):
    ue["bin"] = pd.qcut(r, nb, labels=False)
    gk = ue.groupby("bin").apply(lambda t: int((t.susie >= 0.5).sum()))
    for s, L in lists.items():
        lk = ue.groupby("bin").apply(lambda t: int(t.g.isin(L).sum()))
        ok = int(((gk >= 5) & (lk >= 5)).sum())
        print(f"G2 bins={nb} list={s}: bins meeting >=5 genetic and >=5 list = {ok} of {nb}")
ue["bin"] = pd.qcut(r, 10, labels=False)
print("bin sizes", ue.bin.value_counts().sort_index().tolist())
print("AveExpr ties in U_eval", int(ue.AveExpr.duplicated().sum()))

el = pd.read_csv(E)
m = el.merge(gi[["symbol", "primary_evidence_class"]], left_on="gene", right_on="symbol", how="left")
print("\nElison", len(el), "matched", int(m.primary_evidence_class.notna().sum()), m.primary_evidence_class.value_counts().to_dict())
print(m.loc[m.gene.isin(["PNPLA3", "TM6SF2", "HSD17B13", "MBOAT7"]), ["gene", "primary_evidence_class"]].to_string(index=False))
