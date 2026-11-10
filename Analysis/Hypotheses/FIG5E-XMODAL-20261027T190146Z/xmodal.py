"""Exploratory cross-modal gene table for the Figure 5E redesign (not an analysis result)."""
import pandas as pd, numpy as np, sys
OUT = sys.argv[1]
R = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.tsv.gz"
g = pd.read_csv(R, sep="\t", low_memory=False)
g = g[["ensembl_id","symbol","primary_evidence_class","r8_genetics_call_state","r8_rna_call_state","rna_log2fc","rna_bh_q",
       "rna_canonical_call","coding_cs_state","coding_cs_top_consequence","caqtl_gene_link","zon_mrna_call_yakubovsky2026",
       "zon_protein_call_weiss2026","sc_best_lineage","liver_proteome_any","gtex_postmortem_fat_state","fig3e_F4_vs_F0_call"]]
# protein DE, 46 MASLD vs 12 controls
p = pd.read_csv("figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/protein_de_adjusted.tsv", sep="\t")
print("protein cols", p.columns.tolist())
gcol = [c for c in p.columns if c.lower() in ("gene","genes","symbol")][0]
p = p.rename(columns={gcol:"symbol","logFC":"prot_log2fc","padj":"prot_padj"})[["symbol","prot_log2fc","prot_padj"]]
p = p.sort_values("prot_padj").drop_duplicates("symbol")
g = g.merge(p, on="symbol", how="left")
# COLOC signal in lineage open chromatin (5D)
pairs = pd.read_csv("Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv", sep="\t")
op = pairs[pairs.evidence_state.isin(["replicated_accessible","source_dependent"])]
oa = op.groupby("gene").agg(open_lineages=("lineage", lambda x: ";".join(sorted(set(x)))),
                            open_traits=("trait_class", lambda x: ";".join(sorted(set(x))))).reset_index().rename(columns={"gene":"symbol"})
g = g.merge(oa, on="symbol", how="left")
# disease-shifted peaks linked to gene (descriptive link; peak significance not transferred to genes)
A = "figures/main/fig5_molecular_context/data/atac_context_20261026/"
use = ["cohort","lineage","gene_symbol","peak_id","method","distance_to_tss_bp","logFC","qvalue","lodo_all_same_sign","RNA_logFC","RNA_qvalue"]
pg = pd.read_csv(A+"peak_to_gene_context.tsv.gz", sep="\t", usecols=use)
sig = pg[pg.qvalue < 0.05].copy()
sig["promoter"] = sig.distance_to_tss_bp.abs() <= 2000
da = sig.groupby("gene_symbol").agg(
    da_n_peaks=("peak_id","nunique"), da_n_promoter=("promoter","sum"),
    da_ctx=("cohort", lambda x: ";".join(sorted(set(x)))),
    da_lineage=("lineage", lambda x: ";".join(sorted(set(x)))),
    da_up=("logFC", lambda x: int((x>0).sum())), da_down=("logFC", lambda x: int((x<0).sum()))).reset_index().rename(columns={"gene_symbol":"symbol"})
g = g.merge(da, on="symbol", how="left")
# snRNA hepatocyte disease effects (GSE244832)
r = pd.read_csv(A+"rna_models_rna_effects.tsv.gz", sep="\t")
print("rna cols", r.columns.tolist()[:20])
g.to_csv(OUT+"/xmodal_gene_table.tsv.gz", sep="\t", index=False)
f = pd.DataFrame({"symbol": g.symbol,
  "G": g.r8_genetics_call_state.eq("supported"),
  "B": g.r8_rna_call_state.eq("supported"),
  "P": g.prot_padj.lt(0.05),
  "A_open": g.open_lineages.notna(),
  "A_da": g.da_n_peaks.fillna(0).gt(0),
  "Z": g.zon_mrna_call_yakubovsky2026.isin(["central","portal"])})
for c in "G B P A_open A_da Z".split(): print(c, int(f[c].sum()))
f["n_new"] = f[["P","A_open","A_da"]].sum(1)
f["n_all"] = f[["G","B","P","A_open","A_da"]].sum(1)
print(f.groupby(["G","B","P","A_open","A_da"]).size().sort_values(ascending=False).head(25).to_string())
f.to_csv(OUT+"/xmodal_flags.tsv", sep="\t", index=False)
top = g.merge(f, on="symbol")
top = top[(top.n_all >= 3)].sort_values(["n_all","prot_padj"], ascending=[False, True])
cols = ["symbol","n_all","G","B","P","A_open","A_da","rna_log2fc","prot_log2fc","prot_padj","open_lineages","open_traits","da_n_peaks","da_n_promoter","da_ctx","da_lineage","da_up","da_down","coding_cs_state","caqtl_gene_link","zon_mrna_call_yakubovsky2026","sc_best_lineage"]
top[cols].to_csv(OUT+"/xmodal_ge3.tsv", sep="\t", index=False)
print(len(top)); print(top[cols].head(80).to_string())
for s in ["GNMT","MAT1A","CYP2C19","IGFBP7","BICC1","RORA","PNPLA3","TM6SF2","HSD17B13","GCKR","PPP1R3B","LITAF"]:
    print(top[top.symbol==s][cols].to_string(header=False) if (top.symbol==s).any() else s+" <3")
