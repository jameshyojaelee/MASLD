"""Critic re-derivation: 3G inset/3-lncRNA counts (1,157 protein-coding and 549 lncRNA calls over 21,501 genes)
and the 3C five-cohort union under D3 (12,180 / 3,839 / 1,211 mixed; without GSE126848 2,565 / 346)."""
import pandas as pd
R = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
deg = pd.read_csv(R + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
gm = pd.read_csv(R + "data/gencode_v49_gene_metadata.tsv.gz", sep="\t")
deg["ens"] = deg.gene.str.split(".").str[0]
deg = deg.merge(gm[["ensembl_base", "gene_biotype"]], left_on="ens", right_on="ensembl_base", how="left")
deg["call"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
sub = deg[deg.gene_biotype.isin(["protein_coding", "lncRNA"])]
print("pc+lnc genes:", len(sub), sub.gene_biotype.value_counts().to_dict())
print("calls by biotype:", sub[sub.call].gene_biotype.value_counts().to_dict(), "up:", sub[sub.call & (sub.logFC > 0)].gene_biotype.value_counts().to_dict())
# 3C union
C = pd.read_csv(R + "figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv", sep="\t",
                usecols=["gene_id_versioned", "dataset", "logFC", "FDR"])
C = C[C.dataset != "GSE126848"]
C["call"] = (C.FDR < 0.05) & (C.logFC.abs() > 0.5)
G = pd.read_csv(R + "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/gse126848_filterByExpr_refit_genes.tsv.gz", sep="\t")
G["call"] = (G.FDR < 0.05) & (G.logFC.abs() > 0.5); G["dataset"] = "GSE126848"
print("GSE126848 refit calls (FDR<0.05 & |logFC|>0.5):", int(G.call.sum()), "; column 'call' sum:", int(G["call"].sum()), "; genes", len(G))
def union(df):
    c = df[df.call]
    g = c.groupby("gene_id_versioned").agg(k=("dataset", "nunique"), up=("logFC", lambda v: (v > 0).sum()), dn=("logFC", lambda v: (v < 0).sum()))
    multi = g[g.k >= 2]; mixed = multi[(multi.up > 0) & (multi.dn > 0)]
    return len(c), len(g), len(multi), len(mixed), round(100 * len(mixed) / len(multi), 2)
print("five cohorts (D3):", union(pd.concat([C, G[["gene_id_versioned", "dataset", "logFC", "FDR", "call"]]])))
print("four cohorts:", union(C))
