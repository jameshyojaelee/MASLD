"""Critic check (3C): does any cohort other than GSE126848 show the zero-read artefact signature
(many calls sharing one identical logFC, which is what genes with zero reads in every library of a
cohort produce under voom)? Light proxy for the filterByExpr question; not a substitute for it."""
import pandas as pd
F = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
d = pd.read_csv(F, sep="\t")
d["call"] = (d.FDR < 0.05) & (d.logFC.abs() > 0.5)
d["lfr"] = d.logFC.round(6)
for ds, g in d.groupby("dataset"):
    vc = g.lfr.value_counts()
    top = vc.index[0]; n_top = vc.iloc[0]
    calls_top = int(g[(g.lfr == top) & g.call].shape[0])
    n_shared = int(g[g.lfr.map(vc) >= 10].shape[0]); calls_shared = int(g[(g.lfr.map(vc) >= 10) & g.call].shape[0])
    print(f"{ds}: calls {int(g.call.sum())}; most common logFC {top} shared by {n_top} genes ({calls_top} calls); genes in logFC values shared by >=10 genes: {n_shared} ({calls_shared} calls)")
