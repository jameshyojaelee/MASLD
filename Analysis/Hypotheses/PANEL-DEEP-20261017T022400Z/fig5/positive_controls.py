"""Positive controls for Figure 5.
5B: protein vs mRNA MASLD-vs-control log2FC concordance should be at least as high for abundant hepatocyte proteins
    (a fixed list of canonical secreted / metabolic hepatocyte proteins, chosen before looking at their fold changes)
    and for the top AveExpr quartile of measured proteins, as across the display universe.
5F: Stromal ECM spatial structure: the program should show positive Moran excess before cell-composition adjustment
    (lineage_only design) and most of it should be explained by adjusting for 16 cell-type abundances (all16).
"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
B = P / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
HEP = ["ALB", "APOA1", "APOA2", "APOB", "APOC3", "APOH", "AHSG", "TF", "TTR", "RBP4", "HP", "HPX", "FGA", "FGB", "FGG", "SERPINA1",
       "SERPINC1", "CYP2E1", "CYP3A4", "CYP2C8", "CYP2C9", "CYP1A2", "ARG1", "ASS1", "CPS1", "OTC", "BHMT", "MAT1A", "GNMT", "FABP1",
       "ADH1A", "ADH1B", "ADH4", "ALDOB", "PCK1", "ALDH2", "HMGCS2", "ACADM", "HADHA", "CES1", "GSTA1", "UGT2B7", "SDS", "HGD", "FAH",
       "TAT", "PAH", "HAL", "AGXT", "CRP"]
g = pd.read_csv(B / "recompute/5b_recomputed_genes.tsv", sep="\t")
de = pd.read_csv(P / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv", sep="\t")
g = g.merge(de[["gene", "AveExpr"]], on="gene", how="left")
def summ(x, lab):
    r = stats.spearmanr(x.rna_logFC, x.logFC)
    agree = (np.sign(x.rna_logFC) == np.sign(x.logFC)).mean()
    sig = x[x.padj < 0.05]
    print(f"{lab}: n={len(x)} spearman rho={r.correlation:.3f} (p={r.pvalue:.2g}) sign agreement={agree:.3f}; "
          f"protein-significant n={len(sig)} concordant={int((np.sign(sig.rna_logFC) == np.sign(sig.logFC)).sum())}")
    return {"set": lab, "n": len(x), "rho": r.correlation, "p": r.pvalue, "agree": agree}
res = [summ(g, "5B display universe (all plotted)")]
h = g[g.gene.isin(HEP)]
print("hepatocyte list measured in 5B:", len(h), "of", len(HEP), "missing:", sorted(set(HEP) - set(h.gene)))
res.append(summ(h, "canonical abundant hepatocyte proteins"))
q = g.AveExpr.quantile(0.75)
res.append(summ(g[g.AveExpr >= q], "top AveExpr quartile"))
res.append(summ(g[g.AveExpr < g.AveExpr.quantile(0.25)], "bottom AveExpr quartile"))
print(h.sort_values("rna_logFC")[["gene", "rna_logFC", "rna_padj", "logFC", "padj"]].round(3).to_string())
pd.DataFrame(res).to_csv(B / "recompute/positive_control_5b.tsv", sep="\t", index=False)

f = pd.read_csv(B / "recompute/5f_recomputed_program_level.tsv", sep="\t")
w = f.pivot_table(index=["dataset", "program"], columns="design", values=["excess", "p"]).round(5)
w["fraction_of_excess_removed_by_all16"] = (1 - w[("excess", "all16")] / w[("excess", "lineage_only")]).round(3)
print(w.to_string())
w.to_csv(B / "recompute/positive_control_5f.tsv", sep="\t")
