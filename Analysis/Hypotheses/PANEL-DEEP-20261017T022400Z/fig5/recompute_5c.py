"""Independent recomputation of Figure 5C values from upstream analysis outputs.

Upstream (analysis outputs, not the panel's own source CSVs):
  proteomics/panel4c_adjusted_abundance.tsv, panel4c_metadata.tsv, panel4c_process_enrichment.tsv,
  protein_de_adjusted.tsv (774-unit candidate proteomics dir; identical to canonical for these files)
  RNA 774-unit deg_results.csv; gene-id map panel4c_bulk_gene_ids.tsv; fixed row order fixed_panel4c_rows.tsv
"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
C = BASE / "figures/candidates/bulk-context-774-20261014T140904Z_22236311"
P = C / "proteomics"
OUT = BASE / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute"

rows = pd.read_csv(P / "fixed_panel4c_rows.tsv", sep="\t")
ids = pd.read_csv(C / "panel4c_bulk_gene_ids.tsv", sep="\t")
prot = pd.read_csv(P / "protein_de_adjusted.tsv", sep="\t").set_index("gene")
rna = pd.read_csv(BASE / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv").set_index("gene")

# ---- right column: mRNA and protein log2FC
r = rows.merge(ids, on="gene")
r["mrna_log2fc"] = rna.loc[r.bulk_gene_id, "logFC"].values
r["mrna_padj"] = rna.loc[r.bulk_gene_id, "padj"].values
r["mrna_filled"] = (r.mrna_padj < 0.05) & (r.mrna_log2fc.abs() > 0.5)
r["protein_log2fc"] = prot.loc[r.gene, "logFC"].values
r["protein_padj"] = prot.loc[r.gene, "padj"].values
r["protein_filled"] = r.protein_padj < 0.05
r["sign_agree"] = np.sign(r.mrna_log2fc) == np.sign(r.protein_log2fc)
print(r[["gene", "mrna_log2fc", "mrna_padj", "mrna_filled", "protein_log2fc", "protein_padj", "protein_filled"]].to_string())
print("mRNA filled", r.mrna_filled.sum(), list(r.gene[r.mrna_filled]))
print("protein filled", r.protein_filled.sum(), "open:", list(r.gene[~r.protein_filled]))
print("sign agreement", r.sign_agree.sum(), "of", len(r))
r.to_csv(OUT / "5c_recomputed_right_column.tsv", sep="\t", index=False)

# ---- heatmap z-scores and participant order
ab = pd.read_csv(P / "panel4c_adjusted_abundance.tsv", sep="\t").set_index("gene").loc[rows.gene]
meta = pd.read_csv(P / "panel4c_metadata.tsv", sep="\t").set_index("sample_id").loc[ab.columns]
print("participants", len(meta), meta.group.value_counts().to_dict())
print("batch x group", pd.crosstab(meta.acquisition_batch, meta.group).to_dict())
z = ab.sub(ab.mean(axis=1), axis=0).div(ab.std(axis=1, ddof=1), axis=0)
order = meta.sort_values(["Fibrosis", "NAS"], kind="mergesort").index
z[order].to_csv(OUT / "5c_recomputed_zscores_fib_nas_order.tsv", sep="\t")
meta.loc[order].to_csv(OUT / "5c_recomputed_participant_order.tsv", sep="\t")
print("NA cells per protein", ab.isna().sum()[ab.isna().sum() > 0].to_dict())
print("controls in order positions", [i for i, s in enumerate(order) if meta.loc[s, "group"] == "Control"])
print("z range", np.nanmin(z.values), np.nanmax(z.values), "cells beyond |2|", int((z.abs() > 2).sum().sum()))

# ---- partial Spearman: rank-residual on covariates batch, age, bmi, sex
def partial_spearman(x, y, cov):
    ok = x.notna() & y.notna() & cov.notna().all(axis=1)
    x, y, cv = x[ok], y[ok], cov[ok]
    X = np.column_stack([np.ones(ok.sum()), cv.values.astype(float)])
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    ex = rx - X @ np.linalg.lstsq(X, rx, rcond=None)[0]
    ey = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    rho = np.corrcoef(ex, ey)[0, 1]
    n = int(ok.sum()); k = X.shape[1] - 1
    df = n - k - 2
    t = rho * np.sqrt(df / (1 - rho ** 2))
    p = 2 * stats.t.sf(abs(t), df)
    return rho, n, p, df

cov = pd.DataFrame({"batch": (meta.acquisition_batch == meta.acquisition_batch.iloc[0]).astype(float),
                    "age": meta.age.astype(float), "bmi": meta.bmi.astype(float),
                    "sex": (meta.sex == "Female").astype(float)}, index=meta.index)
feat = ["Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS"]
out = []
for g in rows.gene:
    for f in feat:
        rho, n, p, df = partial_spearman(ab.loc[g], meta[f].astype(float), cov)
        out.append((g, f, rho, n, p, df))
ps = pd.DataFrame(out, columns=["gene", "feature", "rho", "n", "p", "df"])
ps["padj_bh125"] = stats.false_discovery_control(ps.p, method="bh")
ps.to_csv(OUT / "5c_recomputed_partial_spearman.tsv", sep="\t", index=False)
print(ps.groupby("gene").n.first().value_counts().to_dict())
# compare with analysis output
ana = pd.read_csv(P / "panel4c_histology_partial.tsv", sep="\t")
cmp = ps.merge(ana[["gene", "feature", "rho", "n", "pvalue", "padj"]], on=["gene", "feature"], suffixes=("", "_ana"))
print("max |rho diff| vs analysis output", (cmp.rho - cmp.rho_ana).abs().max(), "n mismatches", (cmp.n != cmp.n_ana).sum())
print("max |p diff|", (cmp.p - cmp.pvalue).abs().max())

# within-MASLD only (alternative)
mm = meta.group == "MASLD"
out2 = []
for g in rows.gene:
    for f in feat:
        rho, n, p, df = partial_spearman(ab.loc[g][mm], meta[f][mm].astype(float), cov[mm])
        out2.append((g, f, rho, n, p))
ps2 = pd.DataFrame(out2, columns=["gene", "feature", "rho_masld", "n_masld", "p_masld"])
ps2["padj_masld"] = stats.false_discovery_control(ps2.p_masld, method="bh")
ps2.to_csv(OUT / "5c_recomputed_partial_spearman_masld_only.tsv", sep="\t", index=False)
both = ps.merge(ps2, on=["gene", "feature"])
both["sign_flip"] = np.sign(both.rho) != np.sign(both.rho_masld)
print("within-MASLD: sign flips", both.sign_flip.sum(), "of", len(both), "; |rho| shrink median", (both.rho.abs() - both.rho_masld.abs()).median())
print("q<0.05 all-participant", (both.padj_bh125 < 0.05).sum(), "within-MASLD", (both.padj_masld < 0.05).sum())

# ---- process NES / q
pe = pd.read_csv(P / "panel4c_process_enrichment.tsv", sep="\t")
print(pe[["program", "representable", "NES", "padj", "size", "family"]].to_string())
pe.to_csv(OUT / "5c_process_enrichment_copy.tsv", sep="\t", index=False)
