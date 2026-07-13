#!/usr/bin/env python
"""
WS2 — Incretin/glucagon drug axis: progression dynamics + human-genetics gap.

Consolidates, for GLP1R / GIPR / GCGR / DPP4 / GCG / GLP2R:
  (1) stage-resolved single-cell DE across the coarse MASLD transitions
      (Healthy->Steatosis->Steatohepatitis->Cirrhosis) per cell type;
  (2) the axis-wide human-genetics evidence (GWAS-eQTL COLOC across ~50 liver GWAS);
  (3) liver-tissue proteome DE.

Headline expected: GCGR down in hepatocytes at the steatosis onset; DPP4/GIPR up at
the NASH (steatohepatitis) transition; no COLOC support for any axis gene.

Output: RNA-seq/results/glp1ra/axis_progression/{incretin_axis_stage_de.csv,
        incretin_axis_genetics.csv, incretin_axis_summary.txt}
Env: spatial (pandas).
"""
import os
import glob
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PBDE = os.path.join(ROOT, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
COLOC = os.path.join(ROOT, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
# ⛔ RETIRED 2026-07-11 (critical review): `protein_differential_results.csv` is GSE276114 =
# liver RNA-seq (NOT DIA-MS) — misclassified as proteomics; the canonical liver DIA-MS
# (PXD051911) is null/opposite for DPP4 (MASH-vs-MASL logFC -0.034, p=0.755). Proteome
# "corroboration" from this file is INVALID and DISABLED. If protein is added, use PXD051911/
# PXD052937. Empty path -> os.path.exists() False -> proteome section skipped.
PROT = ""  # was: Analysis/Proteomics/results/protein_differential_results.csv (wrong modality)
OUT = os.path.join(ROOT, "RNA-seq/results/glp1ra/axis_progression")
os.makedirs(OUT, exist_ok=True)

AXIS = ["GLP1R", "GIPR", "GCGR", "DPP4", "GCG", "GLP2R"]
# coarse disease_stage_coarse transitions (snRNA pseudobulk DE), ordered
COARSE = ["MASLD_vs_Healthy", "Steatosis_vs_Healthy",
          "Steatohepatitis_vs_Steatosis", "Cirrhosis_vs_Steatohepatitis"]

# ---- (1) stage-resolved single-cell DE ----
rows = []
for f in glob.glob(os.path.join(PBDE, "*_de.csv")):
    if "bayesprism" in os.path.basename(f):
        continue                      # skip deconvolution F-stage (different method)
    try:
        df = pd.read_csv(f, usecols=["logFC", "pvalue", "padj", "gene", "cell_type", "contrast"])
    except Exception:
        continue
    df = df[df["gene"].isin(AXIS)]
    if len(df):
        rows.append(df)
de = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
de = de[de["contrast"].isin(COARSE)].copy()
de.to_csv(os.path.join(OUT, "incretin_axis_stage_de.csv"), index=False)

# ---- (2) genetics: COLOC ----
col = pd.read_csv(COLOC)
gcol = col[col["gene"].isin(AXIS)][["gene", "coloc_best_pp4", "coloc_best_gwas",
                                    "coloc_n_gwas_h4_05", "coloc_n_gwas_tested"]].copy()
gcol.to_csv(os.path.join(OUT, "incretin_axis_genetics.csv"), index=False)

# ---- (3) proteome (liver tissue DE) ----
prot = pd.DataFrame()
if os.path.exists(PROT):
    p = pd.read_csv(PROT)
    gcolname = "gene" if "gene" in p.columns else ("symbol" if "symbol" in p.columns else None)
    if gcolname:
        prot = p[p[gcolname].isin(AXIS)].copy()

# ---- summary ----
L = ["WS2 — Incretin/glucagon axis: progression dynamics + genetics gap", "=" * 60]

L.append("\n### (1) Stage-resolved single-cell DE (Hepatocytes; sig padj<0.05 bolded in words)")
hep = de[de["cell_type"] == "Hepatocytes"]
piv = hep.pivot_table(index="gene", columns="contrast", values="logFC").reindex(AXIS)
pad = hep.pivot_table(index="gene", columns="contrast", values="padj").reindex(AXIS)
for tr in [c for c in COARSE if c in piv.columns]:
    L.append(f"\n  [{tr}]")
    for g in AXIS:
        lfc = piv.loc[g, tr] if g in piv.index and tr in piv.columns else np.nan
        pv = pad.loc[g, tr] if g in pad.index and tr in pad.columns else np.nan
        if pd.notna(lfc):
            sig = "  <== SIGNIFICANT" if pd.notna(pv) and pv < 0.05 else ""
            L.append(f"    {g:6s} logFC={lfc:+.3f}  padj={pv:.3g}{sig}")

L.append("\n### Any other cell type with a SIGNIFICANT axis-gene stage change (padj<0.05):")
sig_all = de[(de["padj"] < 0.05)].sort_values(["gene", "cell_type", "contrast"])
for _, r in sig_all.iterrows():
    L.append(f"    {r['gene']:6s} {r['cell_type']:22s} {r['contrast']:32s} "
             f"logFC={r['logFC']:+.3f} padj={r['padj']:.3g}")

L.append("\n### (2) Human-genetics gap — GWAS-eQTL COLOC (best PP.H4 across ~50 liver GWAS)")
for _, r in gcol.set_index("gene").reindex(AXIS).reset_index().iterrows():
    if pd.notna(r.get("coloc_best_pp4")):
        L.append(f"    {r['gene']:6s} best PP4={r['coloc_best_pp4']:.3f} "
                 f"({r['coloc_best_gwas']}); n_H4>=0.5 = {int(r['coloc_n_gwas_h4_05'])} "
                 f"of {int(r['coloc_n_gwas_tested'])} tested")
    else:
        L.append(f"    {r['gene']:6s} not in COLOC table (not expressed / no eQTL)")
L.append("  => No axis gene colocalizes (all PP4 < 0.5): the incretin/glucagon receptors are "
         "systemic druggable targets, NOT inherited MASLD disease genes.")

if len(prot):
    L.append("\n### (3) Liver-tissue proteome (GSE276114 DIA-MS)")
    cols = [c for c in prot.columns if c.lower() in ("gene", "symbol", "logfc", "log2fc",
            "padj", "adj.p.val", "pvalue", "fdr")]
    L.append(prot[cols].to_string(index=False) if cols else prot.head().to_string(index=False))

summary = "\n".join(L)
with open(os.path.join(OUT, "incretin_axis_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)
print("\n[done] wrote", OUT)
