#!/usr/bin/env python
"""
Build the GLP-1RA disease-reference table (Analysis B/F scaffold).

Per human gene: disease direction/effect + convergence standing + DEG status.
This is the fixed target scaffold onto which the human and mouse semaglutide
reversal signatures are joined to call reversed-vs-refractory.

Disease direction = sign of the canonical pooled bulk logFC (limma-voom-qw C2).
Tier-1 DEG gate    = treat_fdr < 0.05 (effect-size-aware interval-null FDR).
Convergence set    = 46d convergence_evidence.csv tiers / n_modalities_active.

Output: RNA-seq/results/glp1ra/disease_reference.csv
"""
import os
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CONV = os.path.join(ROOT, "RNA-seq/results/multi_evidence/convergence_evidence.csv")
DEG = os.path.join(ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/"
                         "results/integration/canonical_deg_results.csv")
OUT = os.path.join(ROOT, "RNA-seq/results/glp1ra/disease_reference.csv")

deg = pd.read_csv(DEG)
# symbol column holds the gene symbol; fall back to gene
sym = "symbol" if "symbol" in deg.columns else "gene"
deg = deg.rename(columns={sym: "human_symbol"})
deg_keep = deg[["human_symbol", "logFC", "treat_lfc", "treat_p", "treat_fdr",
                "padj", "AveExpr"]].copy()
deg_keep = deg_keep.dropna(subset=["human_symbol"]).drop_duplicates("human_symbol")
deg_keep["disease_dir"] = np.sign(deg_keep["logFC"])
deg_keep["is_deg_tier1"] = deg_keep["treat_fdr"] < 0.05

conv = pd.read_csv(CONV, low_memory=False)
conv_keep = conv[["human_symbol", "convergence_score", "convergence_rank",
                  "tier", "n_modalities_active", "concordance_state",
                  "sign_S1", "excluded_from_ranking", "druggability_tier",
                  "coloc_genetic_pp4"]].copy()

ref = deg_keep.merge(conv_keep, on="human_symbol", how="outer")

# convenience flags
ref["conv_tier1"] = ref["tier"].astype(str).str.startswith("1")
ref["conv_tier12"] = ref["tier"].astype(str).str.match(r"^[12]")
ref["conv_ge3_modalities"] = ref["n_modalities_active"] >= 3
ref["conv_included"] = ref["excluded_from_ranking"] != True  # noqa: E712

ref.to_csv(OUT, index=False)

print("[info] wrote", OUT, "rows:", len(ref))
print("\n===== target-set sizes =====")
print(f"  genes with a disease logFC (tested):     {ref['logFC'].notna().sum()}")
print(f"  Tier-1 DEGs (treat_fdr<0.05):            {int(ref['is_deg_tier1'].sum())}")
print(f"    ... up / down:                         "
      f"{int(((ref['is_deg_tier1']) & (ref['disease_dir']>0)).sum())} / "
      f"{int(((ref['is_deg_tier1']) & (ref['disease_dir']<0)).sum())}")
print(f"  convergence Tier-1:                      {int(ref['conv_tier1'].sum())}")
print(f"  convergence Tier-1/2:                    {int(ref['conv_tier12'].sum())}")
print(f"  convergence >=3 active modalities:       {int(ref['conv_ge3_modalities'].sum())}")
print("\n===== convergence tier distribution =====")
print(ref["tier"].value_counts(dropna=False).to_string())
print("\n===== incretin-axis genes in reference =====")
for g in ["GLP1R", "GIPR", "GCGR", "GLP2R", "GCG", "DPP4"]:
    r = ref[ref["human_symbol"] == g]
    if len(r):
        r = r.iloc[0]
        print(f"  {g:6s} logFC={r['logFC'] if pd.notna(r['logFC']) else float('nan'):+.3f}  "
              f"treat_fdr={r['treat_fdr'] if pd.notna(r['treat_fdr']) else float('nan'):.3g}  "
              f"tier={r['tier']}  n_mod={r['n_modalities_active']}")
    else:
        print(f"  {g:6s} (not in reference)")
