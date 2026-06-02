#!/usr/bin/env python
"""
330_percell_biological_signatures.py

Analyses H2 + H3 + D1-deep (v1) — Per-cell biological signature scoring.

Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md

Scores per-cell across the integrated atlas (1.2M cells, 16 cell types,
5 conditions) using scanpy.tl.score_genes:

- H2 Senescence-SASP: core senescence, pro-fibrotic SASP, pro-inflammatory
  SASP, growth-arrest, SenMayo
- H3 Cell cycle: G1/S, G2/M
- D1 Macrophage state: Kupffer, LAM, M1, M2
- Metabolic: FAO, Lipogenic, ER_Stress, Ferroptotic

Outputs: per-cell scores table + per-(celltype,condition) summary.

Env: rapids_singlecell
"""
import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATLAS = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")
OUTDIR = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/percell_signatures")
os.makedirs(OUTDIR, exist_ok=True)

print(f"[1] Loading atlas from {ATLAS}...")
t0 = time.time()
adata = sc.read_h5ad(ATLAS)
print(f"  Shape: {adata.shape} (loaded in {time.time()-t0:.0f}s)")

# Ensure data is log-normalized
if adata.X.min() >= 0 and (adata.X.max() > 100 or adata.X.max() > 25):
    # Could be raw counts; check
    if adata.X.max() > 100:
        print("  Normalizing + log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

signatures = {
    "H2_senescence_core": ["CDKN2A","CDKN1A","TP53","GDF15","SERPINE1","IGFBP3",
                            "SERPINB2","MMP3","GLB1","TNFRSF10D"],
    "H2_SASP_profibrotic": ["TIMP1","SERPINE1","TGFB1","TGFB2","CCN2","IL11",
                             "IGFBP3","IGFBP6","SPP1","CTGF","PAI1"],
    "H2_SASP_proinflammatory": ["IL6","IL8","CXCL1","CXCL2","CXCL3","CXCL8",
                                 "IL1A","IL1B","TNF","CCL2","CCL20"],
    "H2_SASP_growth_arrest": ["SPP1","IGFBP3","IGFBP6","IGFBP7","CDKN1A","CDKN2A",
                               "TNFRSF10D","TNFRSF10B"],
    "H2_SenMayo": ["ACVR1B","ANG","AREG","AXL","BMP2","BMP6","C3","CCL2","CCL20",
                    "CCL3","CCL4","CCL5","CCL7","CD55","CD9","CSF1","CSF2","CXCL1",
                    "CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL8","EDN1","EGF",
                    "EGFR","EREG","ESM1","FAS","FGF1","FGF2","FGF7","GDF15","HGF",
                    "HMGB1","ICAM1","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4",
                    "IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18","IL1A",
                    "IL1B","IL2","IL6","IL6ST","IL7","INHA","ITGA2","JUN","KITLG",
                    "MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2","MMP3",
                    "MMP9","PECAM1","PGF","PLAT","PLAU","PLAUR","PTGES","SELPLG",
                    "SERPINB3","SERPINB4","SERPINE1","SERPINE2","SPP1","TIMP2",
                    "TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","VEGFA",
                    "VEGFC","VGF","WNT2"],

    "H3_G1S": ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
                "CDCA7","DTL","PRIM1","UHRF1","CENPU","HELLS","RFC2","RPA2","NASP",
                "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
                "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
                "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"],
    "H3_G2M":  ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
                "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","PIMREG","SMC4","CCNB2",
                "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
                "KIF20B","HJURP","CDCA3","CDC20","TTK","CDC25C","KIF2C","RANGAP1",
                "NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1",
                "ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"],

    "D1_Kupffer": ["CLEC4F","TIMD4","MARCO","VSIG4","CD5L","CD163L1","FCN1","CETP","LYVE1"],
    "D1_LAM":     ["TREM2","CD9","GPNMB","SPP1","LIPA","LGALS3","FABP5","CTSB","CTSD",
                    "CTSL","APOE","APOC1","PLIN2","CD36"],
    "D1_M1":      ["CCL2","TNF","IL6","IL1B","CXCL9","CXCL10","NOS2","IFIT1","MX1","IRF1"],
    "D1_M2":      ["CD163","MRC1","STAB1","CD206","IL10","TGFB1","IGF1"],

    "Metab_FAO":         ["CPT1A","CPT1B","CPT2","ACOX1","ACOX2","HMGCS2","ACADVL",
                            "ACAT1","PPARA","PPARGC1A","HADHA","HADHB"],
    "Metab_Lipogenic":   ["FASN","SCD","ACACA","ACACB","SREBF1","INSIG1","MLXIPL",
                            "DGAT1","DGAT2"],
    "Metab_ER_Stress":   ["HSPA5","ATF4","DDIT3","XBP1","EIF2AK3","ERN1","ATF6"],
    "Metab_Ferroptotic": ["GPX4","SLC7A11","ALOX15","PTGS2","HMOX1","FTH1","FTL","NCOA4"],
}

# Filter signatures to genes present
var_names = set(adata.var_names)
for k in list(signatures.keys()):
    signatures[k] = [g for g in signatures[k] if g in var_names]

print("[2] Scoring signatures per cell...")
for name, genes in signatures.items():
    if len(genes) < 3:
        print(f"  {name}: skip (only {len(genes)} genes)")
        continue
    t0 = time.time()
    sc.tl.score_genes(adata, gene_list=genes, score_name=name,
                      use_raw=False, random_state=42)
    print(f"  {name}: n={len(genes)} ({time.time()-t0:.0f}s)")

score_cols = [c for c in adata.obs.columns if any(c.startswith(k)
              for k in ["H2_","H3_","D1_","Metab_"])]
print(f"  Computed {len(score_cols)} score columns.")

# Per-cell export (slim)
keep_cols = ["sample","dataset","condition","cell_type","cell_type_raw"] + score_cols
cell_df = adata.obs[keep_cols].copy()
cell_df.to_csv(os.path.join(OUTDIR, "percell_signature_scores.csv.gz"),
               compression="gzip", index=True)

# Aggregate: mean per (cell_type, condition)
agg = (cell_df.groupby(["cell_type","condition"], observed=True)[score_cols]
       .mean().reset_index())
agg.to_csv(os.path.join(OUTDIR, "mean_score_by_celltype_condition.csv"), index=False)

# Aggregate: mean per cell_type
agg_ct = cell_df.groupby("cell_type", observed=True)[score_cols].mean().reset_index()
agg_ct.to_csv(os.path.join(OUTDIR, "mean_score_by_celltype.csv"), index=False)

# Disease vs healthy delta per cell type (per signature)
deltas = []
for ct in cell_df["cell_type"].unique():
    sub = cell_df[cell_df.cell_type == ct]
    if len(sub["condition"].unique()) < 2:
        continue
    healthy = sub[sub.condition == "Healthy"]
    masld = sub[sub.condition.isin(["MASLD","NAFLD","NASH","Cirrhotic"])]
    if len(healthy) < 50 or len(masld) < 50:
        continue
    for s in score_cols:
        deltas.append({
            "cell_type": ct, "signature": s,
            "mean_healthy": healthy[s].mean(),
            "mean_disease": masld[s].mean(),
            "delta": masld[s].mean() - healthy[s].mean(),
            "n_healthy": len(healthy), "n_disease": len(masld),
        })
deltas_df = pd.DataFrame(deltas)
deltas_df.to_csv(os.path.join(OUTDIR, "disease_minus_healthy_per_celltype.csv"), index=False)

# Summary
lines = [
    f"Atlas: {adata.shape[0]} cells × {adata.shape[1]} genes",
    f"Signatures scored: {len(score_cols)}",
    "",
    "=== Top 20 disease-enriched (celltype, signature) by delta ===",
    deltas_df.sort_values("delta", ascending=False).head(20).to_string(index=False),
    "",
    "=== Top 20 healthy-enriched by -delta ===",
    deltas_df.sort_values("delta", ascending=True).head(20).to_string(index=False),
]
with open(os.path.join(OUTDIR, "percell_signatures_summary.txt"), "w") as f:
    f.write("\n".join(lines))
print("\n".join(lines[:30]))
print("\nDone. Outputs in:", OUTDIR)
