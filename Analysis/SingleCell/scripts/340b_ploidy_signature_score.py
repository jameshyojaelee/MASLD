#!/usr/bin/env python
"""
340b_ploidy_signature_score.py
Score mechanism-anchored ploidy marker categories per cell on the GSE244832 donor-matched
snRNA hepatocyte subset (the same 18 donors used by Phase A scPloidy).

Outputs:
  Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_percell.csv.gz
  Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_donor_summary.csv

This is the Phase A.5 cross-validation. Markers are from verified PubMed-indexed papers
(Sladky 2020/2022, Gentric 2015) — see data/ploidy_signatures/README.md. No claim of a
"published signature" — these are mechanism-anchored candidates.
"""
import os, sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
os.chdir(ROOT)

ATLAS = ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
MARKERS = ROOT / "data/ploidy_signatures/mechanism_markers.tsv"
OUT_DIR = ROOT / "Analysis/SingleCell/results_gpu_v2/ploidy"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"[{pd.Timestamp.now()}] Loading hepatocyte atlas {ATLAS}", flush=True)
adata = sc.read_h5ad(str(ATLAS))
print(f"  shape: {adata.shape}", flush=True)
print(f"  datasets: {adata.obs['dataset'].value_counts().to_dict()}", flush=True)

# Subset to GSE244832 (the donor-matched scATAC cohort)
gse_mask = adata.obs["dataset"].astype(str).eq("GSE244832")
if gse_mask.sum() == 0:
    # fallback: dataset key may differ - look for it
    print(f"  GSE244832 not found directly; available datasets: {adata.obs['dataset'].unique()}",
          flush=True)
    sys.exit("ERROR: GSE244832 not in atlas — adjust dataset filter")
adata = adata[gse_mask].copy()
print(f"[{pd.Timestamp.now()}] Subset to GSE244832: {adata.shape}", flush=True)
print(f"  samples: {adata.obs['sample'].nunique()}", flush=True)

# Normalize if needed (per project convention; scanpy score_genes wants log-normalized)
# Detect raw counts: if integer-like and max > 100, normalize.
mat_max = adata.X.max() if hasattr(adata.X, "max") else 0
mat_max = float(mat_max if not hasattr(mat_max, "item") else mat_max.item())
if mat_max > 100:
    print(f"  Detected raw counts (max={mat_max:.1f}); normalizing + log1p", flush=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# Cell-cycle gate (Tirosh markers, canonical)
s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7',
           'DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN',
           'WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45',
           'CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B',
           'BRIP1','E2F8']
g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2',
             'CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2',
             'AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','HN1',
             'CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2',
             'KIF23','HMMR','AURKA','PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2',
             'G2E3','GAS2L3','CBX5','CENPA']

# Filter to genes present in the data
present = set(adata.var_names.astype(str))
s_genes_present = [g for g in s_genes if g in present]
g2m_genes_present = [g for g in g2m_genes if g in present]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
print(f"  Cell-cycle gating: phase counts = "
      f"{adata.obs['phase'].value_counts().to_dict()}", flush=True)

# Load mechanism markers
markers = pd.read_csv(MARKERS, sep="\t")
print(f"[{pd.Timestamp.now()}] Markers loaded: {len(markers)} rows, "
      f"categories = {markers['category'].value_counts().to_dict()}", flush=True)

# Score each category separately; missing genes skipped
score_cols = []
for cat, sub in markers.groupby("category"):
    gene_list = [g for g in sub["gene"].tolist() if g in present]
    if len(gene_list) < 2:
        print(f"  [skip] {cat}: only {len(gene_list)} of {len(sub)} genes present", flush=True)
        continue
    score_name = f"sig_{cat}"
    try:
        sc.tl.score_genes(adata, gene_list=gene_list, score_name=score_name, use_raw=False)
        score_cols.append(score_name)
        print(f"  [ok] {cat}: scored {len(gene_list)}/{len(sub)} genes -> {score_name}",
              flush=True)
    except Exception as e:
        print(f"  [fail] {cat}: {e}", flush=True)

# Build per-cell output
keep_obs = ["sample", "condition", "phase", "S_score", "G2M_score",
            "hepatocyte_subtype_label", "axis_periportal", "axis_pericentral"]
keep_obs = [c for c in keep_obs if c in adata.obs.columns]
# library size proxy
if "log1p_total_counts" in adata.obs.columns:
    keep_obs.append("log1p_total_counts")
elif "total_counts" in adata.obs.columns:
    adata.obs["log1p_total_counts"] = np.log1p(adata.obs["total_counts"].astype(float))
    keep_obs.append("log1p_total_counts")
else:
    # Compute on the fly if not present
    counts = (adata.X.sum(axis=1) if hasattr(adata.X, "sum") else None)
    if counts is not None:
        counts = np.asarray(counts).ravel()
        adata.obs["log1p_total_counts"] = np.log1p(counts)
        keep_obs.append("log1p_total_counts")

# Try to map sample -> donor_id using sample_pairing if available
# (GSE244832 samples = MM_xxx / SRR_xxx; ATAC donors = D01..D18)
pairing = ROOT / "Analysis/ATAC/Human_Multiome/scripts/sample_pairing.csv"
sample_to_donor = None
if pairing.exists():
    try:
        pp = pd.read_csv(pairing, comment="#")
        # pp has columns: donor_id, atac_srr, rna_srrs, condition (mostly PLACEHOLDER)
        # We need a robust mapping. For now keep sample id as-is and note.
        pass
    except Exception:
        pass

per_cell = adata.obs[keep_obs + score_cols].copy()
per_cell["cell_barcode"] = adata.obs_names
per_cell = per_cell[["cell_barcode"] + keep_obs + score_cols]
out_csv = OUT_DIR / "signature_scores_percell.csv.gz"
per_cell.to_csv(out_csv, index=False, compression="gzip")
print(f"[{pd.Timestamp.now()}] Wrote {out_csv} ({per_cell.shape[0]} cells)", flush=True)

# Per-sample (donor) summary: mean score per sample, plus cell-cycle mix
gb_cols = score_cols + ["S_score", "G2M_score", "log1p_total_counts"]
gb_cols = [c for c in gb_cols if c in per_cell.columns]
donor_sum = per_cell.groupby("sample").agg(
    {"condition": "first",
     "phase": lambda s: (s == "G1").mean() if "phase" in per_cell.columns else np.nan,
     **{c: "mean" for c in gb_cols}}
).rename(columns={"phase": "frac_G1"}).reset_index()
out_sum = OUT_DIR / "signature_scores_donor_summary.csv"
donor_sum.to_csv(out_sum, index=False)
print(f"[{pd.Timestamp.now()}] Wrote {out_sum} ({donor_sum.shape[0]} samples)", flush=True)
print(donor_sum.head(20).to_string(index=False))
