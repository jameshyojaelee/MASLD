#!/usr/bin/env python
"""
340d_ploidy_signature_full_atlas.py
Score verified mechanism-anchored ploidy markers on the full-gene integrated atlas
(37,533 genes), not the HVG-filtered hep atlas (2,931 genes).

Inputs:
  Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad  (1.23M cells, all genes)
  data/ploidy_signatures/mechanism_markers.tsv
Outputs:
  Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_fullatlas_percell.csv.gz
  Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_fullatlas_sample_summary.csv
  Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_fullatlas_donor_summary.csv

Phase A.5 v2 — addresses the gene-filter issue from the hep-atlas version.
"""
import os, gc, sys
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc, anndata as ad

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
os.chdir(ROOT)
sc.settings.verbosity = 2

ATLAS = ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
MARKERS = ROOT / "data/ploidy_signatures/mechanism_markers.tsv"
PAIRING = ROOT / "data/GSE244832/metadata/donor_pairing.csv"
OUT_DIR = ROOT / "Analysis/SingleCell/results_gpu_v2/ploidy"

print(f"[{pd.Timestamp.now()}] Loading full atlas (backed='r' for speed)...", flush=True)
a = ad.read_h5ad(str(ATLAS), backed="r")
print(f"  shape: {a.shape}", flush=True)
print(f"  datasets: {a.obs['dataset'].value_counts().to_dict()}", flush=True)

# Subset to GSE244832 hepatocytes (donor-matched to scATAC) and other liver hepatocyte sets
# For Phase A.5, focus on GSE244832 to enable scPloidy correlation
gse_mask = a.obs["dataset"].astype(str).eq("GSE244832")
hep_mask = a.obs["cell_type"].astype(str).eq("Hepatocytes")
mask = gse_mask & hep_mask
print(f"  GSE244832 hepatocytes: {mask.sum()}", flush=True)
print(f"[{pd.Timestamp.now()}] Loading subset into memory...", flush=True)
adata = a[mask].to_memory()
del a; gc.collect()
print(f"  Subset shape: {adata.shape}", flush=True)
print(f"  samples: {adata.obs['sample'].nunique()}", flush=True)

# Check counts layer / X (should be log-normalized; raw is also present)
print(f"  X.max(): {float(adata.X.max())}  raw is None: {adata.raw is None}", flush=True)
if adata.raw is not None:
    print(f"  raw.X.max(): {float(adata.raw.X.max())}", flush=True)
# Use the .raw (full gene set, log-normalized in integration) for scoring
if adata.raw is not None:
    print(f"[{pd.Timestamp.now()}] Using adata.raw for scoring (full gene set)", flush=True)
    use_raw = True
else:
    use_raw = False

# Tirosh cell-cycle gating (canonical)
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

# Reference: cell-cycle scoring uses adata.X by default. Use raw if available.
gene_pool = adata.raw.var_names if use_raw else adata.var_names
present_pool = set(gene_pool.astype(str))
s_g = [g for g in s_genes if g in present_pool]
g2m_g = [g for g in g2m_genes if g in present_pool]
print(f"  cell-cycle S genes present: {len(s_g)}/{len(s_genes)};  G2M genes: {len(g2m_g)}/{len(g2m_genes)}",
      flush=True)

# Workaround: scanpy.tl.score_genes_cell_cycle works on adata.X.
# If using raw, we need to swap raw into X temporarily.
if use_raw:
    adata_score = ad.AnnData(adata.raw.X.copy(), obs=adata.obs.copy(),
                              var=adata.raw.var.copy())
else:
    adata_score = adata
# Normalize if X looks like raw counts
mxv = float(adata_score.X.max())
print(f"  scoring matrix max: {mxv:.2f}", flush=True)
if mxv > 50:
    print(f"[{pd.Timestamp.now()}] Normalizing + log1p", flush=True)
    sc.pp.normalize_total(adata_score, target_sum=1e4)
    sc.pp.log1p(adata_score)
    print(f"  post-norm max: {float(adata_score.X.max()):.2f}", flush=True)

sc.tl.score_genes_cell_cycle(adata_score, s_genes=s_g, g2m_genes=g2m_g)
print(f"  cell-cycle phase counts: {adata_score.obs['phase'].value_counts().to_dict()}",
      flush=True)

# Load mechanism markers
markers = pd.read_csv(MARKERS, sep="\t")
print(f"\n[{pd.Timestamp.now()}] Mechanism markers: {len(markers)} rows", flush=True)
score_cols = []
for cat, sub in markers.groupby("category"):
    glist = [g for g in sub["gene"].tolist() if g in present_pool]
    if len(glist) < 2:
        print(f"  [skip] {cat}: only {len(glist)}/{len(sub)} genes present", flush=True)
        continue
    sname = f"sig_{cat}"
    sc.tl.score_genes(adata_score, gene_list=glist, score_name=sname, use_raw=False)
    score_cols.append(sname)
    print(f"  [ok]   {cat}: {len(glist)}/{len(sub)} genes -> {sname}", flush=True)

# Carry scores back to adata.obs and write outputs
keep = ["sample", "condition", "phase", "S_score", "G2M_score"]
keep = [c for c in keep if c in adata_score.obs.columns]
# Library size
counts_raw_max = adata_score.obs.get("total_counts")
if counts_raw_max is None:
    csum = np.asarray(adata.raw.X.sum(axis=1) if use_raw else adata.X.sum(axis=1)).ravel()
    adata_score.obs["log1p_total_counts"] = np.log1p(csum)
keep.append("log1p_total_counts" if "log1p_total_counts" in adata_score.obs.columns
            else "total_counts")
out_cols = ["cell_barcode"] + keep + score_cols
df = adata_score.obs[keep + score_cols].copy()
df.insert(0, "cell_barcode", adata_score.obs_names)
out_pc = OUT_DIR / "signature_scores_fullatlas_percell.csv.gz"
df.to_csv(out_pc, index=False, compression="gzip")
print(f"\n[{pd.Timestamp.now()}] Wrote {out_pc}: {df.shape}", flush=True)

# Per-sample summary
agg = {c: "mean" for c in score_cols + ["S_score", "G2M_score"]
       if c in df.columns}
agg.update({"phase": lambda s: (s == "G1").mean() if "phase" in df.columns else np.nan})
if "log1p_total_counts" in df.columns:
    agg["log1p_total_counts"] = "mean"
per_sample = df.groupby("sample").agg({"condition": "first", **agg}).reset_index()
per_sample.rename(columns={"phase": "frac_G1"}, inplace=True)
out_ss = OUT_DIR / "signature_scores_fullatlas_sample_summary.csv"
per_sample.to_csv(out_ss, index=False)
print(f"  Wrote {out_ss}: {per_sample.shape}", flush=True)

# Per-donor summary using donor_pairing.csv (snRNA SRRs -> donor_id)
pairs = pd.read_csv(PAIRING, comment="#")
pairs_long = pairs[["donor_id", "rna_srrs", "condition"]].copy()
pairs_long["rna_srrs"] = pairs_long["rna_srrs"].str.split(";")
pairs_long = pairs_long.explode("rna_srrs").rename(
    columns={"rna_srrs": "sample"})
merged = per_sample.merge(pairs_long[["donor_id", "sample"]], on="sample", how="left")
print(f"  Merge: {merged['donor_id'].isna().sum()} of {len(merged)} samples unmapped", flush=True)
agg2 = {c: "mean" for c in score_cols + ["S_score", "G2M_score", "frac_G1",
                                          "log1p_total_counts"]
        if c in merged.columns}
per_donor = merged.dropna(subset=["donor_id"]).groupby("donor_id").agg(agg2).reset_index()
# Add donor condition from pairing
per_donor = per_donor.merge(pairs[["donor_id", "condition"]], on="donor_id", how="left")
out_pd = OUT_DIR / "signature_scores_fullatlas_donor_summary.csv"
per_donor.to_csv(out_pd, index=False)
print(f"  Wrote {out_pd}: {per_donor.shape}", flush=True)
print("\nPer-donor table:\n", per_donor.to_string(index=False))
