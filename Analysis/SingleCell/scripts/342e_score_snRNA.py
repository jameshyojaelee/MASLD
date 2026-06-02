#!/usr/bin/env python
"""
342e_score_snRNA.py - Phase 2.1: score all polyploid signatures on the 657K-hepatocyte
integrated atlas across 5 disease cohorts.

Signatures scored:
  - Richter UP / DOWN (mouse snRNA-seq2 2n vs 4n)
  - Katsuda UP / DOWN (rat bulk microarray 2c vs 4c+8c)
  - Yin all-candidates / WT-only / conserved (mouse snRNA-seq2 4n vs 2n, direction-agnostic)
  - Consensus UP / DOWN (≥2 sources)
Composite indices:
  - polyploid_index_richter   = mean(Richter UP) - mean(Richter DOWN)
  - polyploid_index_katsuda   = mean(Katsuda UP) - mean(Katsuda DOWN)
  - polyploid_index_consensus = mean(consensus UP) - mean(consensus DOWN)

Outputs:
  signature_scores_snRNA_percell.csv.gz
  signature_scores_snRNA_sample_summary.csv
  signature_scores_snRNA_donor_summary.csv
"""
import os, gc
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
os.chdir(ROOT)
sc.settings.verbosity = 2

ATLAS = ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
SIG_DIR = ROOT / "data/ploidy_signatures"
OUT_DIR = ROOT / "Analysis/SingleCell/results_gpu_v2/ploidy"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load all signatures
def load_genes(path, col="human_symbol"):
    df = pd.read_csv(path, sep="\t")
    g = df[col].dropna().astype(str).str.upper().unique().tolist()
    return [x for x in g if x and x != "NAN"]

signatures = {
    "richter_up":   load_genes(SIG_DIR / "richter2021_polyploid_up.tsv"),
    "richter_down": load_genes(SIG_DIR / "richter2021_polyploid_down.tsv"),
    "katsuda_up":   load_genes(SIG_DIR / "katsuda2019_polyploid_up.tsv"),
    "katsuda_down": load_genes(SIG_DIR / "katsuda2019_polyploid_down.tsv"),
    "yin_all":      load_genes(SIG_DIR / "yin2024_polyploid_candidates.tsv"),
    "yin_conserved":load_genes(SIG_DIR / "yin2024_conserved_polyploid.tsv"),
    "consensus_up": load_genes(SIG_DIR / "consensus_polyploid_up.tsv"),
    "consensus_down": load_genes(SIG_DIR / "consensus_polyploid_down.tsv"),
}
print(f"Loaded signatures:")
for k, v in signatures.items():
    print(f"  {k}: {len(v)} genes")

# Load atlas (backed mode) then materialize hepatocyte subset
print(f"\n[{pd.Timestamp.now()}] Loading atlas (backed)", flush=True)
a = ad.read_h5ad(str(ATLAS), backed="r")
print(f"  shape: {a.shape}", flush=True)

# Subset to hepatocytes across all datasets
hep_mask = a.obs["cell_type"].astype(str).eq("Hepatocytes")
print(f"  hepatocytes: {hep_mask.sum()}", flush=True)
print(f"  datasets:")
for ds, n in a.obs.loc[hep_mask, "dataset"].value_counts().items():
    print(f"    {ds}: {n}")

print(f"\n[{pd.Timestamp.now()}] Loading hepatocyte subset to memory", flush=True)
adata = a[hep_mask].to_memory()
del a; gc.collect()
print(f"  shape: {adata.shape}", flush=True)

# Use .raw (full gene set) since adata.X may be HVG-filtered in some atlases
if adata.raw is not None:
    print(f"  Using adata.raw (n={adata.raw.X.shape[1]} genes)", flush=True)
    score_adata = ad.AnnData(adata.raw.X.copy(), obs=adata.obs.copy(), var=adata.raw.var.copy())
else:
    print(f"  Using adata.X (n={adata.X.shape[1]} genes)", flush=True)
    score_adata = adata
# Normalize if raw counts
mxv = float(score_adata.X.max())
print(f"  max value: {mxv:.2f}", flush=True)
if mxv > 50:
    print(f"  Normalizing + log1p", flush=True)
    sc.pp.normalize_total(score_adata, target_sum=1e4)
    sc.pp.log1p(score_adata)

# Cell-cycle gating
s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7',
           'DTL','PRIM1','UHRF1','HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN',
           'WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45',
           'CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1','CHAF1B',
           'BRIP1','E2F8']
g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2',
             'CKS1B','MKI67','TMPO','CENPF','TACC3','SMC4','CCNB2','CKAP2L','CKAP2',
             'AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3',
             'CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2',
             'KIF23','HMMR','AURKA','PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2',
             'G2E3','GAS2L3','CBX5','CENPA']
gene_pool = set(score_adata.var_names.astype(str))
s_g  = [g for g in s_genes if g in gene_pool]
g2m_g = [g for g in g2m_genes if g in gene_pool]
sc.tl.score_genes_cell_cycle(score_adata, s_genes=s_g, g2m_genes=g2m_g)
print(f"  cell-cycle phases: {score_adata.obs['phase'].value_counts().to_dict()}", flush=True)

# Score signatures
score_cols = []
for name, genes in signatures.items():
    present = [g for g in genes if g in gene_pool]
    if len(present) < 2:
        print(f"  [skip] {name}: only {len(present)} genes present", flush=True)
        continue
    sc.tl.score_genes(score_adata, gene_list=present, score_name=f"sig_{name}", use_raw=False)
    score_cols.append(f"sig_{name}")
    print(f"  [ok] {name}: scored {len(present)}/{len(genes)} genes", flush=True)

# Composite indices
def composite(adata, name, up_genes, dn_genes):
    pool = set(adata.var_names)
    up_present = [g for g in up_genes if g in pool]
    dn_present = [g for g in dn_genes if g in pool]
    if len(up_present) < 2 or len(dn_present) < 2:
        return None
    up_score_col = f"_tmp_up_{name}"; dn_score_col = f"_tmp_dn_{name}"
    sc.tl.score_genes(adata, up_present, score_name=up_score_col, use_raw=False)
    sc.tl.score_genes(adata, dn_present, score_name=dn_score_col, use_raw=False)
    adata.obs[f"polyploid_index_{name}"] = adata.obs[up_score_col] - adata.obs[dn_score_col]
    score_cols.append(f"polyploid_index_{name}")
    del adata.obs[up_score_col]; del adata.obs[dn_score_col]

for nm, up, dn in [("richter", signatures["richter_up"], signatures["richter_down"]),
                   ("katsuda", signatures["katsuda_up"], signatures["katsuda_down"]),
                   ("consensus", signatures["consensus_up"], signatures["consensus_down"])]:
    composite(score_adata, nm, up, dn)
print(f"\nComposite indices: polyploid_index_richter, _katsuda, _consensus", flush=True)

# Library size + age covariate
total_counts = np.asarray(adata.raw.X.sum(axis=1)).ravel() if adata.raw is not None else \
               np.asarray(adata.X.sum(axis=1)).ravel()
score_adata.obs["log1p_total_counts"] = np.log1p(total_counts)

# Per-cell output
keep_obs = ["sample","dataset","condition","cell_type",
            "phase","S_score","G2M_score","log1p_total_counts",
            "hepatocyte_subtype_label"]
keep_obs = [c for c in keep_obs if c in score_adata.obs.columns]
df = score_adata.obs[keep_obs + score_cols].copy()
df.insert(0, "cell_barcode", score_adata.obs_names)
out_pc = OUT_DIR / "signature_scores_snRNA_percell.csv.gz"
df.to_csv(out_pc, index=False, compression="gzip")
print(f"\n[{pd.Timestamp.now()}] Wrote {out_pc}: {df.shape}", flush=True)

# Per-sample summary
agg = {c: "mean" for c in score_cols + ["S_score","G2M_score","log1p_total_counts"] if c in df.columns}
agg["phase"] = lambda s: (s == "G1").mean()
per_sample = df.groupby("sample").agg({"dataset":"first", "condition":"first", **agg}).reset_index()
per_sample.rename(columns={"phase": "frac_G1"}, inplace=True)
out_ss = OUT_DIR / "signature_scores_snRNA_sample_summary.csv"
per_sample.to_csv(out_ss, index=False)
print(f"  Wrote {out_ss}: {per_sample.shape}", flush=True)

# Per-cell-type summary (sanity check: hep > other?)
# Need to score on non-hep too - skip for now, do as a separate sanity-check script
print(f"\nSample summary head:")
print(per_sample.head(20).to_string(index=False))
