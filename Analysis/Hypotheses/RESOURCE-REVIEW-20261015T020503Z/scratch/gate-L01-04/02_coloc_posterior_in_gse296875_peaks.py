#!/usr/bin/env python3
"""Gate L01-04 step 2 (login node, chunked read).

Use the COLOC shared-variant posteriors (SNP.PP.H4, hg38) from the ATAC v3 replay, restricted
to COLOC r2 primary signal pairs with PP.H4 > 0.5 on the 436 SuSiE-positive genes, and
measure posterior mass inside GSE296875 peaks (authors' 306,706 peaks, hg38). Hepatocyte
openness = donor-pseudobulk fragment sum over hepatocyte nuclei / 46,286.
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-04"
V3 = ROOT / "Analysis/Multimodal_Program_Projection/candidates"
LIFT = V3 / "atac-context-v3-candidate-2026-08-11-r1/genetics/liftover/variant_liftover_raw.tsv.gz"
PRIM = V3 / "atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
LENS = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L01/gene_trait_class_locus.tsv"
RECT = ROOT / "Analysis/MASLD_Model_Benchmark/executions/gse296875-crossmodal-recurrence-20260830T193000Z"
N_HEP = 46286

lens = pd.read_csv(LENS, sep="\t")
prim = pd.read_csv(PRIM, sep="\t")
prim = prim[(prim["lineage"] == "hepatocyte") & (prim["pp_h4"] > 0.5) & prim["ensembl"].isin(lens["ensembl"])]
prim = prim[["gwas_name", "gene", "ensembl", "signal_pair_index", "pp_h4", "any_accessible_mass"]].drop_duplicates()
prim = prim.merge(lens, on="ensembl")
print("r2 primary pairs (hepatocyte rows, PP.H4>0.5, lens genes):", len(prim),
      "genes", prim["ensembl"].nunique(), "loci", prim["locus"].nunique())
x = prim[prim["any_accessible_mass"] >= 0.5]
print("candidate's ATAC v3 count check (any_accessible_mass>=0.5):", len(x), x["locus"].nunique(), x["ensembl"].nunique())

keys = set(zip(prim["gwas_name"], prim["ensembl"], prim["signal_pair_index"].astype(int)))
cols = ["gwas_name", "ensembl", "signal_pair_index", "hg38_chrom", "hg38_position_1based", "SNP.PP.H4",
        "snp", "unique_liftover"]
keep = []
for ch in pd.read_csv(LIFT, sep="\t", usecols=cols, chunksize=2_000_000):
    ch = ch[ch["SNP.PP.H4"] >= 1e-3]
    k = list(zip(ch["gwas_name"], ch["ensembl"], ch["signal_pair_index"].astype(int)))
    ch = ch[[kk in keys for kk in k]]
    keep.append(ch)
post = pd.concat(keep)
print("variant rows kept (SNP.PP.H4>=1e-3):", len(post), "pairs covered:",
      post.groupby(["gwas_name", "ensembl", "signal_pair_index"]).ngroups, "of", len(keys))

pk = pd.read_csv(RECT / "rectangle/atac_peaks.tsv", sep="\t")
don = RECT / "donor"
atac_peaks = [l.strip() for l in open(don / "atac_peaks.txt")]
hep = np.load(don / "donor_atac_hepatocyte.npy")
hep_mean = pd.Series(hep.sum(axis=0) / N_HEP, index=atac_peaks)
pk["hep_frag_per_nucleus"] = pk["peak_id"].map(hep_mean)

hits = []
for chrom, block in post.groupby("hg38_chrom"):
    p = pk[pk["chromosome"] == chrom].sort_values("start_1based")
    s, e = p["start_1based"].to_numpy(), p["end_1based"].to_numpy()
    pos = block["hg38_position_1based"].to_numpy()
    idx = np.searchsorted(s, pos, side="right") - 1
    ok = (idx >= 0) & (pos <= e[np.clip(idx, 0, None)])
    b = block[ok].copy()
    b["peak_id"] = p["peak_id"].to_numpy()[idx[ok]]
    b["hep_frag_per_nucleus"] = p["hep_frag_per_nucleus"].to_numpy()[idx[ok]]
    hits.append(b)
hits = pd.concat(hits)
hits.to_csv(OUT / "coloc_variants_in_gse296875_peaks.tsv", sep="\t", index=False)

g = ["gwas_name", "ensembl", "signal_pair_index"]
mass = hits.groupby(g)["SNP.PP.H4"].sum().rename("mass_in_peak").reset_index()
best = hits.sort_values("SNP.PP.H4", ascending=False).drop_duplicates(g)[g + ["snp", "SNP.PP.H4", "peak_id", "hep_frag_per_nucleus"]]
res = prim.merge(mass, left_on=["gwas_name", "ensembl", "signal_pair_index"], right_on=g, how="left")
res = res.merge(best, on=g, how="left")
res["mass_in_peak"] = res["mass_in_peak"].fillna(0)
res.to_csv(OUT / "coloc_pairs_gse296875_peak_mass.tsv", sep="\t", index=False)

for thr in [0.2, 0.5]:
    y = res[res["mass_in_peak"] >= thr]
    print(f"pairs with COLOC posterior mass in GSE296875 peaks >= {thr}: pairs {len(y)}, loci {y['locus'].nunique()}, genes {y['ensembl'].nunique()}")
for thr in [0.1, 0.2, 0.5]:
    y = res[res["SNP.PP.H4"].fillna(0) >= thr]
    print(f"pairs whose single best in-peak variant has SNP.PP.H4 >= {thr}: pairs {len(y)}, loci {y['locus'].nunique()}, genes {y['ensembl'].nunique()}")

# Locus-level: best peak per locus (by mass); report hepatocyte openness of that peak
lp = hits.merge(prim[["gwas_name", "ensembl", "signal_pair_index", "locus"]], left_on=g, right_on=g)
lp = lp.groupby(["locus", "peak_id"]).agg(max_pp=("SNP.PP.H4", "max"),
                                          hep=("hep_frag_per_nucleus", "first")).reset_index()
lp = lp.sort_values("max_pp", ascending=False).drop_duplicates("locus")
lp.to_csv(OUT / "locus_best_gse296875_peak.tsv", sep="\t", index=False)
for thr in [0.1, 0.2, 0.5]:
    y = lp[lp["max_pp"] >= thr]
    print(f"loci with an in-peak variant SNP.PP.H4 >= {thr}: {len(y)}; peak hepatocyte fragments per nucleus (min, 10%, 50%, 90%):",
          np.round([y['hep'].min(), *np.nanquantile(y['hep'], [0.1, 0.5, 0.9])], 4) if len(y) else "-")
print("all GSE296875 peaks hepatocyte fragments per nucleus (5,25,50,75,95%):",
      np.round(np.nanquantile(pk["hep_frag_per_nucleus"], [0.05, 0.25, 0.5, 0.75, 0.95]), 4))
srt = res[res["gene"] == "SORT1"]
print(srt[["gwas_name", "signal_pair_index", "pp_h4", "any_accessible_mass", "locus", "mass_in_peak", "snp", "SNP.PP.H4", "peak_id", "hep_frag_per_nucleus"]].to_string())
