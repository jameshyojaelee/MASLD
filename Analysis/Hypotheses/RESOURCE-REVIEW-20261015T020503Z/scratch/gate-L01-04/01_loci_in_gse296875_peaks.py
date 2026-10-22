#!/usr/bin/env python3
"""Gate L01-04 step 1 (login node, small files only).

Count COLOC r2 SuSiE-positive loci whose GWAS credible set has posterior mass in a
GSE296875 peak (authors' 306,706-peak set, lifted to hg19 by H22), and how open those
peaks are in GSE296875 hepatocyte nuclei (donor pseudobulk sums / 46,286 nuclei).
Positions come from the `position` column, never from the decimal `locus` column.
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-04"
FM = ROOT / "GWAS/finemapping"
MEMB = FM / "results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv"
PAIRS = FM / "results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz"
TIER = FM / "config/gwas_trait_tier.tsv"
LENS = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L01/gene_trait_class_locus.tsv"
PEAKS19 = ROOT / "Analysis/Hypotheses/H22-20261009T015913Z/peaks/peaks_hg19.tsv"
DON = ROOT / "Analysis/MASLD_Model_Benchmark/executions/gse296875-crossmodal-recurrence-20260830T193000Z/donor"
N_HEP = 46286


def as_bool(s):
    return s.astype(str).str.lower().isin(["true", "1"])


lens = pd.read_csv(LENS, sep="\t")
print("lens genes", len(lens), "loci", lens["locus"].nunique())

tier = pd.read_csv(TIER, sep="\t")
t12 = set(tier.loc[tier["tier"].isin([1, 2]), "study_name"])
pairs = pd.read_csv(PAIRS, sep="\t")
pairs = pairs[as_bool(pairs["eligible"]) & (pairs["PP.H4"] > 0.5) & pairs["gwas"].isin(t12)]
pairs = pairs[pairs["ensembl"].isin(lens["ensembl"])].copy()
pairs["chrom"] = pairs["hit1"].str.split(":").str[0]
pairs["pos"] = pd.to_numeric(pairs["hit1"].str.split(":").str[1])
print("positive signal pairs on lens genes", len(pairs), "genes", pairs["ensembl"].nunique())

memb = pd.read_csv(MEMB, sep="\t", dtype={"chromosome": str},
                   usecols=["credible_set_uid", "study", "variant_id", "chromosome", "position",
                            "normalized_susie_pip", "pip_sum_gate_passed", "consequence_category"])
memb = memb[as_bool(memb["pip_sum_gate_passed"])]

# Map each COLOC (gwas, hit1) to credible sets of the same study containing hit1's position.
key = memb[["credible_set_uid", "study", "chromosome", "position"]].drop_duplicates()
m = pairs.merge(key, left_on=["gwas", "chrom", "pos"], right_on=["study", "chromosome", "position"], how="left")
print("pairs with a matching credible set:", m["credible_set_uid"].notna().sum(), "of", len(pairs))
m = m.dropna(subset=["credible_set_uid"])

# Peaks (hg19) and hepatocyte openness
pk = pd.read_csv(PEAKS19, sep="\t", dtype={"chr": str})
atac_peaks = [l.strip() for l in open(DON / "atac_peaks.txt")]
hep = np.load(DON / "donor_atac_hepatocyte.npy", mmap_mode="r")
print("donor_atac_hepatocyte shape", hep.shape)
hep_sum = np.asarray(hep).sum(axis=0) if hep.shape[1] == len(atac_peaks) else np.asarray(hep).sum(axis=1)
hep_mean = pd.Series(hep_sum / N_HEP, index=atac_peaks)
pk["hep_frag_per_nucleus"] = pk["peak_id"].map(hep_mean)
print("all peaks: hepatocyte fragments per nucleus quantiles (5,25,50,75,95%):",
      np.round(np.nanquantile(pk["hep_frag_per_nucleus"], [0.05, 0.25, 0.5, 0.75, 0.95]), 4))

cs = memb[memb["credible_set_uid"].isin(m["credible_set_uid"])].copy()
cs["position"] = cs["position"].astype(int)
hits = []
for chrom, block in cs.groupby("chromosome"):
    p = pk[pk["chr"] == chrom].sort_values("start")
    if p.empty:
        continue
    starts, ends = p["start"].to_numpy(), p["end"].to_numpy()
    idx = np.searchsorted(starts, block["position"].to_numpy(), side="right") - 1
    ok = (idx >= 0) & (block["position"].to_numpy() <= ends[np.clip(idx, 0, None)])
    b = block[ok].copy()
    b["peak_id"] = p["peak_id"].to_numpy()[idx[ok]]
    b["hep_frag_per_nucleus"] = p["hep_frag_per_nucleus"].to_numpy()[idx[ok]]
    hits.append(b)
hits = pd.concat(hits)
hits.to_csv(OUT / "credible_variants_in_gse296875_peaks.tsv", sep="\t", index=False)

mass = hits.groupby("credible_set_uid")["normalized_susie_pip"].sum().rename("mass_in_peak")
best = hits.sort_values("normalized_susie_pip", ascending=False).drop_duplicates("credible_set_uid")
best = best.set_index("credible_set_uid")[["variant_id", "normalized_susie_pip", "peak_id", "hep_frag_per_nucleus"]]
m = m.merge(lens, on="ensembl", how="left")
m = m.join(mass, on="credible_set_uid").join(best, on="credible_set_uid")
m["mass_in_peak"] = m["mass_in_peak"].fillna(0)
m.to_csv(OUT / "coloc_pairs_peak_mass.tsv", sep="\t", index=False)

for thr in [0.1, 0.2, 0.5]:
    x = m[m["normalized_susie_pip"].fillna(0) >= thr]
    print(f"top in-peak variant PIP >= {thr}: pairs {len(x)}, loci {x['locus'].nunique()}, genes {x['ensembl'].nunique()}")
for thr in [0.2, 0.5]:
    x = m[m["mass_in_peak"] >= thr]
    print(f"credible-set mass in GSE296875 peaks >= {thr}: loci {x['locus'].nunique()}, genes {x['ensembl'].nunique()}")
# Locus-level summary: best in-peak variant per locus
loc = m.sort_values("normalized_susie_pip", ascending=False).drop_duplicates("locus")
loc = loc[loc["normalized_susie_pip"].notna()]
loc[["locus", "gene", "gwas", "variant_id", "normalized_susie_pip", "mass_in_peak", "peak_id",
     "hep_frag_per_nucleus"]].to_csv(OUT / "locus_best_peak.tsv", sep="\t", index=False)
for thr in [0.1, 0.2, 0.5]:
    x = loc[loc["normalized_susie_pip"] >= thr]
    print(f"loci with best in-peak PIP >= {thr}: {len(x)}; their peak hepatocyte fragments/nucleus quantiles (10,50,90%):",
          np.round(np.nanquantile(x["hep_frag_per_nucleus"], [0.1, 0.5, 0.9]), 4) if len(x) else "-")
