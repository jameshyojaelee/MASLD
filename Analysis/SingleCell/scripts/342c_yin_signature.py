#!/usr/bin/env python
"""
342c_yin_signature.py - Build Yin 2024 polyploid signature from pre-computed DE tables.

Yin's DEGs_4n_changes.csv has Gene;Annotation pairs where:
  w_4n_exp     = WT only, 4n vs 2n (expression level)
  w_4n_nuclei  = WT only, 4n vs 2n (nuclei-frequency level)
  h_4n_*       = HNF4 KO only
  c_4n_*       = CEBPA KO only
  hc_4n_*      = HNF4 + CEBPA intersect
  ct_4n_*      = CEBPA + CTCF intersect
The PAPER does NOT specify direction (up/down) in the annotation. We pool all *_4n_*
entries as "ploidy-DE candidates" (direction-agnostic), with a weight equal to the
number of strain backgrounds the gene appears in (a robustness proxy).

We also output a separate WT-only restricted set (most conservative).

NOTE: Yin gene symbols are mouse. We uppercase to get human candidates.
"""
import pandas as pd
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT  = ROOT / "data/ploidy_signatures"

tbl = pd.read_csv(ROOT / "data/external/yin2024/tables/DEGs_4n_changes.csv", sep=";",
                  encoding="utf-8-sig")
print(f"rows: {len(tbl)}")
print(f"annotations:\n{tbl['Annotation'].value_counts()}")
print()

# Per-gene count of distinct annotations (robustness proxy)
gene_counts = tbl.groupby("Gene")["Annotation"].nunique().reset_index(name="n_annotations")
print(f"unique genes: {len(gene_counts)}")
print(f"genes with >=2 annotation tags: {(gene_counts.n_annotations >= 2).sum()}")
print(f"genes with >=3 annotation tags: {(gene_counts.n_annotations >= 3).sum()}")

# Map mouse symbol -> human via uppercase
gene_counts["human_symbol"] = gene_counts["Gene"].str.upper()

# Combined signature: ALL genes with at least one _4n_ annotation
# Sort by number of annotations (most conserved first)
gene_counts = gene_counts.sort_values(["n_annotations", "Gene"], ascending=[False, True])
gene_counts.to_csv(OUT / "yin2024_polyploid_candidates.tsv", sep="\t", index=False)
print(f"\nTOP 20 most-conserved Yin candidate genes:")
print(gene_counts.head(20).to_string(index=False))

# WT-only subset (cleanest)
wt_genes = tbl[tbl["Annotation"].str.startswith("w_")]["Gene"].unique().tolist()
wt_df = pd.DataFrame({"mouse_symbol": wt_genes, "human_symbol": [g.upper() for g in wt_genes]})
wt_df.to_csv(OUT / "yin2024_wt_only_polyploid.tsv", sep="\t", index=False)
print(f"\nWT-only Yin signature: {len(wt_df)} genes")
print(wt_df.to_string(index=False))

# Multi-strain conserved subset (>=2 annotations from different strain backgrounds)
# Each gene can have multiple tags. Count how many DIFFERENT strain backgrounds it appears in.
def strains_in(annots):
    # annots is iterable of strings like "w_4n_exp", "h_4n_nuclei"
    s = set()
    for a in annots:
        prefix = a.split("_")[0]  # w, h, c, t, hc, ct
        for letter in prefix:
            s.add(letter)
    return s

gene_strains = tbl.groupby("Gene")["Annotation"].apply(lambda x: strains_in(x)).reset_index()
gene_strains["n_strains"] = gene_strains["Annotation"].apply(len)
gene_strains["strains"] = gene_strains["Annotation"].apply(lambda x: ",".join(sorted(x)))
conserved = gene_strains[gene_strains["n_strains"] >= 2].copy()
conserved["human_symbol"] = conserved["Gene"].str.upper()
conserved = conserved.sort_values("n_strains", ascending=False)
conserved.to_csv(OUT / "yin2024_conserved_polyploid.tsv", sep="\t", index=False)
print(f"\nMulti-strain conserved Yin signature: {len(conserved)} genes (≥2 strain backgrounds)")
print(conserved[["Gene","strains","n_strains","human_symbol"]].head(20).to_string(index=False))

# Mechanism marker overlap
mm = pd.read_csv(ROOT / "data/ploidy_signatures/mechanism_markers.tsv", sep="\t")
mm_set = set(mm["gene"].str.upper())
all_yin_human = set(gene_counts["human_symbol"])
overlap = mm_set & all_yin_human
print(f"\nMechanism marker overlap (any Yin annotation): {sorted(overlap)} ({len(overlap)} of {len(mm_set)})")
