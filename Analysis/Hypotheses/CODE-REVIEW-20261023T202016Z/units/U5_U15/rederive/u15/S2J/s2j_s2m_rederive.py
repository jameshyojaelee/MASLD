"""Independent retabulation of Figure S2J and S2M from the MVP meSuSiE per-variant table.

S2J: gene maximum of the EUR-specific component PIP and of the 11 shared components, over
     ALT/AST/NAFLD loci, genes linked to variants via Broadaway liver cis-eQTL pairs
     (allele-orientation agnostic). Also the per-locus EUR-component maximum over all variants.
S2M: meSuSiE credible sets in the ALT/AST/NAFLD loci by ancestry sharing.
No randomness.
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = R / "GWAS/finemapping/results"
BW = R / "data/broadaway_eqtl"
OUTJ = Path(__file__).resolve().parent
OUTM = OUTJ.parent / "S2M"
TRAITS = ["ALT", "AST", "NAFLD"]
ANC = ["EUR", "EAS", "AFR", "AMR"]

mv = pd.read_csv(FM / "mesusie_mvp/mesusie_variant_summary_mvp.csv")
print("meSuSiE variant rows (all traits):", len(mv), "loci:", mv["locus_id"].nunique())
print("trait_pair values:", mv["trait_pair"].value_counts().to_dict())
mv = mv[mv["trait_pair"].isin(TRAITS)].copy()
print("ALT/AST/NAFLD variant rows:", len(mv), "loci:", mv["locus_id"].nunique())
print("loci by trait:", mv.groupby("trait_pair")["locus_id"].nunique().to_dict())

loc = pd.read_csv(FM / "mesusie_mvp/mesusie_locus_summary_mvp.csv")
locs = loc[loc["trait_pair"].isin(TRAITS)]
print("locus summary ALT/AST/NAFLD loci:", len(locs), "converged:", locs["converged"].value_counts().to_dict())
print("locus summary n_cs total:", locs["n_cs"].sum())

# ---------------- S2M ----------------
print("\n=== S2M")
print("cs_category values (variant rows in CS):",
      mv.loc[mv["in_cs"] == True, "cs_category"].value_counts(dropna=False).to_dict())
cs = mv[mv["in_cs"] == True].groupby(["locus_id", "cs_index"]).agg(
    cats=("cs_category", lambda s: ";".join(sorted(set(map(str, s))))),
    n_var=("variant_id", "size")).reset_index()
print("credible sets (locus, cs_index) with in_cs:", len(cs), "in loci:", cs["locus_id"].nunique())
print("CS with >1 category label:", (cs["cats"].str.contains(";")).sum())


def sharing(c):
    toks = c.split("_")
    if len(toks) >= 2:
        return "Model-shared (>=2 ancestries)"
    return f"{toks[0]}-restricted"


cs["sharing"] = cs["cats"].map(sharing)
tab = cs["sharing"].value_counts().rename_axis("sharing").reset_index(name="n")
tab["pct"] = (100 * tab["n"] / tab["n"].sum()).round(1)
print(tab.to_string(index=False))
# cross-check with locus summary n_cs_<combo> columns
ncols = [c for c in locs.columns if c.startswith("n_cs_")]
by_combo = locs[ncols].sum()
shared_ls = sum(v for k, v in by_combo.items() if len(k[5:].split("_")) >= 2)
print("locus-summary totals: shared", shared_ls,
      {k: int(v) for k, v in by_combo.items() if len(k[5:].split("_")) == 1})
# converged-only view
conv_ids = set(locs.loc[locs["converged"] == True, "locus_id"])
print("CS in converged loci only:", cs[cs["locus_id"].isin(conv_ids)].shape[0])
OUTM.mkdir(exist_ok=True)
tab.to_csv(OUTM / "S2M_rederived.csv", index=False)
cs.to_csv(OUTM / "S2M_credible_sets.csv", index=False)

# ---------------- S2J ----------------
print("\n=== S2J")
shared_cols = ["pip_" + "_".join(c) for k in (2, 3, 4) for c in itertools.combinations(ANC, k)]
assert all(c in mv.columns for c in shared_cols), shared_cols
print("shared component columns:", len(shared_cols))
mv["shared_max"] = mv[shared_cols].max(axis=1)

# per-locus EUR-component maximum over all variants of the locus
lmax = mv.groupby("locus_id").agg(eur_max=("pip_EUR", "max"), shared_max=("shared_max", "max")).reset_index()
lmax["eur_max_4dp"] = lmax["eur_max"].round(4)
in_band = lmax[(lmax["eur_max_4dp"] >= 0.4217) & (lmax["eur_max_4dp"] <= 0.4218)]
print("loci:", len(lmax), "; loci with EUR-component max 0.4217-0.4218:", len(in_band))
print(in_band.to_string(index=False))
lmax.to_csv(OUTJ / "S2J_locus_eur_max.csv", index=False)

# Broadaway join (allele-orientation agnostic); heavy (reads ~14 GB), skipped on the login node
import os
if os.environ.get("SKIP_BROADAWAY") == "1":
    raise SystemExit("SKIP_BROADAWAY=1: stopping before the Broadaway gene join")
mv["chr"] = mv["chr"].astype(int)
mv["pos"] = mv["pos"].astype(int)
mv["akey"] = [":".join(sorted([a, b])) for a, b in zip(mv["a1"].astype(str), mv["a2"].astype(str))]
need = mv[["chr", "pos"]].drop_duplicates()
parts = []
for ch in sorted(need["chr"].unique()):
    pos_set = set(need.loc[need["chr"] == ch, "pos"])
    kept = []
    for chunk in pd.read_csv(BW / f"chr{ch}_marginal_summary_results.tsv", sep="\t",
                             usecols=["CHR", "POS", "NEA", "EA", "GeneSymbol", "ENSG"],
                             chunksize=2_000_000, dtype={"CHR": int, "POS": int}):
        kept.append(chunk[chunk["POS"].isin(pos_set)])
    d = pd.concat(kept)
    print(f"  chr{ch}: {len(d)} gene-variant rows at needed positions")
    parts.append(d)
eq = pd.concat(parts, ignore_index=True)
eq["akey"] = [":".join(sorted([str(a), str(b)])) for a, b in zip(eq["NEA"], eq["EA"])]
eq = eq.rename(columns={"CHR": "chr", "POS": "pos"})[["chr", "pos", "akey", "GeneSymbol", "ENSG"]].drop_duplicates()
n_before = len(mv)
mj = mv.merge(eq, on=["chr", "pos", "akey"], how="inner")
print("variant rows before join:", n_before, "; gene-variant rows after join:", len(mj),
      "; distinct variants linked:", mj[["locus_id", "variant_id"]].drop_duplicates().shape[0])
mj = mj[mj["GeneSymbol"].notna() & (mj["GeneSymbol"] != "")]
g = mj.groupby("GeneSymbol").agg(eur_pip=("pip_EUR", "max"), joint_pip=("shared_max", "max"),
                                 n_loci=("locus_id", "nunique"),
                                 traits=("trait_pair", lambda s: ";".join(sorted(set(s))))).reset_index()
print("genes with a linked variant:", len(g))
CANON = ["GCKR", "PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", "TRIB1", "GPAM", "HNF1A", "HNF1B",
         "GATAD2A", "MLXIPL", "ABO", "HFE", "RECQL4", "GGT1", "SERPINA1", "MARC1", "MTARC1",
         "APOE", "PPP1R3B", "TOR1B"]
gc = g[g["GeneSymbol"].isin(CANON)].sort_values("joint_pip", ascending=False)
gc["eur_4dp"] = gc["eur_pip"].round(4)
gc["joint_4dp"] = gc["joint_pip"].round(4)
print(gc.to_string(index=False))
print("canonical genes absent:", sorted(set(CANON) - set(g["GeneSymbol"])))
g.to_csv(OUTJ / "S2J_gene_maxima.csv", index=False)
