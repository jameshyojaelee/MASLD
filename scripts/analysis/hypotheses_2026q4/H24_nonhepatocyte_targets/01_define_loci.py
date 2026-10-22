#!/usr/bin/env python3
"""H24 step 1: direct-disease loci and their bulk-eQTL colocalization status.

Inputs only; no H24 outcome is read here.

Loci are the operational 1-Mb clusters used for the corrected Figure 2
protein-altering analysis (scripts/figures/fig2_multiplicity_correction.py):
gated credible sets from the DNA-only fine-mapping substrate, clustered by
single linkage (gap > 1 Mb starts a new cluster) on sentinel positions pooled
across the direct and enzyme trait scopes. A cluster is a direct-disease locus
if it holds at least one gated direct (tier-1 NAFLD/NASH/PDFF) credible set.
The expected count is 24 (docs/RESULTS.md, corrected 2026-09-23).

Bulk-eQTL colocalization comes from the completed September COLOC r2 candidate
(GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/; not adopted, the
corrected 50-study rerun is active). A gene is colocalized in a cluster window
(direct member span +/- 500 kb, hg19) when a tier-1/2 study gives
eligible SuSiE PP.H4 > 0.5 (hit1 in window) or ABF PP.H4 > 0.5 (top_snp in
window), the same > 0.5 rule as the release tier-1/2 sets.

Usage: 01_define_loci.py <OUT_DIR>
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = ROOT / "GWAS/finemapping"
SUB = FM / "results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
ARCH = SUB / "credible_set_pip_architecture.tsv"
MEMB = SUB / "credible_set_members_annotated.tsv"
REL = FM / "results/susie_coloc_r2_20260924T053244Z"
PAIRS = REL / "susie_coloc_signal_pairs.tsv.gz"
ALLG = REL / "susie_coloc_all_gwas.csv"
TIER = FM / "config/gwas_trait_tier.tsv"

GAP = 1_000_000
WINDOW = 500_000
DIRECT = "tier1_direct_masld_pdff"
ENZYME = "tier2_liver_enzyme"
MHC = ("6", 25_000_000, 34_000_000)
EXPECTED_DIRECT_CLUSTERS = 24


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1"])


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    arch = pd.read_csv(ARCH, sep="\t", dtype={"locus": str})
    memb = pd.read_csv(MEMB, sep="\t", dtype={"locus": str, "chromosome": str},
                       usecols=["credible_set_uid", "study", "trait_scope", "ancestry",
                                "locus", "variant_id", "chromosome", "position",
                                "allele1", "allele2", "normalized_susie_pip",
                                "consequence_category", "pip_sum_gate_passed"])
    gate_arch = arch.loc[as_bool(arch["pip_sum_gate_passed"]), "credible_set_uid"]
    gate_memb = memb.groupby("credible_set_uid")["pip_sum_gate_passed"].apply(
        lambda s: as_bool(s).all())
    keep = set(gate_arch) & set(gate_memb[gate_memb].index)
    arch = arch[arch["credible_set_uid"].isin(keep)].copy()
    memb = memb[memb["credible_set_uid"].isin(keep)].copy()

    arch["chrom"] = arch["locus"].str.split(".").str[0]
    arch["pos"] = pd.to_numeric(arch["locus"].str.split(".").str[1])
    arch = arch.sort_values(["chrom", "pos"]).copy()
    keys = []
    for chrom, block in arch.groupby("chrom", sort=False):
        new = np.concatenate([[True], np.diff(block["pos"].to_numpy()) > GAP])
        keys.append(pd.Series([f"{chrom}:{i}" for i in np.cumsum(new)], index=block.index))
    arch["cluster_id"] = pd.concat(keys).reindex(arch.index)

    direct_cl = sorted(arch.loc[arch["trait_scope"] == DIRECT, "cluster_id"].unique())
    if len(direct_cl) != EXPECTED_DIRECT_CLUSTERS:
        raise SystemExit(f"direct clusters {len(direct_cl)} != {EXPECTED_DIRECT_CLUSTERS}")

    dmemb = memb[memb["trait_scope"] == DIRECT].merge(
        arch[["credible_set_uid", "cluster_id"]], on="credible_set_uid")
    dmemb = dmemb[dmemb["cluster_id"].isin(direct_cl)]

    tier = pd.read_csv(TIER, sep="\t")
    t12 = tier[tier["tier"].isin([1, 2])].set_index("study_name")["tier"].to_dict()

    pairs = pd.read_csv(PAIRS, sep="\t")
    pairs = pairs[as_bool(pairs["eligible"]) & (pairs["PP.H4"] > 0.5) & pairs["gwas"].isin(t12)]
    pairs = pairs.assign(src="susie", pos=pd.to_numeric(pairs["hit1"].str.split(":").str[1]),
                         chrom=pairs["hit1"].str.split(":").str[0], tier=pairs["gwas"].map(t12))
    allg = pd.read_csv(ALLG, usecols=["gwas_name", "gene", "ensembl", "PP.H4.abf", "top_snp"])
    allg = allg[(allg["PP.H4.abf"] > 0.5) & allg["gwas_name"].isin(t12) & allg["top_snp"].notna()]
    allg = allg.assign(src="abf", gwas=allg["gwas_name"],
                       chrom=allg["top_snp"].astype(str).str.split(":").str[0],
                       pos=pd.to_numeric(allg["top_snp"].astype(str).str.split(":").str[1],
                                         errors="coerce"),
                       tier=allg["gwas_name"].map(t12))
    coloc = pd.concat([pairs[["gene", "ensembl", "gwas", "tier", "src", "chrom", "pos"]],
                       allg[["gene", "ensembl", "gwas", "tier", "src", "chrom", "pos"]]])
    coloc["gene"] = coloc["gene"].fillna(coloc["ensembl"])

    rows = []
    for cl in direct_cl:
        m = dmemb[dmemb["cluster_id"] == cl]
        a = arch[(arch["cluster_id"] == cl) & (arch["trait_scope"] == DIRECT)]
        chrom = cl.split(":")[0]
        lo, hi = int(m["position"].min()), int(m["position"].max())
        w = coloc[(coloc["chrom"] == chrom) & (coloc["pos"] >= lo - WINDOW) & (coloc["pos"] <= hi + WINDOW)]
        d_genes = sorted(set(w.loc[w["tier"] == 1, "gene"]))
        e_genes = sorted(set(w.loc[w["tier"] == 2, "gene"]))
        top = m.sort_values("normalized_susie_pip", ascending=False).iloc[0]
        coding = (a["protein_altering_pip_mass"] + a["canonical_splice_pip_mass"]).mean()
        mhc = chrom == MHC[0] and hi >= MHC[1] and lo <= MHC[2]
        rows.append({
            "cluster_id": cl, "chrom": chrom, "span_start_hg19": lo, "span_end_hg19": hi,
            "n_direct_cs": a["credible_set_uid"].nunique(),
            "n_direct_members": m["variant_id"].nunique(),
            "direct_studies": ";".join(sorted(a["study"].unique())),
            "ancestries": ";".join(sorted(a["ancestry"].unique())),
            "top_member": top["variant_id"], "top_member_pip": round(float(top["normalized_susie_pip"]), 4),
            "mean_coding_pip_mass": round(float(coding), 4),
            "has_enzyme_cs": bool((arch.loc[arch["cluster_id"] == cl, "trait_scope"] == ENZYME).any()),
            "direct_coloc_genes": ";".join(d_genes), "n_direct_coloc_genes": len(d_genes),
            "enzyme_coloc_genes": ";".join(e_genes), "n_enzyme_coloc_genes": len(e_genes),
            "mhc": mhc,
            "h24_primary": (len(d_genes) == 0) and not mhc,
            "h24_strict_no_any_coloc": (len(d_genes) == 0) and (len(e_genes) == 0) and not mhc,
        })
    loci = pd.DataFrame(rows)
    loci.to_csv(out / "loci_direct_clusters.tsv", sep="\t", index=False)
    dmemb.to_csv(out / "direct_members_hg19.tsv", sep="\t", index=False)
    summary = {
        "release_state": "COLOC r2 candidate susie_coloc_r2_20260924T053244Z; not adopted; corrected 50-study rerun active",
        "direct_clusters": len(direct_cl),
        "h24_primary_loci": int(loci["h24_primary"].sum()),
        "h24_strict_loci": int(loci["h24_strict_no_any_coloc"].sum()),
        "mhc_clusters": int(loci["mhc"].sum()),
        "gated_direct_credible_sets": int(arch.loc[arch["trait_scope"] == DIRECT, "credible_set_uid"].nunique()),
        "inputs_sha256": {str(p): sha256(p) for p in [ARCH, MEMB, PAIRS, TIER]},
    }
    (out / "loci_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(loci.drop(columns=["enzyme_coloc_genes"]).to_string())


if __name__ == "__main__":
    main(Path(sys.argv[1]))
