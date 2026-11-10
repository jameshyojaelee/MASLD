#!/usr/bin/env python3
"""A2 step 1: gated SuSiE credible sets, their window variants, and COLOC r2 variant posteriors.

Writes (work/):
  cs_members.parquet   gated, deduplicated CS members with normalized PIP (hg19)
  cs_table.tsv         one row per CS (study, locus, cs, lead, class, gate details)
  windows.parquet      every variant SuSiE tested in each study-locus that holds a kept CS
  coloc_vars.parquet   SNP.PP.H4 rows of the 767 COLOC r2 primary signal pairs
  lift_in.bed          unique hg19 positions of window variants for liftover
No statistics are computed here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
A2 = ROOT / "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage"
W = A2 / "work"
FM = ROOT / "GWAS/finemapping/results"
CAND = ROOT / "Analysis/Multimodal_Program_Projection/candidates"
PRIMARY = CAND / "atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
AUDIT = CAND / "atac-context-v3-candidate-2026-08-11-r1/genetics/context/variant_liftover_audit.tsv"
TIER = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
MHC = (6, 28_477_797, 33_448_354)

CLASS_BY_TRAIT = {"NAFLD": "direct", "NASH": "direct", "PDFF": "direct",
                  "ALT": "liver_enzyme", "AST": "liver_enzyme", "GGT": "liver_enzyme",
                  "CHRONLIVER": "sequela", "CIRRHOSIS": "sequela",
                  "ALBUMIN": "calib_albumin", "PLATELET": "calib_platelet"}

# ---- credible sets -------------------------------------------------------------------
cs = pd.read_csv(FM / "credible_sets.csv", dtype={"locus": str}, low_memory=False)
n_rows = len(cs)
cs = cs[cs.susie_cs > 0].copy()
cs["trait_u"] = cs.trait.str.upper()
cs["trait_class"] = cs.trait_u.map(CLASS_BY_TRAIT)
assert cs.trait_class.notna().all(), cs.loc[cs.trait_class.isna(), "trait"].unique()
key = ["study", "locus", "susie_cs"]
g = cs.groupby(key)
tab = g.agg(n_members=("susie_pip", "size"), raw_pip_sum=("susie_pip", "sum"),
            converged=("susie_converged", "first"), reliable=("susie_reliable", "first"),
            trait=("trait_u", "first"), trait_class=("trait_class", "first"),
            ancestry=("ancestry", "first")).reset_index()
lead = cs.loc[g.susie_pip.idxmax(), key + ["variant_id", "chromosome", "position", "susie_pip"]]
lead = lead.rename(columns={"variant_id": "lead_variant", "chromosome": "lead_chr",
                            "position": "lead_pos", "susie_pip": "lead_pip"})
tab = tab.merge(lead, on=key, validate="1:1")
tab["gate_pip_sum"] = tab.raw_pip_sum.between(0.90, 1.10)
tab["gate"] = tab.gate_pip_sum & (tab.converged == True) & (tab.reliable == True)
tab["mhc"] = (tab.lead_chr == MHC[0]) & tab.lead_pos.between(MHC[1], MHC[2])
# collapse CSs sharing a lead variant within study: keep highest lead PIP, then first locus id
tab = tab.sort_values(["study", "lead_variant", "lead_pip", "locus", "susie_cs"],
                      ascending=[True, True, False, True, True])
tab["dup_lead"] = False
gidx = tab.index[tab.gate]
tab.loc[gidx, "dup_lead"] = tab.loc[gidx].duplicated(["study", "lead_variant"]).to_numpy()
tab["kept"] = tab.gate & ~tab.dup_lead & ~tab.mhc
tab["cs_uid"] = tab.study + "|" + tab.locus + "|" + tab.susie_cs.astype(str)
tab.to_csv(W / "cs_table.tsv", sep="\t", index=False)
print("CS rows", n_rows, "CS total", len(tab), "gated", int(tab.gate.sum()),
      "dup", int(tab.dup_lead.sum()), "mhc", int((tab.gate & tab.mhc).sum()),
      "kept", int(tab.kept.sum()), flush=True)

kept = tab[tab.kept]
mem = cs.merge(kept[key + ["cs_uid", "raw_pip_sum"]], on=key, validate="m:1")
mem["pip_norm"] = mem.susie_pip / mem.raw_pip_sum
mem = mem[["cs_uid", "study", "locus", "susie_cs", "trait_class", "ancestry", "variant_id",
           "chromosome", "position", "allele1", "allele2", "susie_pip", "pip_norm"]]
mem.to_parquet(W / "cs_members.parquet", index=False)
print("members", len(mem), flush=True)

# ---- window variants -----------------------------------------------------------------
want = set(zip(kept.study, kept.locus))
parts, n_seen = [], 0
cols = ["chromosome", "position", "allele1", "allele2", "locus", "PIP", "study"]
for ch in pd.read_csv(FM / "susie_all_results.csv", usecols=cols, chunksize=4_000_000,
                      dtype={"locus": str, "study": str, "allele1": str, "allele2": str,
                             "chromosome": str}):
    n_seen += len(ch)
    m = pd.Series(list(zip(ch.study, ch.locus))).isin(want).to_numpy()
    if m.any():
        sub = ch.loc[m].copy()
        sub["chromosome"] = pd.to_numeric(sub.chromosome, errors="coerce")
        sub = sub[sub.chromosome.between(1, 22)]
        sub["chromosome"] = sub.chromosome.astype(np.int8)
        sub["position"] = sub.position.astype(np.int64)
        for c in ["allele1", "allele2", "study", "locus"]:
            sub[c] = sub[c].astype("string[pyarrow]")
        sub["PIP"] = sub.PIP.astype(np.float32)
        parts.append(sub)
    print("read", n_seen, flush=True)
win = pd.concat(parts, ignore_index=True)
del parts
win = win.drop_duplicates(["study", "locus", "chromosome", "position", "allele1", "allele2"])
win.to_parquet(W / "windows.parquet", index=False)
got = set(zip(win.study.astype(str), win.locus.astype(str)))
print("window rows", len(win), "study-loci wanted", len(want), "found", len(want & got),
      flush=True)
assert want <= got, list(want - got)[:5]

u = win[["chromosome", "position"]].drop_duplicates()
mu = mem[["chromosome", "position"]].drop_duplicates()
u = pd.concat([u, mu]).drop_duplicates().sort_values(["chromosome", "position"])
bed = pd.DataFrame({"chrom": "chr" + u.chromosome.astype(str),
                    "start0": u.position - 1, "end": u.position})
bed.to_csv(W / "lift_in.bed", sep="\t", header=False, index=False)
print("unique positions for liftover", len(bed), flush=True)

# ---- COLOC r2 primary-pair variant posteriors ----------------------------------------
prim = pd.read_csv(PRIMARY, sep="\t")
keys = prim[["gwas_name", "gene", "ensembl", "trait_class", "ancestry",
             "signal_pair_index", "pp_h4"]].drop_duplicates()
assert len(keys) == 767
keys["signal_pair_index"] = keys.signal_pair_index.astype(str)
keyset = set(zip(keys.gwas_name, keys.ensembl, keys.signal_pair_index))
acols = ["gwas_name", "ensembl", "signal_pair_index", "hg19_variant_id", "allele1", "allele2",
         "snp_pp_h4", "mapping_status", "hg38_chrom", "hg38_position_1based"]
parts = []
for ch in pd.read_csv(AUDIT, sep="\t", usecols=acols, chunksize=2_000_000,
                      dtype={"gwas_name": str, "ensembl": str, "signal_pair_index": str,
                             "hg19_variant_id": str, "mapping_status": str, "hg38_chrom": str,
                             "allele1": str, "allele2": str}):
    k = list(zip(ch.gwas_name, ch.ensembl, ch.signal_pair_index))
    m = np.fromiter((x in keyset for x in k), bool, len(k))
    if m.any():
        parts.append(ch.loc[m])
cv = pd.concat(parts, ignore_index=True)
del parts
cv = cv.merge(keys, on=["gwas_name", "ensembl", "signal_pair_index"], validate="m:1")
hg19 = cv.hg19_variant_id.str.split(":", expand=True)
cv["chromosome"] = pd.to_numeric(hg19[0], errors="coerce")
cv["position"] = pd.to_numeric(hg19[1], errors="coerce")
cv.to_parquet(W / "coloc_vars.parquet", index=False)
print("coloc variant rows", len(cv), "pairs", cv.groupby(["gwas_name", "ensembl"]).ngroups,
      flush=True)
