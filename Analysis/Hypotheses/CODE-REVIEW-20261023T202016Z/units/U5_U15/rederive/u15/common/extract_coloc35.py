"""Extract the 35 primary strata from the adopted COLOC r2 master (independent code).

Writes coloc35.csv.gz with one row per master row in the 35 'main' strata and a log of
row counts. No randomness.
"""
import hashlib
import sys
from pathlib import Path

import pandas as pd

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
MASTER = R / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
TIER = R / "GWAS/finemapping/config/gwas_trait_tier.tsv"
OUT = Path(__file__).resolve().parent

tier = pd.read_csv(TIER, sep="\t")
main = tier[tier["placement"] == "main"]
print("main strata:", len(main), "tier counts:", main["tier"].value_counts().to_dict())
assert len(main) == 35

cols = ["gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H4.susie", "method",
        "top_snp", "top_snp_PP", "ancestry", "ld_panel", "n_cs_pairs", "eligibility_source"]
m = pd.read_csv(MASTER, usecols=cols, dtype={"chr": str})
print("master rows:", len(m), "strata:", m["gwas_name"].nunique())
x = m[m["gwas_name"].isin(main["study_name"])].copy()
print("rows in 35 strata:", len(x), "strata present:", x["gwas_name"].nunique())
missing = set(main["study_name"]) - set(x["gwas_name"])
print("strata with no rows:", sorted(missing))
x = x.merge(main[["study_name", "trait", "tier"]], left_on="gwas_name", right_on="study_name",
            how="left", validate="many_to_one").drop(columns="study_name")
assert x["tier"].notna().all()
dup = x.duplicated(["gwas_name", "gene", "ensembl"]).sum()
print("duplicate (gwas, gene, ensembl) rows:", dup)
print("duplicate (gwas, gene) rows:", x.duplicated(["gwas_name", "gene"]).sum())
print("rows with empty gene symbol:", x["gene"].isna().sum())
print("method values:", x["method"].value_counts(dropna=False).to_dict())
x.to_csv(OUT / "coloc35.csv.gz", index=False)

h = hashlib.sha256()
with open(MASTER, "rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 24), b""):
        h.update(chunk)
print("master sha256:", h.hexdigest())
print("python", sys.version)
