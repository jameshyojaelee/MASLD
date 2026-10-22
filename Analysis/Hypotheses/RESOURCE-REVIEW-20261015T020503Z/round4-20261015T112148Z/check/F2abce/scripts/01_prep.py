"""Blind check F2abce, step 1: main-strata COLOC r2 rows and fine-mapping key set.
Reads only release files. Writes into the check dir."""
import os, sys
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CK = sys.argv[1]
OUT = os.path.join(CK, "out")
R2 = f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
tier = pd.read_csv(f"{ROOT}/GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
main = set(tier.loc[tier.placement == "main", "study_name"])
print("main strata", len(main))

usecols = ["gwas_name", "gene", "chr", "PP.H4.susie", "PP.H4.abf", "top_snp"]
d = pd.read_csv(R2, usecols=usecols, dtype={"gene": str, "top_snp": str, "gwas_name": str},
                keep_default_na=True, na_values=[""])
d["row"] = range(len(d))
print("r2 rows", len(d))
d = d[d.gwas_name.isin(main)]
print("main rows", len(d), "strata present", d.gwas_name.nunique())
# gene key: keep the literal string; pandas NaNs empty only (na_values default also
# includes 'NA','None' etc; record if any gene symbol could be lost that way)
d2 = pd.read_csv(R2, usecols=["gene"], dtype=str, keep_default_na=False)
lost = set(d2.gene.unique()) & {"NA", "None", "NaN", "nan", "null", "NULL", "N/A", "n/a"}
print("gene strings that default NA parsing would drop:", lost)
d = d[d.gene.notna() & (d.gene != "")]
d.to_csv(os.path.join(OUT, "coloc_main.tsv.gz"), sep="\t", index=False)
k = d[d.top_snp.notna() & (d.top_snp != "")].copy()
k[["c", "p"]] = k.top_snp.str.split(":", expand=True).iloc[:, :2]
keys = k[["gwas_name", "c", "p"]].drop_duplicates()
keys.to_csv(os.path.join(OUT, "keys.csv"), index=False, header=False)
print("unique study:chr:pos keys", len(keys))
