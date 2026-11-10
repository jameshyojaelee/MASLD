"""Join RV5 (consistency class), RV8 (coding_neighbor) and RV9 (replication_state) into one candidate
Gene Catalog annotation table keyed by unversioned Ensembl ID. Candidate only; the r8 Catalog is unchanged.
Usage: python -I rvcat_build.py <new_out_dir>"""
import sys
from pathlib import Path
import pandas as pd

H = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses")
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=False)

rv5 = pd.read_csv(H / "RV5-20261102T144853Z/out-20261102T144930Z/call_heterogeneity.tsv", sep="\t")
rv5["ensembl"] = rv5.gene.str.replace(r"\.\d+$", "", regex=True)
rv5 = rv5.rename(columns={"symbol": "symbol_rv5", "class": "consistency_class", "tau": "re_tau", "I2": "re_I2",
                          "pi_lo": "re_pi_lo", "pi_hi": "re_pi_hi", "loo_fail_cohorts": "loo_loss_cohorts"})
rv5 = rv5[["ensembl", "symbol_rv5", "consistency_class", "re_tau", "re_I2", "re_pi_lo", "re_pi_hi", "loo_loss_cohorts"]]
assert len(rv5) == 1854 and rv5.ensembl.is_unique

rv8 = pd.read_csv(H / "RV8-20261102T150226Z/out-20261102T150342Z/gene_flags.tsv", sep="\t")
rv8 = rv8[rv8.analysis == "primary"]
LAB = {"anchor_is_other_coding": "anchor_is_other_coding", "not_separable": "not_separable", "moderate_LD": "moderate_LD",
       "moderate_ld": "moderate_LD", "own_gene_coding": "own_gene_coding", "ld_missing": "not_assessed", "none": "none"}
rv8 = rv8.assign(coding_neighbor=rv8.cls.map(LAB))
assert rv8.coding_neighbor.notna().all(), rv8.cls.unique()
rv8 = rv8.rename(columns={"gene": "symbol_rv8", "coding_gene": "coding_neighbor_gene", "coding_variant": "coding_neighbor_variant",
                          "r2": "coding_neighbor_r2"})
rv8 = rv8[["ensembl", "symbol_rv8", "coding_neighbor", "coding_neighbor_gene", "coding_neighbor_variant", "coding_neighbor_r2"]]
assert len(rv8) == 436 and rv8.ensembl.is_unique

rv9 = pd.read_csv(H / "RV9-20261102T151012Z/out-20261102T151110Z/gene_states.tsv", sep="\t")
rv9 = rv9.rename(columns={"gene": "symbol_rv9", "state": "replication_state", "discovery": "replication_discovery_gwas",
                          "replication": "replication_partner_gwas", "disc_pp4": "replication_discovery_pp4",
                          "rep_pp3": "replication_partner_pp3", "rep_pp4": "replication_partner_pp4"})
rv9 = rv9.drop(columns="direct")
assert len(rv9) == 436 and rv9.ensembl.is_unique

gen = rv8.merge(rv9, on="ensembl", how="outer", validate="1:1")
assert len(gen) == 436
out = rv5.merge(gen, on="ensembl", how="outer", validate="1:1")
out.insert(1, "symbol", out.symbol_rv5.fillna(out.symbol_rv8).fillna(out.symbol_rv9))
out = out.drop(columns=["symbol_rv5", "symbol_rv8", "symbol_rv9"])
FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
tier = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t").set_index("study_name").tier
a = pd.read_csv(FM / "results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv", usecols=["gwas_name", "ensembl", "PP.H4.susie"])
best = a[a.gwas_name.map(tier) <= 2].groupby("ensembl")["PP.H4.susie"].max().rename("best_tier12_pp4")
out = out.merge(best, left_on="ensembl", right_index=True, how="left")
assert (out.loc[out.coding_neighbor.notna(), "best_tier12_pp4"] > 0.5).all()
# provisional: colocalized gene whose best tier 1-2 SuSiE PP.H4 is below 0.8 and that does not replicate
out["provisional_coloc"] = (out.replication_state.notna() & (out.replication_state != "replicated")
                            & (out.best_tier12_pp4 < 0.8))
out.to_csv(OUT / "catalog_rv_annotations.tsv", sep="\t", index=False)

with open(OUT / "summary.txt", "w") as fh:
    fh.write(f"rows {len(out)}; bulk calls {out.consistency_class.notna().sum()}; colocalized genes {out.coding_neighbor.notna().sum()}; both {(out.consistency_class.notna() & out.coding_neighbor.notna()).sum()}\n\n")
    for c in ["consistency_class", "coding_neighbor", "replication_state", "provisional_coloc"]:
        fh.write(out[c].value_counts(dropna=False).to_string() + "\n\n")
print("done")
