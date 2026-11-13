# Verify program_uid -> name and example-gene membership for the three Figure 4D programs.
import sys, pandas as pd
P = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/passport_program_membership.parquet"
G = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/figures/example_genes.tsv"
out = sys.argv[1]
m = pd.read_parquet(P)
print(m.columns.tolist()); print(m.head(3).to_string())
genes = pd.read_csv(G, sep="\t")
pref = ["553ce773", "e65842a3", "07f7af82", "f05c535a", "48f39dd4"]
uid_col = [c for c in m.columns if "uid" in c.lower()][0]
sel = m[m[uid_col].astype(str).str.split("_").str[-1].str[:8].isin(pref)]
print(sel.groupby(uid_col).size())
cols = list(m.columns)
print(sel.drop_duplicates(subset=[c for c in cols if "uid" in c.lower() or "name" in c.lower() or "label" in c.lower()]).iloc[:, :8].to_string())
sym_cols = [c for c in cols if c.lower() in ("symbol", "gene_symbol", "gene_name", "gene")]
ens_cols = [c for c in cols if "ensembl" in c.lower()]
hit = sel[sel[ens_cols[0]].astype(str).str.split(".").str[0].isin(genes.ensembl_id)] if ens_cols else sel[sel[sym_cols[0]].isin(genes.symbol)]
print(hit.to_string())
hit.to_csv(out, sep="\t", index=False)
