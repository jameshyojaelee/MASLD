import pandas as pd, sys
pd.set_option("display.width",250); pd.set_option("display.max_columns",80); pd.set_option("display.max_colwidth",120)
C="RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/"
g=pd.read_parquet(C+"passport_gene_index.parquet"); e=pd.read_parquet(C+"passport_evidence_long.parquet"); cov=pd.read_parquet(C+"passport_coverage_long.parquet")
for n,d in [("gene",g),("evidence",e),("coverage",cov)]:
    print("==",n,d.shape); print(d.dtypes.to_string())
for c in g.columns:
    if g[c].nunique()<30: print(g[c].value_counts(dropna=False).to_string()); print()
for c in e.columns:
    if e[c].nunique()<30: print(e[c].value_counts(dropna=False).to_string()); print()
for c in cov.columns:
    if cov[c].nunique()<30: print(cov[c].value_counts(dropna=False).to_string()); print()
h=["MTARC1","HKDC1","GLP1R","THRB"]
print(g[g.symbol.isin(h)].T.to_string())
ids=g[g.symbol.isin(h)].passport_id
print(e[e.passport_id.isin(ids)].T.to_string())
print(cov[cov.passport_id.isin(ids)].T.to_string())
