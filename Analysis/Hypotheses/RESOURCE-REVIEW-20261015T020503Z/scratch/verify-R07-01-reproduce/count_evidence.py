import pandas as pd
C="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/"
ev=pd.read_parquet(C+"passport_evidence_long.parquet")
print(ev.shape); print(ev.columns.tolist())
for c in ev.columns:
    if ev[c].nunique()<30: print(c, ev[c].value_counts(dropna=False).to_dict())
cov=pd.read_parquet(C+"passport_coverage_long.parquet")
print("coverage",cov.shape, cov.columns.tolist())
for c in cov.columns:
    if cov[c].nunique()<30: print(" ",c, cov[c].value_counts(dropna=False).to_dict())
pe=pd.read_parquet(C+"portal_export/passport_evidence.parquet")
print("portal evidence",pe.shape)
for c in pe.columns:
    if pe[c].nunique()<30: print(" ",c, pe[c].value_counts(dropna=False).to_dict())
