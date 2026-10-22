import pandas as pd
C="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/"
for f in ["passport_program_index.parquet","passport_program_membership.parquet","passport_program_context.parquet","passport_gene_index.parquet"]:
    d=pd.read_parquet(C+f); print("==",f,d.shape); print(d.columns.tolist())
    for c in d.columns:
        try:
            if d[c].nunique()<15: print("  ",c,d[c].value_counts(dropna=False).to_dict())
        except TypeError: pass
m=pd.read_parquet(C+"passport_program_membership.parquet")
gcol=[c for c in m.columns if c in("symbol","gene_symbol","gene")]
print("membership gene col",gcol, m[gcol[0]].nunique() if gcol else None)
if gcol:
    print("GNMT in membership", (m[gcol[0]]=="GNMT").sum(), "HKDC1", (m[gcol[0]]=="HKDC1").sum())
