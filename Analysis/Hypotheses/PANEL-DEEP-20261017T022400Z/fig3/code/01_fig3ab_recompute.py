"""Recompute Fig 3A (cohort composition, metadata availability) and 3B (fibrosis x NAS grid)
from the adopted 774-unit fit metadata (meta_matched.rds dumped by 00_dump_meta.R)."""
import pandas as pd, json
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
m=pd.read_csv(O+"meta_774.tsv",sep="\t",keep_default_na=False,na_values=["NA",""])
# 3A
rows=[]
for d,g in m.groupby('dataset'):
    n=len(g)
    def state(col, exact=True):
        k=g[col].notna().sum()
        return f"{k}/{n}"
    rows.append(dict(dataset=d,n=n,control=(g.group_binary=='Control').sum(),disease=(g.group_binary=='Disease').sum(),
        sex_annotated=(g.sex_source=='annotated').sum(), sex_inferred_nonNA=g.inferred_sex.notna().sum(),
        nas_nonNA=g.nas_score.notna().sum(), fib_nonNA=g.fibrosis_stage.notna().sum(),
        fib_binned=g.condition.str.contains('F0F1|F3F4|F2F3').sum()))
a=pd.DataFrame(rows).sort_values('n',ascending=False); a.to_csv(O+"fig3a_recomputed.tsv",sep="\t",index=False); print(a.to_string())
# 3B : exact fibrosis (not GSE213621 bins), NAS present
s=m[m.dataset.isin(['GSE130970','GSE135251','GSE162694']) & m.fibrosis_stage.notna()]
print("staged units (exact F0-F4):",len(s), s.fibrosis_stage.value_counts().sort_index().to_dict())
b=s[s.nas_score.notna()]
grid=pd.crosstab(b.nas_score.astype(int),b.fibrosis_stage.astype(int),margins=True)
print(grid); grid.to_csv(O+"fig3b_recomputed_grid.tsv",sep="\t")
print("staged without NAS:",s.nas_score.isna().sum(), s[s.nas_score.isna()].groupby(['dataset','fibrosis_stage']).size().to_dict())
print("controls in grid:", b.groupby(['dataset','group_binary']).size().to_dict())
print("controls by cell:", b[b.group_binary=='Control'].groupby(['nas_score','fibrosis_stage']).size().to_dict())
