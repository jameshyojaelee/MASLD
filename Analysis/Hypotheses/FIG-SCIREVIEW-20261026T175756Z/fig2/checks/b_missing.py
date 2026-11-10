# Fig 2B: why are multi-signal genes absent from the scatter?
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
x=pd.read_csv(R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/out/fig2b_same_row_x_by_rule.tsv",sep="\t")
b=pd.read_csv(R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/2B.tsv",sep="\t")
hi=x[x.pp4>0.5]; print("genes with max SuSiE PP.H4>0.5 (by symbol):",len(hi))
nd=hi[~hi.gene.isin(b.gene)]
mvp=nd.y_gwas.isin(["MVP_ALT_EUR","MVP_AST_EUR","MVP_NAFLD_EUR"])
print("not drawn",len(nd),"of which max in MVP EUR",mvp.sum(),"other",(~mvp).sum())
print("other not-drawn by y_gwas:",nd[~mvp].y_gwas.value_counts().to_dict())
print("other not-drawn: any fine-mapping row?",nd[~mvp].n_rows.notna().sum(),"finite pip",np.isfinite(nd[~mvp].x_max).sum())
