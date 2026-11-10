# Fig 2B: which fine-mapping run supplies the x (PIP) value of each plotted gene
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
b=pd.read_csv(R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/2B.tsv",sep="\t")
x=pd.read_csv(R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/out/fig2b_same_row_x_by_rule.tsv",sep="\t")
m=b.merge(x[["gene","y_gwas","sources","x_max","x_min"]],left_on=["gene","gwas"],right_on=["gene","y_gwas"],how="left")
print("rows",len(m),"unmatched",m.sources.isna().sum())
print("all 866 by source:",m.sources.value_counts().to_dict())
print("coloured 327 by source:",m[m.shown_rule].sources.value_counts().to_dict())
hi=m.susie_pp4>0.5
print("PP.H4>0.5 plotted:",hi.sum(),"by source",m[hi].sources.value_counts().to_dict())
print("PP.H4<=0.5 plotted:",(~hi).sum(),"by source",m[~hi].sources.value_counts().to_dict())
print("coloc_targeted-only rows by GWAS:",m[m.sources=="coloc_targeted"].gwas.value_counts().head(10).to_dict())
print("x plotted = x_max?", (abs(m.susie_pip_max-m.x_max)<1e-9).mean(), " x_min col equal?", (abs(m.x_min_x-m.x_min_y)<1e-9).mean() if 'x_min_x' in m else '')
ct=pd.read_csv(R+"GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv",usecols=["study","locus"])
print("coloc_targeted studies:",ct.study.value_counts().to_dict(), "loci", ct.groupby('study').locus.nunique().to_dict())
