import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
ct=pd.read_pickle("ct.pkl")
b=pd.read_csv(R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/2B.tsv",sep="\t")
b["chr"]=b.lead_variant.str.split(":").str[0].astype(int); b["pos"]=b.lead_variant.str.split(":").str[1].astype(int)
g=ct.groupby(["study","chromosome","position"]).agg(ct_pip_max=("susie_pip","max"),ct_pip_min=("susie_pip","min"),ct_rec_max=("recommended_pip","max"),ct_rows=("susie_pip","size")).reset_index()
m=b.merge(g,left_on=["gwas","chr","pos"],right_on=["study","chromosome","position"],how="left")
inct=m.ct_rows.notna()
print("2B rows",len(b),"lead present in coloc-targeted run (same study):",inct.sum())
eq=(abs(m.susie_pip_max-m.ct_pip_max)<1e-9)
print("x(susie_pip_max)==ct susie_pip max:",eq.sum(), " among present:", eq[inct].sum())
print("x>ct max (another run supplies larger PIP):",(m.susie_pip_max>m.ct_pip_max+1e-9).sum())
hi=m.susie_pp4>0.5
print("PP.H4>0.5 points:",hi.sum()," of which lead in ct:",(hi&inct).sum(), " x equals ct:",(hi&eq).sum())
print("coloured (shown_rule) points:",m.shown_rule.sum()," lead in ct:",(m.shown_rule&inct).sum(), " x equals ct:",(m.shown_rule&eq).sum())
print("ct-present by gwas:",m[inct].gwas.value_counts().to_dict())
for gname in ["RORA","FABP1","HSD17B13"]:
    print(m[m.gene==gname][["gene","gwas","lead_variant","susie_pp4","susie_pip_max","x_min","ct_pip_max","ct_rows","shown_rule","class"]].to_string())
m.to_pickle("b_m.pkl")
