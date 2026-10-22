import sys; sys.path.insert(0, sys.argv[1]); from refit import *
full=family(d)
ref=pd.read_csv(REF,sep="\t")
m=full.merge(ref[["cell_type","module","beta","qvalue","hc3_qvalue","n_donors"]],on=["cell_type","module"])
print("max|beta diff|",np.nanmax(abs(m.beta_x-m.beta_y)),"max|q diff|",np.nanmax(abs(m.q-m.qvalue)),"max|hc3 q diff|",np.nanmax(abs(m.q_hc3-m.hc3_qvalue)),"n mismatch",(m.n!=m.n_donors).sum())
duc=(full.cell_type=="hepatocytes")&(full.module==20)
print("FULL Ductular:",full[duc].to_string(index=False)); print("FULL n q<0.05:",(full.q<0.05).sum(),"n hc3 q<0.05:",(full.q_hc3<0.05).sum())
hep=d[(d.cell_type=="hepatocytes")&(d.module==20)]
print("Ductular donors by dataset x stage:"); print(pd.crosstab(hep.dataset,hep.disease_stage_coarse))
print(hep[hep.dataset=="GSE189600"][["donor","disease_stage_coarse","score"]].to_string(index=False))
# drop GSE189600
nog=family(d[d.dataset!="GSE189600"])
print("DROP GSE189600 Ductular:",nog[duc].to_string(index=False)); print("DROP GSE189600 n q<0.05:",(nog.q<0.05).sum(),"hc3:",(nog.q_hc3<0.05).sum())
print(nog.sort_values("q").head(4).to_string(index=False))
# single-donor deletions (donors in Ductular fit)
res=[]
for don in sorted(hep.donor.unique()):
    f=family(d[d.donor!=don]); r=f[duc].iloc[0]
    res.append((don,hep.set_index("donor").loc[don,"dataset"],hep.set_index("donor").loc[don,"disease_stage_coarse"],r.beta,r.q,r.q_hc3,(f.q<0.05).sum()))
L=pd.DataFrame(res,columns=["donor","dataset","stage","beta","q","q_hc3","n_sel"])
L.to_csv(sys.argv[1]+"/loo_donor_ductular.tsv",sep="\t",index=False)
print("LOO donors:",len(L),"q>=0.05:",(L.q>=0.05).sum(),"q range",L.q.min(),L.q.max())
print(L[L.q>=0.05].to_string(index=False))
print(pd.crosstab(L[L.q>=0.05].dataset,L[L.q>=0.05].stage))
