import pandas as pd, numpy as np
x=pd.read_pickle("t12.pkl")
s=x[x["PP.H4.susie"].notna()].copy()
s["eur"]=s.reg_anc=="EUR"
g=s.groupby(["ensembl","eur"])["PP.H4.susie"].max().unstack()
g.columns=["nonEUR","EUR"] if list(g.columns)==[False,True] else g.columns
pos=g[(g.max(axis=1)>0.5)]
eur_only=pos[(pos.EUR>0.5)&~(pos.nonEUR>0.5)]; shared=pos[(pos.EUR>0.5)&(pos.nonEUR>0.5)]; non=pos[~(pos.EUR>0.5)&(pos.nonEUR>0.5)]
print("EUR only",len(eur_only)," tested nonEUR",eur_only.nonEUR.notna().sum()," untested",eur_only.nonEUR.isna().sum())
print("shared",len(shared),"nonEUR only",len(non)," of which EUR untested",non.EUR.isna().sum())
sym=x.drop_duplicates("ensembl").set_index("ensembl").gene
print(sym.reindex(non[non.EUR.isna()].index).tolist())
# ABF in nonEUR for EUR-only
a=x[x.reg_anc!="EUR"].groupby("ensembl")["PP.H4.abf"].max()
eo=eur_only.join(a.rename("abf_non"))
print("EUR-only with nonEUR ABF>0.5: tested",((eo.nonEUR.notna())&(eo.abf_non>0.5)).sum()," untested",((eo.nonEUR.isna())&(eo.abf_non>0.5)).sum())
# EUR ABF for the 3
ae=x[x.reg_anc=="EUR"].groupby("ensembl")["PP.H4.abf"].max()
print(ae.reindex(non[non.EUR.isna()].index).round(3).tolist())
