# Fig 2E and 2F re-derived from COLOC r2 tier-1/2 rows
import pandas as pd, numpy as np
x=pd.read_csv("t12_coloc_r2.tsv.gz",sep="\t",low_memory=False)
x["key"]=x.ensembl.fillna(x.gene)
x["s"]=pd.to_numeric(x["PP.H4.susie"],errors="coerce"); x["a"]=pd.to_numeric(x["PP.H4.abf"],errors="coerce")
x["eur"]=x.reg_anc.eq("EUR")
g=x.groupby("key")
ms=g.s.max(); multi=set(ms[ms>0.5].index); print("multi-signal genes",len(multi))
pos=x[x.s>0.5]; e=set(pos[pos.eur].key); n=set(pos[~pos.eur].key)
eo=multi&e-n; sh=multi&e&n; no=multi&n-e
print("EUR only",len(eo),"shared",len(sh),"nonEUR only",len(no))
tested_ne=set(x[(~x.eur)&x.s.notna()].key)
print("EUR-only tested in nonEUR",len(eo&tested_ne),"untested",len(eo-tested_ne))
abf_ne=set(x[(~x.eur)&(x.a>0.5)].key)
print("ABF nonEUR>0.5 within tested",len(eo&tested_ne&abf_ne),"within untested",len(eo-tested_ne-set()) and len((eo-tested_ne)&abf_ne))
# 2F
x["pp"]=np.fmax(x.s.fillna(-1),x.a.fillna(-1))
ga=x.groupby(["key","reg_anc"]).agg(pp=("pp","max"),s=("s","max"),sym=("gene","first")).reset_index()
ga["sup"]=ga.pp>0.5; ga["ssup"]=ga.s>0.5
cnt=ga.groupby("key").agg(n=("sup","sum"),ns=("ssup","sum"),sym=("sym","first"))
f=cnt[cnt.n>=4].sort_values(["ns","sym"],ascending=[False,True]); print("2F genes:",len(f)); print(f.to_string())
print("all five:",(cnt.n==5).sum())
print(ga[ga.key.isin(f.index)].pivot(index="sym",columns="reg_anc",values="pp").round(2).to_string())
print(ga[ga.key.isin(f.index)].pivot(index="sym",columns="reg_anc",values="s").round(2).to_string())
