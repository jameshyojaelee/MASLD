# Results-text counts for Figure 2: MVP-only support, direct-trait support, per-pair rates
import pandas as pd, numpy as np
x=pd.read_csv("t12_coloc_r2.tsv.gz",sep="\t",low_memory=False)
x["key"]=x.ensembl.fillna(x.gene); x["s"]=pd.to_numeric(x["PP.H4.susie"],errors="coerce")
mvp={"MVP_ALT_EUR","MVP_AST_EUR","MVP_NAFLD_EUR"}
x["cls"]=np.where(x.trait.isin({"NAFLD","NASH","PDFF"}),"direct","enzyme"); pos=x[x.s>0.5]; sup=pos.groupby("key").gwas_name.apply(set)
print("multi",len(sup),"MVP-EUR-only support",sum(v<=mvp for v in sup))
direct={"NAFLD","NASH","PDFF"}; x["cls"]=np.where(x.trait.isin(direct),"direct","enzyme")
pc=pos.groupby("key").apply(lambda d:set(np.where(d.trait.isin(direct),"direct","enzyme")))
nd=sum("direct" in v for v in pc); both=sum(v=={"direct","enzyme"} for v in pc); print("direct-supported",nd,"both",both,"direct only",nd-both)
dsup=pos[pos.trait.isin(direct)].groupby("key").gwas_name.apply(set); print("MVP_NAFLD_EUR sole direct support",sum(v=={"MVP_NAFLD_EUR"} for v in dsup))
fin=x[x.s.notna()]
print("finite pairs by class",fin.groupby("cls").size().to_dict(),"positives",fin[fin.s>0.5].groupby("cls").size().to_dict())
tested=fin.groupby("key").cls.apply(set); tb=set(tested[tested.apply(lambda v:v=={"direct","enzyme"})].index); print("genes SuSiE-tested in both classes",len(tb))
pb=pos[pos.key.isin(tb)].groupby("key").cls.apply(set)
print(" of these enzyme-supported",sum("enzyme" in v for v in pb),"direct-supported",sum("direct" in v for v in pb))
dd=pos[pos.trait.isin(direct)]; print("direct-only genes SuSiE-tested in enzyme GWAS:", sum(1 for k,v in pc.items() if v=={"direct"} and k in set(fin[fin.cls=='enzyme'].key)))
best=dd.sort_values("s",ascending=False).drop_duplicates("key"); print("MVP_NAFLD_EUR strongest direct",(best.gwas_name=="MVP_NAFLD_EUR").sum())
