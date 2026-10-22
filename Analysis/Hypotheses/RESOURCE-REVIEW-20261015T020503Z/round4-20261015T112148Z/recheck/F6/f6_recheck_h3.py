import json, sys
import numpy as np, pandas as pd
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CAT=ROOT+"/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
COLOC=ROOT+"/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
OUT=sys.argv[1]; R={}
b=lambda s: s.astype(str).str.split(".").str[0]
e=pd.read_parquet(CAT+"/passport_evidence_long.parquet"); gen=e[e.evidence_domain=="genetics"].set_index("ensembl_id")
tier=pd.read_csv(ROOT+"/GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t").set_index("study_name").tier
a=pd.read_csv(COLOC+"/susie_coloc_all_gwas.csv",usecols=["gwas_name","ensembl","PP.H3.susie","PP.H4.susie","PP.H3.abf","PP.H4.abf","n_cs_pairs","n_cs_pairs_eligible","method","eligibility_source"])
a["ens"]=b(a.ensembl); a["tier"]=a.gwas_name.map(tier)
a12=a[a.tier.isin([1,2])]
R["method_t12"]={str(k):int(v) for k,v in a12.method.value_counts(dropna=False).items()}
R["eligibility_source_t12"]={str(k):int(v) for k,v in a12.eligibility_source.value_counts(dropna=False).items()}
sp=pd.read_csv(COLOC+"/susie_coloc_signal_pairs.tsv.gz",sep="\t"); sp["ens"]=b(sp.ensembl)
# compare all_gwas PP.H3.susie with signal-pair level for rerun gene-studies
m=sp.groupby(["ens","gwas"]).agg(maxh3=("PP.H3","max"),maxh4=("PP.H4","max")).reset_index()
h4top=sp.sort_values("PP.H4",ascending=False).drop_duplicates(["ens","gwas"])[["ens","gwas","PP.H3"]].rename(columns={"PP.H3":"h3_at_maxh4"})
m=m.merge(h4top,on=["ens","gwas"]).merge(a[["ens","gwas_name","PP.H3.susie","PP.H4.susie"]].rename(columns={"gwas_name":"gwas"}),on=["ens","gwas"],how="left")
mm=m.dropna(subset=["PP.H3.susie"])
R["n_compare"]=len(mm)
R["allg_h3_eq_max_h3"]=int(np.isclose(mm["PP.H3.susie"],mm.maxh3,atol=1e-9).sum())
R["allg_h3_eq_h3_at_maxh4"]=int(np.isclose(mm["PP.H3.susie"],mm.h3_at_maxh4,atol=1e-9).sum())
R["allg_h4_eq_max_h4"]=int(np.isclose(mm["PP.H4.susie"],mm.maxh4,atol=1e-9).sum())
s=a12.dropna(subset=["PP.H4.susie"])
gg=s.groupby("ens").agg(h3=("PP.H3.susie","max"),h4=("PP.H4.susie","max"))
gg["call"]=[gen.call_state.get(x,"absent") for x in gg.index]
for thr in [0.8,0.9]:
    h=gg[(gg.h3>thr)&(gg.h4<=0.5)]
    R[f"allg_h3gt{thr}_h4le05_by_call"]={str(k):int(v) for k,v in h.call.value_counts().items()}
R["allg_genes_by_call"]={str(k):int(v) for k,v in gg.call.value_counts().items()}
# H3+H4 > 0.8 genes with susie (both signals present) to show denominator
open(OUT+"/recheck_h3.json","w").write(json.dumps(R,indent=1))
print("h3 done")
