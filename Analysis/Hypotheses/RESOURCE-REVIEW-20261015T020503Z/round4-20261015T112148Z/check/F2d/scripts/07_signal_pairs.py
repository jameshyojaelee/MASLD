import pandas as pd
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
R=ROOT+"/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
tier=pd.read_csv(ROOT+"/GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
main=tier[tier.placement=="main"].set_index("study_name")
sp=pd.read_csv(R+"/susie_coloc_signal_pairs.tsv.gz",sep="\t")
print("pairs",len(sp),"eligible",sp.eligible.value_counts().to_dict())
sp=sp[sp.gwas.isin(main.index)].copy()
sp["tier"]=sp.gwas.map(main.tier)
el=sp[sp.eligible.astype(str).str.upper()=="TRUE"]
for t in (1,2):
    a=sp[sp.tier==t]; e=el[el.tier==t]
    print(f"tier {t}: all pairs {len(a)} eligible {len(e)} hits_eligible {(e['PP.H4']>0.5).sum()} rate {(e['PP.H4']>0.5).mean():.4f} hits_all {(a['PP.H4']>0.5).sum()} studies {e.gwas.nunique()} gene-study rows {e.groupby(['gwas','gene']).ngroups}")
    print("   eligible pairs by study:", e.groupby('gwas').size().to_dict())
    print("   hits by study:", e[e['PP.H4']>0.5].groupby('gwas').size().to_dict())
