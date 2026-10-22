import sys, numpy as np, pandas as pd
from scipy import stats
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; O=sys.argv[1]
INST=f"{R}/figures/main/fig4_singlecell_programs/source_tables/corrected_774/stage_continuum_comparison"
ATLAS=f"{R}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/overlays/all_node_stage_continuum_overlays.tsv.gz"
feat=pd.read_csv(f"{O}/feature_betas_recheck.tsv.gz",sep="\t")
print(feat.groupby(['roster']).agg(n1=('n1','max'),n2=('n2','max')))
ip=pd.read_csv(f"{INST}/pathway_stage_vs_continuum.tsv.gz",sep="\t"); ih=pd.read_csv(f"{INST}/hotspot_stage_vs_continuum.tsv",sep="\t")
testable=set(ip.loc[ip.estimable.fillna(False).astype(bool),"set_id"])|set(ih.loc[ih.estimable.fillna(False).astype(bool),"program_uid"])
# reproduction of installed feature betas on 435
inst=pd.concat([ip[['set_id','beta','continuum_beta']].rename(columns={'set_id':'feature_id'}),ih[['program_uid','beta','continuum_beta']].rename(columns={'program_uid':'feature_id'})])
x=feat[feat.roster=='as_published_435'].merge(inst,on='feature_id')
print('435 repro: n',len(x),'max|stage diff|',np.nanmax(abs(x.s_beta-x.beta)),'max|cont diff|',np.nanmax(abs(x.c_beta-x.continuum_beta)))
rows=[]
for (rn,coll),g in feat[feat.feature_id.isin(testable)].groupby(['roster','collection']):
    g=g.dropna(subset=['s_beta','c_beta'])
    rows.append((rn,coll,len(g),round(stats.spearmanr(g.s_beta,g.c_beta).correlation,4),int((np.sign(g.s_beta)==np.sign(g.c_beta)).sum())))
lay=pd.DataFrame(rows,columns=['roster','collection','n','rho','same']); lay.to_csv(f"{O}/layer_rho_recheck.tsv",sep="\t",index=False); print(lay.to_string())
atlas=pd.read_csv(ATLAS,sep="\t",usecols=["feature_id","collection","community_id","component_role"]); atlas=atlas[atlas.component_role=="main_system_component"]
syss=[]
for rn in feat.roster.unique():
    f=feat[feat.roster==rn].merge(atlas,on=["feature_id","collection"])
    fam=f.groupby(["community_id","collection"]).agg(s=("s_beta","median"),c=("c_beta","median")).reset_index()
    sy=fam.groupby("community_id").agg(stage_beta=("s","mean"),continuum_beta=("c","mean")).reset_index(); sy['roster']=rn; syss.append(sy)
    rho=stats.spearmanr(sy.stage_beta,sy.continuum_beta).correlation; same=int((np.sign(sy.stage_beta)==np.sign(sy.continuum_beta)).sum())
    q98=np.quantile(np.abs(f.c_beta.dropna()),0.98,method="median_unbiased")
    s01=sy.set_index('community_id').loc['system_01']
    print(f"SYSTEM {rn}: n_sys={len(sy)} n_nodes={len(f)} rho={rho:.6f} same={same} q98={q98:.6f} S01 stage={s01.stage_beta:.6f} cont={s01.continuum_beta:.6f}")
sy=pd.concat(syss); sy.to_csv(f"{O}/system_betas_recheck.tsv",sep="\t",index=False)
cand=pd.read_csv(f"{R}/figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/fig4h_disease_only_source.tsv",sep="\t")
y=sy[sy.roster=='disease_only_367'].merge(cand,on='community_id',suffixes=('_mine','_cand'))
print('vs candidate 367: max|stage diff|',np.abs(y.stage_beta_mine-y.stage_beta_cand).max(),'max|cont diff|',np.abs(y.continuum_beta_mine-y.continuum_beta_cand).max())
print('candidate printed rho',stats.spearmanr(cand.stage_beta,cand.continuum_beta).correlation,'same',int((np.sign(cand.stage_beta)==np.sign(cand.continuum_beta)).sum()), 'n',len(cand))
