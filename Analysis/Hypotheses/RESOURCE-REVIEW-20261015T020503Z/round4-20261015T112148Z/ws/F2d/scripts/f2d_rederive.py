# F2d re-derivation: Fig 2D trait-subclass split (R14-02), cluster concentration (R14-06),
# and colocalization yield per eligible pair (R14-03). Seed 20260923, 10,000 physical-cluster draws.
import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, json, sys
from scipy.stats import fisher_exact
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT=sys.argv[1]
d=pd.read_csv(ROOT+'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv',sep='\t')
assert len(d)==6415 and d.pip_sum_gate_passed.all()
d['coding']=d.protein_altering_pip_mass+d.canonical_splice_pip_mass
d=d.sort_values(['chrom','pos']).reset_index(drop=True)
cid=[];c=-1;prev=None
for ch,p in zip(d.chrom,d.pos):
    if prev is None or ch!=prev[0] or p-prev[1]>1_000_000: c+=1
    cid.append(c); prev=(ch,p)
d['cluster']=cid
direct=d.trait_scope=='tier1_direct_masld_pdff'
d['cls']=np.where(~direct,'ENZ',np.where(d.trait=='PDFF','PDFF','NAFLD_NASH'))
def colmean(mask):
    return d[mask].groupby('cluster').coding.mean()
cols={
 'ENZ':colmean(~direct),'DIRECT':colmean(direct),
 'PDFF':colmean(direct&(d.trait=='PDFF')),'NAFLD_NASH':colmean(direct&(d.trait!='PDFF')),
 'NAFLD_NASH_EUR':colmean(direct&(d.trait!='PDFF')&(d.ancestry=='EUR')),
 'NAFLD_AMR':colmean(direct&(d.ancestry=='AMR')),
 'ENZ_EUR':colmean(~direct&(d.ancestry=='EUR')),'ENZ_AMR':colmean(~direct&(d.ancestry=='AMR')),
}
tab=pd.DataFrame(cols).reindex(range(d.cluster.max()+1))
ncl=len(tab); A=tab.values; names=list(tab.columns)
rng=np.random.default_rng(20260923); B=10000
idx=rng.integers(0,ncl,(B,ncl))
# column means over resampled clusters ignoring NaN
boot={}
for j,n in enumerate(names):
    v=A[:,j][idx]; boot[n]=np.nanmean(v,axis=1)
obs={n:np.nanmean(A[:,j]) for j,n in enumerate(names)}
rows=[]
for num,den in [('DIRECT','ENZ'),('PDFF','ENZ'),('NAFLD_NASH','ENZ'),('NAFLD_NASH_EUR','ENZ_EUR'),('PDFF','ENZ_EUR'),('NAFLD_AMR','ENZ_AMR')]:
    r=boot[num]/boot[den]; df=boot[num]-boot[den]; ok=np.isfinite(r)&np.isfinite(df)
    rows.append(dict(direct_subset=num,enzyme_subset=den,n_clusters_direct=int(tab[num].notna().sum()),n_clusters_enzyme=int(tab[den].notna().sum()),
      n_cs_direct=None,frac_direct=obs[num],frac_enzyme=obs[den],ratio=obs[num]/obs[den],ratio_lo95=np.percentile(r[ok],2.5),ratio_hi95=np.percentile(r[ok],97.5),
      difference=obs[num]-obs[den],diff_lo95=np.percentile(df[ok],2.5),diff_hi95=np.percentile(df[ok],97.5),n_finite=int(ok.sum())))
res=pd.DataFrame(rows)
ncs={'DIRECT':direct.sum(),'PDFF':(direct&(d.trait=='PDFF')).sum(),'NAFLD_NASH':(direct&(d.trait!='PDFF')).sum(),
     'NAFLD_NASH_EUR':(direct&(d.trait!='PDFF')&(d.ancestry=='EUR')).sum(),'NAFLD_AMR':(direct&(d.ancestry=='AMR')).sum()}
res['n_cs_direct']=res.direct_subset.map(ncs)
res.to_csv(OUT+'/fig2d_trait_subclass_split.tsv',sep='\t',index=False)
print(res.to_string())
print('direct study composition'); print(d[direct].groupby(['trait','ancestry']).agg(n_cs=('coding','size'),studies=('study','nunique'),clusters=('cluster','nunique')).to_string())
# R14-06 cluster concentration of direct coding fraction
dm=cols['DIRECT'].sort_values(ascending=False)
share=(dm/dm.sum())
lab=d[direct].groupby('cluster').agg(chrom=('chrom','first'),pos_min=('pos','min'),n_cs=('coding','size'),studies=('study',lambda s:';'.join(sorted(set(s)))))
conc=pd.DataFrame({'cluster_coding_mean':dm,'share_of_direct_sum':share,'cum_share':share.cumsum()}).join(lab)
enz=obs['ENZ']; drop=[]
for k in range(0,6):
    keep=dm.index[k:]; drop.append(dict(top_k_removed=k,direct_frac=dm.loc[keep].sum()/19,  # note: denominator fixed below
       ))
conc.to_csv(OUT+'/fig2d_direct_cluster_concentration.tsv',sep='\t')
print(conc.head(8).to_string())
for k in range(6):
    keep=dm.index[k:]; f=dm.loc[keep].mean()
    print('drop top',k,'direct frac',round(f,4),'ratio',round(f/enz,3))
# PNPLA3 cluster
p=d[direct&(d.chrom==22)&(d.pos.between(43_000_000,45_500_000))]
print('chr22 PNPLA3 region direct cs',len(p),'cluster ids',p.cluster.unique(),'coding mean',round(p.groupby('cluster').coding.mean().iloc[0],4))
# R14-03 colocalization yield per eligible pair
R=ROOT+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tier=pd.read_csv(ROOT+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t'); t=dict(zip(tier.study_name,tier.tier))
sp=pd.read_csv(R+'susie_coloc_signal_pairs.tsv.gz',sep='\t'); sp['tier']=sp.gwas.map(t); sp=sp[sp.tier.isin([1,2])]
el=sp.eligible if sp.eligible.dtype==bool else sp.eligible.astype(str).str.lower().isin(['true','1'])
sp=sp[el]; sp['hit']=sp['PP.H4']>0.5
yt=sp.groupby('tier').agg(eligible_pairs=('hit','size'),hit_pairs=('hit','sum'),studies_with_pairs=('gwas','nunique'))
yt['gated_cs']=[int((d.tier==k).sum()) for k in yt.index]
hit=sp[sp.hit]; g1=set(hit[hit.tier==1].ensembl); g2=set(hit[hit.tier==2].ensembl)
gl=pd.read_csv(R+'gene_level_coloc_tier12.csv'); h=gl[gl.coloc_best_susie_pp4>0.5]
tabf=[[int(yt.loc[1,'hit_pairs']),int(yt.loc[1,'eligible_pairs']-yt.loc[1,'hit_pairs'])],[int(yt.loc[2,'hit_pairs']),int(yt.loc[2,'eligible_pairs']-yt.loc[2,'hit_pairs'])]]
fo=fisher_exact(tabf)
summ=dict(yield_by_tier=yt.reset_index().to_dict('records'),genes_tier1_hit=len(g1),genes_tier2_hit=len(g2),genes_both=len(g1&g2),genes_union=len(g1|g2),
  genes_tier1_only=len(g1-g2),driving_tier=h.driving_tier.value_counts().to_dict(),n_susie_gt05=len(h),fisher_OR=fo[0],fisher_p_nominal=fo[1])
print(json.dumps(summ,indent=1,default=str))
json.dump(summ,open(OUT+'/coloc_yield_by_tier.json','w'),indent=1,default=str)
