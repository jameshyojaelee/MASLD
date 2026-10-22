"""Independent re-derivation of Fig 2D coding PIP shares (region- and instance-weighted)
from the release credible-set tables. Positions come from member rows, never from `locus`."""
import pandas as pd, numpy as np, sys
base='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/'
out=sys.argv[1]
arch=pd.read_csv(base+'credible_set_pip_architecture.tsv',sep='\t',dtype={'locus':str})
mem=pd.read_csv(base+'credible_set_members_annotated.tsv',sep='\t',dtype={'locus':str,'chromosome':str})
ctx=pd.read_csv(base+'credible_set_regulatory_context.tsv',sep='\t',dtype={'locus':str})
print('arch rows',len(arch), arch.pip_sum_gate_passed.value_counts().to_dict(), arch.trait_scope.value_counts().to_dict())
a=arch[(arch.pip_sum_gate_passed.astype(str).str.lower()=='true')&arch.trait_scope.isin(['tier1_direct_masld_pdff','tier2_liver_enzyme'])].copy()
print('gated tier1/2',len(a), a.trait_scope.value_counts().to_dict())
# sentinel: `locus` is chrom.position written as a number, so trailing zeros are lost.
# Restore by appending zeros until the position falls inside the member coordinate span.
bounds=mem.groupby('credible_set_uid').agg(mchrom=('chromosome','first'),mmin=('position','min'),mmax=('position','max'))
a=a.merge(bounds,left_on='credible_set_uid',right_index=True,how='left')
def restore(row):
    chrom,frag=row.locus.split('.')
    hits=[]
    for k in range(0,4):
        p=int(frag+'0'*k)
        if row.mmin<=p<=row.mmax: hits.append((k,p))
    if len(hits)!=1:
        # fall back: closest candidate to the member span
        cands=[int(frag+'0'*k) for k in range(0,4)]
        p=min(cands,key=lambda x:min(abs(x-row.mmin),abs(x-row.mmax)))
        return pd.Series({'chrom':chrom,'pos':p,'n_hits':len(hits),'k':-1})
    return pd.Series({'chrom':chrom,'pos':hits[0][1],'n_hits':1,'k':hits[0][0]})
rp=a.apply(restore,axis=1); a=pd.concat([a,rp],axis=1)
print('restore n_hits',a.n_hits.value_counts().to_dict(),'zeros added',a.k.value_counts().to_dict(),'chrom mismatch',int((a.chrom!=a.mchrom.astype(str)).sum()))
# naive parse (the trap)
a['naive_pos']=a.locus.str.split('.').str[1].astype(int)
print('naive differs',int((a.naive_pos!=a.pos).sum()))
# 1-Mb grouping: sort sentinels per chrom (all trait classes together), chain-link if gap <=1e6
def cluster(df,poscol):
    df=df.sort_values(['chrom',poscol]).copy(); cid=[];c=-1;prev_ch=None;prev=None
    for ch,p in zip(df.chrom,df[poscol]):
        if ch!=prev_ch or p-prev>1_000_000: c+=1
        cid.append(c); prev_ch=ch; prev=p
    df['cluster']=cid; return df
mass=['protein_altering_pip_mass','canonical_splice_pip_mass','synonymous_or_utr_pip_mass','other_noncoding_pip_mass']
a['tot']=a[mass].sum(axis=1)
print('per-instance total mass summary',a.tot.describe().round(4).to_dict())
a['coding']=a.protein_altering_pip_mass+a.canonical_splice_pip_mass
ctx_g=ctx[ctx.credible_set_uid.isin(a.credible_set_uid)][['credible_set_uid','promoter_proximal_pip_mass','abc_enhancer_pip_mass','lineage_accessible_pip_mass','no_deposited_context_pip_mass','unresolved_context_pip_mass']]
a=a.merge(ctx_g,on='credible_set_uid',how='left')
rows=[]
for grouping in ['joint_all_classes','within_class']:
  for scope,g0 in a.groupby('trait_scope'):
    pass
for scope in ['tier1_direct_masld_pdff','tier2_liver_enzyme']:
    g=a[a.trait_scope==scope]
    inst_pool=g.coding.sum()/g.tot.sum(); inst_mean=(g.coding/g.tot).mean()
    # region: cluster within this class's sentinels (one-Mb chaining)
    gc=cluster(g,'pos'); reg=gc.groupby('cluster').apply(lambda d:(d.coding/d.tot).mean()); 
    gn=cluster(g.assign(pos2=g.naive_pos),'pos2'); regn=gn.groupby('cluster').apply(lambda d:(d.coding/d.tot).mean())
    # noncoding context fractions (region-weighted), denominator = noncoding mass per instance
    nc=g.synonymous_or_utr_pip_mass+g.other_noncoding_pip_mass
    gc2=gc.assign(nc=(gc.synonymous_or_utr_pip_mass+gc.other_noncoding_pip_mass))
    gc2=gc2[gc2.nc>0]
    r_unres=gc2.groupby('cluster').apply(lambda d:(d.unresolved_context_pip_mass/d.nc).mean()).mean()
    r_nodep=gc2.groupby('cluster').apply(lambda d:(d.no_deposited_context_pip_mass/d.nc).mean()).mean()
    i_unres=(gc2.unresolved_context_pip_mass/gc2.nc).mean(); i_nodep=(gc2.no_deposited_context_pip_mass/gc2.nc).mean()
    rows.append(dict(trait_scope=scope,n_instances=len(g),n_regions=reg.size,n_regions_naive=regn.size,
        coding_instance_pooled=inst_pool,coding_instance_mean=inst_mean,coding_region_weighted=reg.mean(),coding_region_weighted_naive=regn.mean(),
        nc_n_instances=len(gc2),nc_n_regions=gc2.cluster.nunique(),
        unresolved_ctx_region=r_unres,no_deposited_ctx_region=r_nodep,sum_region=r_unres+r_nodep,
        unresolved_ctx_instance=i_unres,no_deposited_ctx_instance=i_nodep,
        n_studies=g.study.nunique(),largest_study_share=g.study.value_counts(normalize=True).iloc[0],largest_study=g.study.value_counts().index[0]))
res=pd.DataFrame(rows)
# joint clustering across both classes (shared physical clusters)
j=cluster(a,'pos'); print('joint physical clusters',j.cluster.nunique(), 'direct',j[j.trait_scope.str.startswith('tier1')].cluster.nunique(),'enzyme',j[j.trait_scope.str.startswith('tier2')].cluster.nunique(),
  'shared',len(set(j[j.trait_scope.str.startswith('tier1')].cluster)&set(j[j.trait_scope.str.startswith('tier2')].cluster)))
for scope in ['tier1_direct_masld_pdff','tier2_liver_enzyme']:
    d=j[j.trait_scope==scope]; print(scope,'joint-cluster region-weighted coding',d.groupby('cluster').apply(lambda x:(x.coding/x.tot).mean()).mean())
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
print(res.T)
r=res.set_index('trait_scope')
print('ratio region',r.coding_region_weighted.iloc[0]/r.coding_region_weighted.iloc[1],'ratio instance pooled',r.coding_instance_pooled.iloc[0]/r.coding_instance_pooled.iloc[1],'ratio instance mean',r.coding_instance_mean.iloc[0]/r.coding_instance_mean.iloc[1])
res.to_csv(out+'/r13_02_recomputed.tsv',sep='\t',index=False)

# ---- noncoding context, from members (noncoding DNA = other_noncoding plus non-synonymous UTR) ----
m=mem[mem.credible_set_uid.isin(a.credible_set_uid)].copy()
cons=m.consequence.astype(str)
is_utr=cons.str.contains(r'(?:^|,)[35]_prime_UTR_variant(?:$|,)'); is_syn=cons.str.contains(r'(?:^|,)synonymous_variant(?:$|,)')
nc=m[(m.consequence_category=='other_noncoding')|((m.consequence_category=='synonymous_or_utr')&is_utr&~is_syn)].copy()
fl=lambda c: nc[c].astype(str).str.lower()=='true'
pr=nc.promoter_context_status=='resolved'; lr=nc.lineage_accessibility_status=='resolved'
w=nc.normalized_susie_pip
nc['ncm']=w; nc['unres']=w*(~pr|~lr); nc['nodep']=w*(pr&lr&~fl('promoter_proximal')&~fl('abc_enhancer_overlap')&~fl('lineage_accessible'))
nc['unres_with_abc']=w*((~pr|~lr)&fl('abc_enhancer_overlap'))
ps=nc.groupby('credible_set_uid')[['ncm','unres','nodep','unres_with_abc']].sum().reset_index()
ps=ps.merge(j[['credible_set_uid','trait_scope','cluster']],on='credible_set_uid')
crow=[]
for scope,d in ps.groupby('trait_scope'):
    pl=d.groupby('cluster')[['ncm','unres','nodep','unres_with_abc']].mean()
    crow.append(dict(trait_scope=scope,n_instances=len(d),n_regions=len(pl),
      unresolved_region=pl.unres.sum()/pl.ncm.sum(),no_deposited_region=pl.nodep.sum()/pl.ncm.sum(),
      sum_region=(pl.unres.sum()+pl.nodep.sum())/pl.ncm.sum(),
      unresolved_instance=d.unres.sum()/d.ncm.sum(),no_deposited_instance=d.nodep.sum()/d.ncm.sum(),
      unresolved_with_abc_region=pl.unres_with_abc.sum()/pl.ncm.sum()))
cr=pd.DataFrame(crow); print(cr.T); cr.to_csv(out+'/r13_02_context_recomputed.tsv',sep='\t',index=False)
# status breakdown of unresolved mass
nc2=nc
u=nc2[(nc2.promoter_context_status!='resolved')|(nc2.lineage_accessibility_status!='resolved')]
print(u.groupby(['trait_scope','promoter_context_status','lineage_accessibility_status']).normalized_susie_pip.sum().round(2))
