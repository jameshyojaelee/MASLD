"""Independent recheck of the Fig 5D locus-matched expectation. Peak membership by
count of peaks started minus peaks ended (no interval merging). Seed 20260923."""
import sys, numpy as np, pandas as pd
from scipy import stats
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/'
OUT=sys.argv[1]
prim=pd.read_csv(B+'atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv',sep='\t')
prim['key']=prim.gwas_name+'|'+prim.ensembl+'|'+prim.signal_pair_index.astype(str)
keys=set(prim.key)
parts=[]
for ch in pd.read_csv(B+'atac-context-v3-candidate-2026-08-11-r1/genetics/context/variant_liftover_audit.tsv',sep='\t',
        usecols=['gwas_name','ensembl','signal_pair_index','snp_pp_h4','mapping_status','hg38_chrom','hg38_position_1based'],
        dtype=str,chunksize=3_000_000):
    k=ch.gwas_name+'|'+ch.ensembl+'|'+ch.signal_pair_index
    ch=ch[k.isin(keys)].copy(); ch['key']=k[k.isin(keys)]
    parts.append(ch)
v=pd.concat(parts,ignore_index=True); del parts
v['pp']=v.snp_pp_h4.astype(float)
v['mapped']=v.mapping_status=='mapped'
v['pos']=pd.to_numeric(v.hg38_position_1based,errors='coerce')
print('variants',len(v),'pairs',v.key.nunique(),flush=True)
pk=pd.read_csv(B+'atac-context-v3-candidate-2026-08-11-r1/consensus_peak_manifest.tsv',sep='\t',
        usecols=['lineage','chrom','start0','end','supported_gse244832','supported_gse281367'])
for c in ['supported_gse244832','supported_gse281367']: pk[c]=pk[c].astype(str).str.upper()=='TRUE'
def member(sub):
    res=np.zeros(len(v),bool)
    for chrom,g in sub.groupby('chrom'):
        s=np.sort(g.start0.values); e=np.sort(g.end.values)
        sel=(v.hg38_chrom==chrom).values & v.mapped.values
        p=v.pos.values[sel]
        res[sel]=(np.searchsorted(s,p,'left')-np.searchsorted(e,p,'left'))>0
    return res
rows=[]
for lin in ['hepatocyte','stellate','macrophage','cholangiocyte','t_nk']:
    pl=pk[pk.lineage==lin]
    a=member(pl[pl.supported_gse244832]); b=member(pl[pl.supported_gse281367])
    d=pd.DataFrame({'key':v.key,'pp':v.pp,'m':v.mapped,'a':a,'b':b}); d['any']=d.a|d.b
    d['ppa']=d.pp*d.a; d['ppb']=d.pp*d.b; d['ppany']=d.pp*d['any']
    g=d.groupby('key').agg(nm=('m','sum'),nany=('any','sum'),oa=('ppa','sum'),ob=('ppb','sum'),oany=('ppany','sum')).reset_index()
    g['lineage']=lin; rows.append(g); print(lin,flush=True)
r=pd.concat(rows); r['exp']=r.nany/r.nm
r=prim.merge(r,on=['key','lineage'],how='left',validate='1:1')
for o,t in [('oa','gse244832_accessible_mass'),('ob','gse281367_accessible_mass'),('oany','any_accessible_mass')]:
    print('max abs diff vs release',o,(r[o]-r[t]).abs().max())
cand=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv',sep='\t')
cand['key']=cand.gwas_name+'|'+cand.ensembl+'|'+cand.signal_pair_index.astype(str)
cm=r.merge(cand[['key','lineage','exp_any','obs_any']],on=['key','lineage'])
print('rows matched',len(cm),'max |exp diff| vs candidate',(cm.exp-cm.exp_any).abs().max(),'max |obs diff|',(cm.oany-cm.obs_any).abs().max())
r['diff']=r.oany-r.exp
gb=r.groupby(['trait_class','lineage','ensembl']).agg(obs=('oany','mean'),exp=('exp','mean'),diff=('diff','mean')).reset_index()
rng=np.random.default_rng(20260923); res=[]
for (tc,lin),x in gb.groupby(['trait_class','lineage']):
    idx=rng.integers(0,len(x),(2000,len(x)))
    rat=x.obs.values[idx].mean(1)/x.exp.values[idx].mean(1)
    res.append(dict(trait_class=tc,lineage=lin,n_genes=len(x),obs=x.obs.mean(),exp=x.exp.mean(),ratio=x.obs.mean()/x.exp.mean(),
        lo=np.percentile(rat,2.5),hi=np.percentile(rat,97.5),n_obs_gt_exp=int((x['diff']>0).sum()),
        wilcoxon_p=stats.wilcoxon(x['diff']).pvalue))
res=pd.DataFrame(res); res['bh10']=stats.false_discovery_control(res.wilcoxon_p)
pd.set_option('display.width',250); print(res.to_string())
res.to_csv(OUT+'/r11_02_background_indep_summary.tsv',sep='\t',index=False)
gb.to_csv(OUT+'/r11_02_gene_obs_exp_indep.tsv',sep='\t',index=False)
