# Re-derive donor-first system stage coefficients in GSE213621 with and without the
# 68 controls coded as fibrosis_stage 0; read-only inputs.
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, norm
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/histology_anchored_continuum/'
S=B+'molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables/'
sc=pd.read_csv(S+'participant_system_scores.tsv.gz',sep='\t')
sc=sc[sc.aggregation_method=='collection_balanced']
man=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv',sep='\t',keep_default_na=False)
print(pd.crosstab(man[man.dataset=='GSE213621'].group_binary, man[man.dataset=='GSE213621'].fibrosis_stage))
sc=sc.merge(man[['sample_id','group_binary']],on='sample_id',how='left')
ref=pd.read_csv(S+'system_fibrosis_stage_cohort_models.tsv',sep='\t')
ref=ref[ref.aggregation_method=='collection_balanced']
out=[]
for ds in ['GSE162694','GSE213621']:
    for label,sub in [('as_published',None),('disease_only','Disease')]:
        d=sc[(sc.dataset==ds)&sc.fibrosis_stage.notna()].copy()
        if sub: d=d[d.group_binary==sub]
        for cid,g in d.groupby('community_id'):
            g=g.copy()
            g['fibrosis_stage_z']=(g.fibrosis_stage-g.fibrosis_stage.mean())/g.fibrosis_stage.std()
            m=smf.ols('system_score_z ~ fibrosis_stage_z + inferred_sex',data=g).fit()
            out.append(dict(dataset=ds,roster=label,community_id=cid,n=len(g),beta=m.params['fibrosis_stage_z'],se=m.bse['fibrosis_stage_z']))
o=pd.DataFrame(out)
chk=o[o.roster=='as_published'].merge(ref[['dataset','community_id','n_participants','beta']],on=['dataset','community_id'])
print('reproduction max abs diff beta:',(chk.beta_x-chk.beta_y).abs().max(),'n match:',(chk.n==chk.n_participants).all(), chk.n.unique())
# fixed-effect meta per roster (GSE162694 identical in both)
res={}
for label in ['as_published','disease_only']:
    a=o[(o.dataset=='GSE162694')&(o.roster=='as_published')].set_index('community_id')
    b=o[(o.dataset=='GSE213621')&(o.roster==label)].set_index('community_id')
    w1=1/a.se**2; w2=1/b.se**2
    beta=(a.beta*w1+b.beta*w2)/(w1+w2); se=np.sqrt(1/(w1+w2))
    p=2*norm.sf(np.abs(beta/se))
    from statsmodels.stats.multitest import multipletests
    q=multipletests(p,method='fdr_bh')[1]
    conc=np.sign(a.beta)==np.sign(b.beta)
    res[label]=pd.DataFrame(dict(beta=beta,q=q,conc=conc,b213=b.beta))
    print(label,'GSE213621 n=',b.n.iloc[0],' stage-supported (BH q<0.05 & concordant):',((q<0.05)&conc).sum())
r=res['as_published'].join(res['disease_only'],lsuffix='_pub',rsuffix='_dz')
print('Spearman meta stage beta published vs disease-only:',spearmanr(r.beta_pub,r.beta_dz))
print('GSE213621 beta ratio median disease-only/published:',np.median(r.b213_dz/r.b213_pub))
print(r.describe().to_string())
# compare to continuum meta beta
cm=pd.read_csv(S+'system_continuum_meta_analysis.tsv',sep='\t')
cm=cm[(cm.aggregation_method=='collection_balanced')&(cm.axis_id=='fixed_projection')].set_index('community_id')
for label in ['pub','dz']:
    x=r['beta_'+label]; y=cm.loc[r.index,'fixed_beta']
    print(label,'donor-first stage vs continuum: rho',round(spearmanr(x,y)[0],4),'same direction',int((np.sign(x)==np.sign(y)).sum()))
o.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R05/system_stage_betas_with_without_controls.tsv',sep='\t',index=False)
