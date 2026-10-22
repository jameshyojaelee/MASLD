# Donor-level system stage models (system_score_z ~ stage_z + sex per cohort, fixed-effect meta, BH over 43) with and without GSE213621 controls.
import numpy as np, pandas as pd
from scipy.stats import spearmanr, norm
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=ROOT+'RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables/'
man=pd.read_csv(ROOT+'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv',sep='\t')
sc=pd.read_csv(S+'participant_system_scores.tsv.gz',sep='\t')
sc=sc[(sc.aggregation_method=='collection_balanced')&sc.dataset.isin(['GSE162694','GSE213621'])].merge(man[['sample_id','group_binary']],on='sample_id')
ref=pd.read_csv(S+'system_fibrosis_stage_meta_analysis.tsv',sep='\t'); ref=ref[ref.aggregation_method=='collection_balanced'].set_index('community_id')
cm=pd.read_csv(S+'system_continuum_meta_analysis.tsv',sep='\t'); cm=cm[(cm.aggregation_method=='collection_balanced')&(cm.axis_id=='fixed_projection')].set_index('community_id')
def bh(p):
    p=np.asarray(p); n=len(p); o=np.argsort(p); q=np.empty(n); q[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; return np.minimum(q,1)
def fit(d):
    st=d.fibrosis_stage.values; st=(st-st.mean())/st.std(ddof=1)
    X=np.column_stack([np.ones(len(d)),st,(d.inferred_sex.values=='M').astype(float)]); y=d.system_score_z.values
    XtXi=np.linalg.inv(X.T@X); b=XtXi@X.T@y; r=y-X@b; se=np.sqrt((r@r)/(len(y)-3)*XtXi[1,1]); return b[1],se,len(y)
for roster in ['published','disease_only']:
    rows=[]
    for cid,g in sc.groupby('community_id'):
        out={}
        for ds in ['GSE162694','GSE213621']:
            d=g[(g.dataset==ds)&g.fibrosis_stage.notna()&g.inferred_sex.notna()]
            if roster=='disease_only' and ds=='GSE213621': d=d[d.group_binary=='Disease']
            out[ds]=fit(d)
        (b1,s1,n1),(b2,s2,n2)=out['GSE162694'],out['GSE213621']
        w1,w2=1/s1**2,1/s2**2; b=(w1*b1+w2*b2)/(w1+w2); se=np.sqrt(1/(w1+w2))
        rows.append(dict(community_id=cid,beta=b,p=2*norm.sf(abs(b/se)),conc=np.sign(b1)==np.sign(b2),n1=n1,n2=n2))
    t=pd.DataFrame(rows).set_index('community_id'); t['q']=bh(t.p.values)
    if roster=='published': print('max|diff| vs stored fixed_beta',np.abs(t.beta-ref.loc[t.index,'fixed_beta']).max())
    x=t.beta; y=cm.loc[t.index,'fixed_beta']
    print(roster,'n',t.n1.iloc[0],t.n2.iloc[0],'stage-supported (BH43 q<0.05 & concordant)',int(((t.q<0.05)&t.conc).sum()),
          'rho vs continuum',round(spearmanr(x,y)[0],4),'same dir',int((np.sign(x)==np.sign(y)).sum()))
