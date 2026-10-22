# Blinded power check: null spread of Spearman(TWAS z, bulk DE t) under circular genomic shifts.
# The observed (shift 0) statistic and shifts within 500 genes of 0 are NOT computed.
import pandas as pd, numpy as np, glob
from scipy.stats import rankdata
B='RNA-seq/results/causal_inference/otters_broadaway'
d=pd.read_csv('RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
d['ens']=d.gene.str.split('.').str[0]
rng=np.random.default_rng(20261015)
for t in ['ukbb_alt','ghodsian']:
    parts=[]
    for f in glob.glob(f'{B}/{t}/chr*/P0.05.txt'):
        parts.append(pd.read_csv(f,sep='\t'))
    o=pd.concat(parts).drop_duplicates('TargetID')
    m=o.merge(d[['ens','t']],left_on='TargetID',right_on='ens').dropna(subset=['FUSION_Z','t'])
    m=m.sort_values(['CHROM','GeneStart']).reset_index(drop=True)
    n=len(m); zr=rankdata(m.FUSION_Z.values); tr=rankdata(m.t.values)
    zr=(zr-zr.mean())/zr.std(); tr=(tr-tr.mean())/tr.std()
    shifts=rng.integers(500,n-500,size=2000)
    null=np.array([np.mean(np.roll(zr,s)*tr) for s in shifts])
    print(t,'genes',n,'chromosomes',m.CHROM.nunique(),'null mean %.4f null SD %.4f  null 2.5/97.5%% %.4f %.4f  independence SD %.4f'%(null.mean(),null.std(),*np.percentile(null,[2.5,97.5]),1/np.sqrt(n)))
