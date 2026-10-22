# Orientation check between TWAS traits (no bulk DE involved): NAFLD (Ghodsian) vs UKBB ALT TWAS z.
import pandas as pd, glob
from scipy.stats import spearmanr
B='RNA-seq/results/causal_inference/otters_broadaway'
def load(t):
    return pd.concat([pd.read_csv(f,sep='\t') for f in glob.glob(f'{B}/{t}/chr*/P0.05.txt')]).drop_duplicates('TargetID').set_index('TargetID').FUSION_Z
a=load('ukbb_alt'); n=load('ghodsian'); j=pd.concat([a,n],axis=1,keys=['alt','nafld']).dropna()
print('shared genes',len(j),'spearman all %.3f'%spearmanr(j.alt,j.nafld)[0])
s=j[(j.alt.abs()>4)&(j.nafld.abs()>4)]
print('genes |z|>4 in both',len(s),'same sign',int((s.alt*s.nafld>0).sum()))
