import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(20261015)
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
L=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01'
deg=pd.read_csv(B+'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
deg['ens']=deg.gene.str.split('.').str[0]
can=deg[(deg.padj<0.05)&(deg.logFC.abs()>0.5)]
g=pd.read_csv(B+'data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz',sep='\t',skiprows=2)
g['ens']=g.Name.str.split('.').str[0]; g=g.drop_duplicates('ens').set_index('ens')
samp=[c for c in g.columns if c.startswith('GTEX-')]
inter=can[can.ens.isin(g.index)].copy()
meta=pd.read_csv(f'{L}/liver_rna_meta.tsv',sep='\t').merge(pd.read_csv(f'{L}/note_grades.tsv',sep='\t')[['donor','grade','ord']],on='donor')
meta=meta[meta.DTHHRDY.notna()].copy()
col={('-'.join(c.split('-')[:2])):c for c in samp}; meta['col']=meta.donor.map(col); meta=meta[meta.col.notna()]
T=g.loc[inter.ens.values,meta.col.values].astype(float)
Y=np.log2(T.values.T+1)
X=pd.get_dummies(meta[['AGE']].astype(str),drop_first=True).astype(float)
X=pd.concat([X,pd.get_dummies(meta.DTHHRDY.astype(int).astype(str),prefix='H',drop_first=True).astype(float)],axis=1)
X['sex']=(meta.SEX==2).astype(float).values; X['isch']=meta.SMTSISCH.values; X['rin']=meta.SMRIN.values; X.insert(0,'int',1.0)
Q,_=np.linalg.qr(X.values); n,p=X.shape
x=meta.ord.copy(); x[meta.grade=='no_mention']=0; x[meta.grade=='unspecified']=2; x=x.values.astype(float)
strata=meta.DTHHRDY.astype(int).values
h=pd.read_csv(B+'Analysis/Hypotheses/H8-20261009T015529Z/run-20261009T020920Z/meta_steatosis.tsv.gz',sep='\t').set_index('gene_id')
beta_full=inter.ens.map(h.mu).values  # NaN where H8 did not test
lfc=inter.logFC.values
var_ok=Y.std(0)>0
def perm():
    xp=x.copy()
    for s in np.unique(strata):
        i=np.where(strata==s)[0]; xp[i]=rng.permutation(x[i])
    return xp
def cross(Ys,xv,mask):
    Yr=Ys-Q@(Q.T@Ys); xr=xv-Q@(Q.T@xv); b=xr@Yr/(xr@xr)
    return stats.spearmanr(lfc[mask],b[mask])[0]
for label,mask,bvec in [('all 1,587 (beta 0 where H8 untested)',var_ok,np.nan_to_num(beta_full)),
                        ('723 with H8 estimate',var_ok&~np.isnan(beta_full),np.nan_to_num(beta_full)),
                        ('median TPM>=1 genes',(np.median(T.values,1)>=1),np.nan_to_num(beta_full))]:
    null=[];alt={1.0:[],0.5:[],0.25:[]}
    for i in range(400):
        xp=perm(); null.append(cross(Y,xp,mask))
        for a in alt:
            alt[a].append(cross(Y+np.outer(xp,a*bvec),xp,mask))
    null=np.array(null); lo,hi=np.quantile(null,[0.025,0.975])
    out=f'{label}: genes {mask.sum()}; null sd {null.std():.3f} [{lo:.3f},{hi:.3f}]'
    for a in alt:
        v=np.array(alt[a]); out+=f'; atten {a}: mean rho {v.mean():.3f} power {np.mean((v>hi)|(v<lo)):.2f}'
    print(out)
