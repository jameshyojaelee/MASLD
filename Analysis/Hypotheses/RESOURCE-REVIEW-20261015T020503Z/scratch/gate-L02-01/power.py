import pandas as pd, numpy as np, gzip
from scipy import stats
np.random.seed(20261015)
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
L=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01'
deg=pd.read_csv(B+'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
deg['ens']=deg.gene.str.split('.').str[0]
can=deg[(deg.padj<0.05)&(deg.logFC.abs()>0.5)]
print('deg rows',len(deg),'canonical',len(can),'up',(can.logFC>0).sum(),'down',(can.logFC<0).sum())
g=pd.read_csv(B+'data/external/allelic_refs/gtex_v8/gene_tpm_2017-06-05_v8_liver.gct.gz',sep='\t',skiprows=2)
print('gct shape',g.shape, list(g.columns[:3]))
g['ens']=g.Name.str.split('.').str[0]
print('gct dup ens',g.ens.duplicated().sum())
g=g.drop_duplicates('ens').set_index('ens')
samp=[c for c in g.columns if c.startswith('GTEX-')]
print('gct samples',len(samp))
inter=can[can.ens.isin(g.index)]
print('canonical in gct',len(inter),'up',(inter.logFC>0).sum(),'down',(inter.logFC<0).sum())
meta=pd.read_csv(f'{L}/liver_rna_meta.tsv',sep='\t')
gr=pd.read_csv(f'{L}/note_grades.tsv',sep='\t')
meta=meta.merge(gr[['donor','grade','ord','pct','fib_or_nash']],on='donor')
meta=meta[meta.DTHHRDY.notna()].copy()
col={('-'.join(c.split('-')[:2])):c for c in samp}
meta['col']=meta.donor.map(col)
meta=meta[meta.col.notna()]
print('model n',len(meta))
E=np.log2(g.loc[inter.ens.values,meta.col.values].astype(float)+1)
medtpm=g.loc[inter.ens.values,meta.col.values].median(axis=1)
print('canonical in gct: median TPM>=1:',(medtpm>=1).sum(),' >=0.1:',(medtpm>=0.1).sum(),' ==0:',(medtpm==0).sum())
# design
X=pd.get_dummies(meta[['AGE']].astype(str),drop_first=True).astype(float)
X=pd.concat([X,pd.get_dummies(meta.DTHHRDY.astype(int).astype(str),prefix='H',drop_first=True).astype(float)],axis=1)
X['sex']=(meta.SEX==2).astype(float).values
X['isch']=meta.SMTSISCH.values; X['rin']=meta.SMRIN.values
X.insert(0,'int',1.0)
Xm=X.values; n,p=Xm.shape; print('design n',n,'columns incl intercept',p)
Q,_=np.linalg.qr(Xm)
Y=E.values.T  # n x genes
Yr=Y-Q@(Q.T@Y)
dfres=n-p-1
sdY=np.sqrt((Yr**2).sum(0)/dfres)
# exposure for null: ordinal with no_mention->0, unspecified->2
x=meta.ord.copy(); x[meta.grade=='no_mention']=0; x[meta.grade=='unspecified']=2
x=x.values.astype(float)
strata=meta.DTHHRDY.astype(int).values
lfc=inter.logFC.values
def stat(xv):
    xr=xv-Q@(Q.T@xv)
    b=xr@Yr/(xr@xr)
    res=Yr-np.outer(xr,b)
    se=np.sqrt((res**2).sum(0)/dfres/(xr@xr))
    return b,b/se
nul_b=[];nul_t=[]
for i in range(1000):
    xp=x.copy()
    for s in np.unique(strata):
        idx=np.where(strata==s)[0]; xp[idx]=np.random.permutation(x[idx])
    b,t=stat(xp)
    nul_b.append(stats.spearmanr(lfc,b)[0]); nul_t.append(stats.spearmanr(lfc,t)[0])
nul_b=np.array(nul_b); nul_t=np.array(nul_t)
print('NULL cross-gene spearman (slope): sd %.3f  2.5/97.5%% %.3f %.3f'%(nul_b.std(),*np.quantile(nul_b,[.025,.975])))
print('NULL cross-gene spearman (t):     sd %.3f  2.5/97.5%% %.3f %.3f'%(nul_t.std(),*np.quantile(nul_t,[.025,.975])))
print('MDE approx (1.96+0.84)*sd: slope %.3f t %.3f'%(2.8*nul_b.std(),2.8*nul_t.std()))
# null if genes were independent: sd ~ 1/sqrt(G-1)
print('independent-gene sd would be %.3f'%(1/np.sqrt(len(lfc)-1)))
# PLIN2 power
ip=np.where(inter.symbol.values=='PLIN2')[0]
print('PLIN2 in canonical∩gct:',len(ip)>0, inter[inter.symbol=='PLIN2'][['logFC','padj']].to_dict('records'))
pl=np.log2(g.loc['ENSG00000147872',meta.col.values].astype(float)+1).values
plr=pl-Q@(Q.T@pl); sdpl=np.sqrt((plr**2).sum()/dfres)
for lab,xv,mask in [('ordinal no_mention=0, all',x,np.ones(n,bool)),('graded only',meta.ord.values,meta.ord.notna().values)]:
    Qm,_=np.linalg.qr(Xm[mask]); xx=xv[mask].astype(float); xr=xx-Qm@(Qm.T@xx)
    yy=pl[mask]; yr=yy-Qm@(Qm.T@yy); df=mask.sum()-p-1; sdy=np.sqrt((yr**2).sum()/df)
    for beta in [0.168,0.168*0.5]:
        se=sdy/np.sqrt(xr@xr); z=beta/se
        pw=stats.norm.sf(1.96-z)+stats.norm.cdf(-1.96-z)
        print(f'PLIN2 {lab}: n {mask.sum()} resid SD log2TPM {sdy:.3f} SD(resid x) {np.sqrt(xr@xr/df):.3f} beta {beta:.3f} -> z {z:.2f} power(two-sided .05) {pw:.2f}; MDE beta {2.8*se:.3f}')
print('PLIN2 median TPM',np.median(2**pl-1).round(1))
# per-gene MDE partial r
for a in [0.05/len(lfc),0.005,0.05]:
    z=stats.norm.isf(a/2)+0.84
    print('per-gene alpha %.2e MDE partial r ~ %.3f'%(a,np.tanh(z/np.sqrt(n-p-3))))
# plausible effect: H8 steatosis meta vs bulk logFC
h=pd.read_csv(B+'Analysis/Hypotheses/H8-20261009T015529Z/run-20261009T020920Z/meta_steatosis.tsv.gz',sep='\t')
hh=inter.merge(h[['gene_id','mu','q_meta']],left_on='ens',right_on='gene_id')
print('H8 steatosis meta ∩ canonical∩gct:',len(hh),' spearman(bulk logFC, biopsy steatosis-grade logFC) = %.3f p %.2g'%stats.spearmanr(hh.logFC,hh.mu))
print('  of these, H8 steatosis q<0.05:',(hh.q_meta<0.05).sum(), ' same sign as disease among q<0.05:',(np.sign(hh.logFC)==np.sign(hh.mu))[hh.q_meta<0.05].sum())
print('  median |mu| per grade among canonical: %.3f'%hh.mu.abs().median())
# expression residual SD of canonical genes
print('median resid SD (log2 TPM+1) across canonical genes: %.3f'%np.median(sdY))
np.save(f'{L}/null_spearman_slope.npy',nul_b)
