import pandas as pd, numpy as np
from scipy.stats import fisher_exact, hypergeom
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/'
e=pd.read_parquet(D+'passport_evidence_long.parquet')
w=e.pivot_table(index='ensembl_id',columns='evidence_domain',values='call_state',aggfunc='first')
J=w[w.genetics.isin(['supported','indeterminate'])&w.transcriptomics.isin(['supported','indeterminate'])].copy()
J['G']=J.genetics=='supported'; J['R']=J.transcriptomics=='supported'
N=len(J); a=int((J.G&J.R).sum()); nG=int(J.G.sum()); nR=int(J.R.sum())
print('Catalog joint N',N,'G',nG,'R',nR,'both',a)
for NN in [14682,14704,14720]:
    E=nG*nR/NN
    tab=[[a,nG-a],[nR-a,NN-nG-nR+a]]
    orr,p=fisher_exact(tab); pl=hypergeom.cdf(a,NN,nG,nR)
    print(NN,'expected %.4f'%E,'OR %.4f'%orr,'two-sided Fisher p %.4f'%p,'lower-tail p %.4f'%pl, 'nonoverlap of G %.4f'%(1-a/nG), 'expected nonoverlap %.4f'%(1-nR/NN))
c=pd.read_csv(B+'coefficient_table.tsv.gz',sep='\t'); c['ens']=c.gene.str.split('.').str[0]
J=J.join(c.set_index('ens')[['AveExpr']])
print('missing AveExpr',J.AveExpr.isna().sum())
J['dec']=pd.qcut(J.AveExpr,10,labels=False)
st=J.groupby('dec').agg(n=('G','size'),g=('G','sum'),r=('R','sum'),b=('G',lambda s: int((s & J.loc[s.index,'R']).sum())))
st['E']=st.g*st.r/st.n; print(st); print('expression-decile stratified expected',round(st.E.sum(),3),'observed',a)
rng=np.random.default_rng(20261015); B_=10000; cnt=np.empty(B_)
grp=[(J.loc[J.dec==d,'G'].values, J.loc[J.dec==d,'R'].values) for d in range(10)]
for i in range(B_):
    cnt[i]=sum((rng.permutation(gv)&rv).sum() for gv,rv in grp)
print('stratified permutation: mean %.3f, lower-tail p (<=obs) %.4f, two-sided %.4f, B=%d seed 20261015'%(cnt.mean(), (1+(cnt<=a).sum())/(B_+1), min(1,2*min((1+(cnt<=a).sum())/(B_+1),(1+(cnt>=a).sum())/(B_+1))), B_))
