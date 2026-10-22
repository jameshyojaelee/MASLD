# H72 power from SEs and the sex-combined (both_sexes) effect only. Reads se (col 8) of the
# female and male files and beta/se of the both_sexes file. Never reads female or male beta.
import gzip, math, statistics as st
from multiprocessing import Pool
from scipy.stats import norm
import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
W=R+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72'
L=R+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L06'
S=R+'/GWAS/MR_Data/PanUKBB/sex_stratified'
PC={'ALT':'30620','AST':'30650','GGT':'30730'}
best=pd.read_csv(W+'/best_enzyme_row.tsv',sep='\t')
best=best[best.pp>=0.5]
def scan(t):
    want={}
    for r in best[best.trait==t].itertuples():
        c,p=r.top_snp.split(':')[:2]; want[(c,int(p))]=('gene',r.ensembl,r.pp)
    for l in open(f'{L}/ukbb_{t.lower()}_clumped500kb.tsv'):
        a=l.split('\t'); want.setdefault((a[0],int(a[1])),('lead',a[0]+':'+a[1],None))
        if want[(a[0],int(a[1]))][0]=='gene': want[('L',a[0],int(a[1]))]=('lead',a[0]+':'+a[1],None)
    al={}
    for k in want:
        if k[0]!='L': al[k]=None
    # alleles: COLOC source for genes, clumped file for leads
    src={}
    with open(f'{R}/GWAS/finemapping/data/sumstats/UKBB_{t}_reformatted_hg19.tsv') as f:
        next(f)
        for l in f:
            a=l.split('\t',4); k=(a[0],int(a[1]))
            if k in al: src[k]=frozenset((a[2],a[3]))
    for l in open(f'{L}/ukbb_{t.lower()}_clumped500kb.tsv'):
        a=l.split('\t'); src.setdefault((a[0],int(a[1])),frozenset((a[2],a[3])))
    fs=[gzip.open(f'{S}/{PC[t]}_irnt.gwas.imputed_v3.{s}.varorder.tsv.bgz','rt') for s in ('female','male','both_sexes')]
    for f in fs: next(f)
    rec={}
    for lf,lm,lb in zip(*fs):
        a=lf.split('\t',10); b=lm.split('\t',10); c=lb.split('\t',10)
        assert a[0]==b[0]==c[0]
        v=a[0].split(':'); k=(v[0],int(v[1]))
        if k in src and frozenset((v[2],v[3]))==src[k] and a[3]!='true' and b[3]!='true':
            try: rec[k]=(float(a[8]),float(b[8]),float(c[7]),float(c[8]))
            except ValueError: pass
    rows=[]
    for k,(kind,id_,pp) in want.items():
        kk=k[1:] if k[0]=='L' else k
        if kk in rec:
            sF,sM,bB,sB=rec[kk]; rows.append((t,kind,id_,pp,sF,sM,bB,sB))
    return rows
with Pool(3) as p: rows=[r for x in p.map(scan,list(PC)) for r in x]
d=pd.DataFrame(rows,columns=['trait','kind','id','pp','sF','sM','bB','sB'])
d['sD']=(d.sF**2+d.sM**2)**0.5
def p_equiv(r,lo=2/3,hi=1.5,z=1.959964):
    # P(95% Fieller CI of bF/bM inside [lo,hi]) if bF=bM=bB, normal approx on log scale
    b=abs(r.bB); se=math.sqrt((r.sF/b)**2+(r.sM/b)**2)
    if abs(r.bB)/r.sM<z: return 0.0
    # CI inside bound iff |log ratio_hat| + z*se <= log(1.5); log ratio_hat ~ N(0,se)
    h=math.log(hi)-z*se
    return max(0.0,2*norm.cdf(h/se)-1) if h>0 else 0.0
d['p_equiv']=d.apply(p_equiv,axis=1)
d.to_csv(W+'/mde_rows.tsv',sep='\t',index=False)
out=open(W+'/mde.out','w')
for kind,fam in [('gene',None),('lead',None)]:
    x=d[d.kind==kind]; m=len(x)
    zc=norm.isf(0.05/m/2); 
    print(f'{kind}: n {m}; by trait {dict(x.trait.value_counts())}; median SE(diff) {x.sD.median():.5f} SD; BH single-hit z {zc:.3f}; MDE80 {(zc+norm.isf(0.2))*x.sD.median():.4f} SD',file=out)
    print(f'{kind}: median |beta_both| {x.bB.abs().median():.4f} SD; n with |beta_both| >= MDE80 {(x.bB.abs()>=(zc+norm.isf(0.2))*x.sD).sum()}',file=out)
    print(f'{kind}: expected homogeneous-within-[0.67,1.5] calls if every effect were equal: {x.p_equiv.sum():.1f} of {m}; n with P>=0.8 {(x.p_equiv>=0.8).sum()}',file=out)
    if kind=='gene':
        y=x[x.pp>=0.8]; print(f'gene PP>=0.8: n {len(y)}; expected homogeneous {y.p_equiv.sum():.1f}; n with P>=0.8 {(y.p_equiv>=0.8).sum()}',file=out)
out.close()
