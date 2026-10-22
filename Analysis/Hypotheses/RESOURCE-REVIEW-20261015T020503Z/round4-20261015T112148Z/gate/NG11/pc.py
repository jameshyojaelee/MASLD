# positive control: SLC2A9 urate (Neale 30880 irnt female vs male), hg19 chr4:9.0-10.5 Mb
import gzip, math
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/public_r4/neale_30880'
G='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11'
F=gzip.open(f'{D}/30880_irnt.gwas.imputed_v3.female.varorder.tsv.bgz','rt'); M=gzip.open(f'{D}/30880_irnt.gwas.imputed_v3.male.varorder.tsv.bgz','rt')
next(F); next(M); best=None; nF=nM=None
for lf,lm in zip(F,M):
    a=lf.rstrip('\n').split('\t'); b=lm.rstrip('\n').split('\t')
    if nF is None: nF,nM=a[4],b[4]
    c,pos=a[0].split(':')[:2]
    if c!='4': 
        if best and c not in ('4',): break
        continue
    pos=int(pos)
    if not (9000000<=pos<=10500000): continue
    assert a[0]==b[0]
    if a[3]!='false' or b[3]!='false' or 'NaN' in (a[7],b[7]): continue
    bF,sF,bM,sM=map(float,(a[7],a[8],b[7],b[8]))
    zc=(bF/sF**2+bM/sM**2)/math.sqrt(1/sF**2+1/sM**2)  # inverse-variance combined z: lead chosen sex-blind
    if best is None or abs(zc)>abs(best[0]): best=(zc,a[0],bF,sF,bM,sM,a[2])
zc,v,bF,sF,bM,sM,maf=best
zd=(bF-bM)/math.sqrt(sF**2+sM**2)
open(G+'/positive_control.txt','w').write(f'nF {nF} nM {nM}\nlead {v} MAF {maf} combined z {zc:.1f}\nbeta_F {bF} se_F {sF}\nbeta_M {bM} se_M {sM}\nz_diff {zd:.2f}  |beta_F|>|beta_M|: {abs(bF)>abs(bM)}\n')
