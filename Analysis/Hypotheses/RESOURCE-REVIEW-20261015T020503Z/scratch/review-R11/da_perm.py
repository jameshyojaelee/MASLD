# Donor-label permutation null for the ATAC v3 cross-cohort directional replication statistic.
# Statistic: number of GSE281367 BH-significant peaks (q<0.05) whose GSE244832 MASH-minus-NORMAL
# logFC has the same sign. Reported: 849/1016 hepatocyte, 15/15 stellate, binomial q 2.64e-110 / 4.58e-5.
import numpy as np, gzip, csv, itertools
from scipy.io import mmread
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1'
# saved joint results
res={}
with gzip.open(f'{C}/da/da_peak_results.tsv.gz','rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        res.setdefault(r['lineage'],[]).append(r)
qc={}
for r in csv.DictReader(open(f'{C}/da/da_model_qc.tsv'),delimiter='\t'):
    if r['cohort']=='GSE244832': qc[(r['lineage'],r['donor_id'])]=float(r['effective_library_size'])
for L in ['hepatocyte','stellate']:
    rows=res[L]
    disc=[r for r in rows if float(r['qvalue_gse281367'])<0.05]
    sign281=np.array([np.sign(float(r['logFC_gse281367'])) for r in disc])
    saved244=np.array([np.sign(float(r['logFC_gse244832'])) for r in disc])
    print(L,'discoveries',len(disc),'saved same-sign',int((saved244==sign281).sum()))
    donors=list(csv.DictReader(open(f'{C}/counts/GSE244832/{L}.donors.tsv'),delimiter='\t'))
    peaks=[r['peak_coordinate'] for r in csv.DictReader(open(f'{C}/counts/GSE244832/{L}.peaks.tsv'),delimiter='\t')]
    M=mmread(gzip.open(f'{C}/counts/GSE244832/{L}.mtx.gz')).tocsc()
    pidx={p:i for i,p in enumerate(peaks)}
    keep=[i for i,d in enumerate(donors) if d['contrast_eligible']=='TRUE' and d['condition'] in ('NORMAL','MASH')]
    ids=[donors[i]['donor_id'] for i in keep]; cond=np.array([donors[i]['condition']=='MASH' for i in keep])
    cols=[pidx[r['peak_coordinate']] for r in disc]
    X=M[keep][:,cols].toarray().astype(float)
    lib=np.array([qc[(L,i)] for i in ids])
    logcpm=np.log2((X+0.5)/(lib[:,None]+1.0)*1e6)
    def stat(lab):
        d=logcpm[lab].mean(0)-logcpm[~lab].mean(0)
        return int((np.sign(d)==sign281).sum())
    obs=stat(cond)
    n=len(ids); k=int(cond.sum())
    null=[]
    for comb in itertools.combinations(range(n),k):
        lab=np.zeros(n,bool); lab[list(comb)]=True; null.append(stat(lab))
    null=np.array(null)
    print(f'  donors n={n} MASH={k} NORMAL={n-k}; observed same-sign (unweighted log-CPM mean difference) = {obs}/{len(disc)}')
    print(f'  exact donor-label permutation: {len(null)} relabelings; null median {np.median(null):.0f}, p95 {np.percentile(null,95):.0f}, max {null.max()}; draws >= observed {int((null>=obs).sum())}; one-sided p = {(null>=obs).mean():.4f}; floor 1/{len(null)} = {1/len(null):.2e}')
