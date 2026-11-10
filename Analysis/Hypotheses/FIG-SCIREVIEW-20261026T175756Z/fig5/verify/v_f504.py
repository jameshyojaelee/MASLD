# F5-04: do the 5 displayed rows drive the curated-set NES? Own weighted-KS GSEA (p=1), gene-set permutation null, seed 20261026.
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
de=pd.read_csv(R+"Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv",sep="\t")
de=de[np.isfinite(de.t)].drop_duplicates("gene")
de=de.sort_values("t",ascending=False).reset_index(drop=True)
genes=de.gene.values; t=de.t.values; N=len(t); absr=np.abs(t)
lipid=["PLIN1","PLIN2","PLIN3","PLIN4","PLIN5","ABHD5","CIDEA","CIDEB","CIDEC","G0S2","HILPDA","ACLY","FASN","ACACA","SCD","DGAT1","DGAT2","LPIN1"]
immune=["HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","CD74","CD68","CD163","AIF1","TYROBP","FCER1G","LGALS3","CTSS","LAPTM5","MARCO","VSIG4","CAPG","SPP1","CASP1","PYCARD","S100A8","S100A9","PTPRC","LCP1","CORO1A","ITGAL","C1QA","C1QB","C1QC","C3"]
shown={"lipid":["PLIN2","ABHD5","PLIN4","PLIN1","ACLY"],"immune":["CAPG","CASP1","CORO1A","CD74","HLA-DQA1"]}
idx={g:i for i,g in enumerate(genes)}
def es_from_pos(pos):
    pos=np.sort(pos); k=len(pos)
    w=absr[pos]; NR=w.sum()
    hit=np.cumsum(w)/NR   # running sum at each hit (after)
    miss_before=(pos-np.arange(k))/(N-k)  # misses before each hit
    up=hit-miss_before
    down=(hit-w/NR)-miss_before  # value just before each hit
    mx=up.max(); mn=min(down.min(),0)
    return mx if abs(mx)>=abs(mn) else mn
rng=np.random.default_rng(20261026)
def nes(S,nperm=20000):
    pos=np.array([idx[g] for g in S if g in idx]); k=len(pos)
    es=es_from_pos(pos)
    null=np.array([es_from_pos(rng.choice(N,k,replace=False)) for _ in range(nperm)])
    same=null[null>0] if es>0 else null[null<0]
    n_=es/abs(same.mean())
    p=( (same>=es).sum()+1)/(len(same)+1) if es>0 else ((same<=es).sum()+1)/(len(same)+1)
    return k,es,n_,p
for name,S in [("lipid",lipid),("immune",immune)]:
    k,es,n_,p=nes(S); print(f"{name} full: size {k} ES {es:.3f} NES {n_:.2f} p {p:.2e}")
    S2=[g for g in S if g not in shown[name]]
    k,es,n_,p=nes(S2); print(f"{name} minus 5 displayed: size {k} ES {es:.3f} NES {n_:.2f} p {p:.2e}")
    tt=de.set_index("gene").t
    print("  displayed t:",{g:round(tt[g],2) for g in shown[name]})
    print("  other measured t:",{g:round(tt[g],2) for g in S2 if g in tt.index})
