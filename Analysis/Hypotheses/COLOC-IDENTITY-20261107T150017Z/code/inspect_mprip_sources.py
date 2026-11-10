#!/usr/bin/env python3
"""Recover source REF/ALT and full-precision LD without choosing by association."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse

RUN=Path(__file__).resolve().parents[1];ROOT=RUN.parents[2]
OUT=RUN/'mprip_sources';OUT.mkdir(exist_ok=False)
raw=ROOT/'GWAS/MR_Data/MVP/raw/labs/MVP_R4.1000G_AGR.ALT_Mean_INT.EUR.GIA.dbGaP.txt.gz'
hits=[];n=0
for d in pd.read_csv(raw,sep='\t',chunksize=500000):
    hits.append(d[(d.SNP_ID=='rs56073511') | ((d.chrom.astype(str)=='17') & d.pos.between(17261520,17261535))])
    n+=len(d)
    if n%5000000==0:print('MVP source records read',n,flush=True)
pd.concat(hits).to_csv(OUT/'MVP_original_repeat_records.tsv',sep='\t',index=False)
eqtl=ROOT/'data/broadaway_eqtl/chr17_marginal_summary_results.tsv'
gene=[];repeat=[]
for d in pd.read_csv(eqtl,sep='\t',chunksize=500000):
    gene.append(d[d.ENSG=='ENSG00000133030']);repeat.append(d[d.POS==17164840])
pd.concat(gene).to_csv(OUT/'MPRIP_original_eqtl.tsv',sep='\t',index=False)
pd.concat(repeat).to_csv(OUT/'Broadaway_repeat_all_genes.tsv',sep='\t',index=False)
prior=ROOT/'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/pilot_22811698/exports'
aligned=pd.read_csv(prior/'MVP_ALT_EUR__ENSG00000133030__aligned_variants.tsv',sep='\t')
prefix=ROOT/'GWAS/finemapping/data/ld_ref/polyfun_eur/raw/chr17_15000001_18000001'
meta=pd.read_csv(str(prefix)+'.gz',sep=r'\s+')
meta['SNP']=meta.chromosome.astype(str)+':'+meta.position.astype(str)+':'+meta.allele1+':'+meta.allele2
assert meta.SNP.is_unique
lookup=pd.Series(np.arange(len(meta)),index=meta.SNP)
assert aligned.SNP.isin(lookup.index).all()
ix=lookup.loc[aligned.SNP].to_numpy()
raw_lower=sparse.load_npz(str(prefix)+'.npz').tocsr()
assert raw_lower.shape==(len(meta),len(meta))
small=raw_lower[ix,:][:,ix].toarray().astype(np.float64)
matrix=small+small.T
assert np.allclose(np.diag(matrix),1,atol=1e-6)
assert np.isfinite(matrix).all() and np.allclose(matrix,matrix.T,atol=1e-12)
matrix.ravel(order='F').tofile(OUT/'MPRIP_raw_single_window_LD_float64.bin')
aligned.to_csv(OUT/'MPRIP_raw_LD_order.tsv',sep='\t',index=False)
meta[meta.position==17164840].to_csv(OUT/'PolyFun_original_repeat_metadata.tsv',sep='\t',index=False)
rows=np.flatnonzero(aligned.position.eq(17164840))
pd.DataFrame(matrix[np.ix_(rows,rows)],index=aligned.SNP.iloc[rows],columns=aligned.SNP.iloc[rows]).to_csv(OUT/'PolyFun_repeat_original_correlations.tsv',sep='\t')
summary=dict(MVP_source_records=n,aligned_variants=len(aligned),reference_rows=len(meta),
 raw_source=str(prefix),matrix_shape=list(matrix.shape),matrix_dtype='float64_column_major',
 correction_role='Original single-window, unrounded LD; no row selection or fitted result in this step',
 reference_source_matched_to_MVP=False)
(OUT/'source_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
