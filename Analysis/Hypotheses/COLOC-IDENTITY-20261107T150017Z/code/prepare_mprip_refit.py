#!/usr/bin/env python3
"""Prepare full-precision LD for every original eQTL-fitted allele row."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

RUN = Path(__file__).resolve().parents[1]
ROOT = RUN.parents[2]
OUT = RUN / 'mprip_refit_inputs'
OUT.mkdir(exist_ok=False)
eq = pd.read_csv(RUN/'mprip_sources/MPRIP_original_eqtl.tsv', sep='\t')
eq = eq[eq.POS.notna() & eq.Beta.notna() & eq.SE.gt(0)]
prefix = ROOT/'GWAS/finemapping/data/ld_ref/polyfun_eur/raw/chr17_15000001_18000001'
meta = pd.read_csv(str(prefix)+'.gz', sep=r'\s+')
meta['SNP'] = meta.chromosome.astype(str)+':'+meta.position.astype(str)+':'+meta.allele1+':'+meta.allele2
assert meta.SNP.is_unique
direct = dict(zip(eq.CHR.astype(str)+':'+eq.POS.astype(str)+':'+eq.EA+':'+eq.NEA, zip(eq.Beta, eq.SE)))
reverse = dict(zip(eq.CHR.astype(str)+':'+eq.POS.astype(str)+':'+eq.NEA+':'+eq.EA, zip(-eq.Beta, eq.SE)))
lookup = {**reverse, **direct}
ix = np.flatnonzero(meta.SNP.isin(lookup))
aligned = meta.iloc[ix].copy()
aligned['beta'] = [lookup[k][0] for k in aligned.SNP]
aligned['standard_error'] = [lookup[k][1] for k in aligned.SNP]
aligned.to_csv(OUT/'eqtl_aligned.tsv', sep='\t', index=False)
raw_lower = sparse.load_npz(str(prefix)+'.npz').tocsr()
small = raw_lower[ix,:][:,ix].toarray().astype(np.float64)
matrix = small + small.T
assert np.allclose(np.diag(matrix), 1, atol=1e-6)
assert np.isfinite(matrix).all()
matrix.ravel(order='F').tofile(OUT/'eqtl_raw_LD_float64.bin')
print('Prepared full eQTL matrix', matrix.shape, flush=True)
