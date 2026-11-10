from pathlib import Path
import json
import numpy as np
import pandas as pd
run=Path(__file__).resolve().parents[1]; p=run/'gram_inputs'
x=np.load(p/'standardized_dosages.npz')['X']
ids=pd.read_csv(p/'genotype_LD_order.tsv',sep='\t').SNP
r=np.memmap(p/'genotype_LD_float64.bin',dtype='float64',mode='r',shape=(len(ids),len(ids)),order='F')
assert x.shape==(388,3926) and len(ids)==len(set(ids))
assert np.max(abs(x.mean(axis=0)))<1e-14
assert np.max(abs(np.sum(x*x,axis=0)-1))<1e-12
ix=np.unique(np.r_[np.linspace(0,len(ids)-1,39,dtype=int),np.flatnonzero(ids.str.startswith('17:17164840:'))])
assert np.max(abs(r[np.ix_(ix,ix)]-x[:,ix].T@x[:,ix]))<1e-12
assert np.max(abs(np.diag(r)-1))<1e-12
meta=pd.read_csv(p/'mapped_reference_events.tsv',sep='\t')
assert meta.reciprocal_unique.all()
people=pd.read_csv(p/'reference_participants.tsv',sep='\t',dtype=str)
assert people.FamilyID.is_unique and people.SampleID.is_unique
assert people.Superpopulation.eq('EUR').all() and people.FatherID.eq('0').all() and people.MotherID.eq('0').all()
coverage=pd.read_csv(p/'variant_coverage.tsv',sep='\t')
repeat=coverage[coverage.position.eq(17164840)]
assert len(repeat)==2 and repeat.status.eq('retained').all()
repeat.to_csv(run/'independent_repeat_genotype_check.tsv',sep='\t',index=False)
result=dict(status='PASS',genotype_rows=388,variants=3926,independently_checked_LD_entries=len(ix)**2,
            exact_reciprocal_anchor_mappings=int(meta.reciprocal_unique.sum()),pedigree_sampling_verified=True,
            mapping_remains_reference_conditional=True)
(run/'genotype_input_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
