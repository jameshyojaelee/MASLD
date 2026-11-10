"""Independent identity, donor-roster, BH and saved GLS checks."""
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
TASK = ROOT / 'Analysis/Hypotheses/NUCLEUS-FIX-20261107'


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p); n = len(p)
    answer = np.empty(n)
    answer[order] = np.minimum(1, np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
    return answer


def main():
    assert os.environ.get('SLURM_JOB_ID')
    completion=json.loads((TASK/'unique_counts/completion.json').read_text())
    assert not completion['technical_pilot']
    with h5py.File(TASK/'unique_counts/unique_counts.h5ad') as handle:
        obs=read_elem(handle['obs']);var=read_elem(handle['var'])
        shape=handle['raw/X'].attrs['shape']
        assert tuple(shape)==(len(obs),len(var))
        assert handle['raw/X/indptr'][-1]==len(handle['raw/X/data'])
        assert handle['raw/X/data'].dtype.kind in 'iu'
    assert len(obs)==completion['n_cells'] and obs.index.is_unique
    keys=obs.capture.astype(str)+'|'+obs.index.str.rsplit('_',n=1).str[-1]
    assert keys.is_unique
    manifest=pd.read_csv(TASK/'capture_manifest.tsv',sep='\t').fillna('')
    excluded=set(';'.join(manifest.excluded_runs).split(';')) - {''}
    assert not set(obs['sample'].astype(str)).intersection(excluded)
    md=pd.read_csv(TASK/'unique_counts/metadata/library_metadata_source_F.tsv',sep='\t')
    census=obs.groupby('sample',observed=True).size()
    assert np.array_equal(md.n_cells,census.reindex(md['sample']))
    rosters=[]; families=[]
    for arm in ['binary','binary_strict','stage']:
        assert list((TASK/arm).glob('*_de.csv')), arm
        for path in sorted((TASK/arm).glob('*_donor_roster.tsv')):
            r=pd.read_csv(path,sep='\t');assert r['sample'].is_unique
            if 'n_cells' in r: assert r.n_cells.ge(50).all()
            rosters.append(dict(arm=arm,file=path.name,n=len(r),unique_donors=True))
        for path in sorted((TASK/arm).glob('*_de.csv')):
            t=pd.read_csv(path)
            assert np.isfinite(t.pvalue).all()
            q=bh(t.pvalue)
            assert np.allclose(q,t.padj,rtol=1e-10,atol=1e-14)
            families.append(dict(arm=arm,file=path.name,genes=len(t),bh_calls=int((q<.05).sum())))
    assert families and rosters
    allcell=pd.read_csv(TASK/'allcell/allcell_pseudobulk_de.csv')
    assert np.allclose(bh(allcell.pvalue),allcell.padj,atol=1e-14)
    assert pd.read_csv(TASK/'allcell/donors_in_model.tsv',sep='\t')['sample'].is_unique
    composition=pd.read_csv(TASK/'composition/lineages_results.tsv',sep='\t')
    assert len(composition)==16 and np.allclose(bh(composition['P.Value']),composition['adj.P.Val'],atol=1e-14)
    assert pd.read_csv(TASK/'composition/lineages_donor_counts.tsv',sep='\t').donor.is_unique
    for path,expected in json.loads((TASK/'checks/preserved_files.json').read_text()).items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    result=dict(unique_cells=len(obs),unique_capture_barcodes=True,
        sample_census_exact=True,frozen_program_files_unchanged=True,
        historical_severity_package_unchanged=True,rosters_checked=len(rosters),families_checked=len(families))
    program_file=TASK/'programs_rerun-20261109T2230Z/results.tsv'
    assert program_file.is_file() and (TASK/'programs_rerun-20261109T2230Z/COMPLETED').is_file()
    if program_file.exists():
        program=pd.read_csv(program_file,sep='\t')
        scores=pd.read_csv(TASK/'fixed_programs/donor_scores.tsv',sep='\t')
        assert len(program)==117 and np.allclose(bh(program.p),program.q,atol=1e-13)
        errors=[]
        for r in program[program.testable].itertuples():
            d=scores[scores.cell_type.eq(r.cell_type)&scores.module.eq(r.module)].sort_values(['dataset','donor'])
            assert d.donor.is_unique and len(d)==r.n
            second=d.dataset.eq('GSE244832').astype(float).to_numpy()
            design=np.column_stack([np.ones(len(d)),d.stage.to_numpy(),second])
            sd=np.where(second,r.sigma_GSE244832,r.sigma_GSE202379)
            precision=1/sd**2
            inv=np.linalg.inv(design.T@(precision[:,None]*design))
            beta=inv@(design.T@(precision*d.score.to_numpy()))
            error=abs(beta[1]-r.beta); errors.append(error)
            assert error<1e-7 and abs(np.sqrt(inv[1,1])-r.se)<1e-7
        result.update(program_BH117=True,GLS_independent_WLS_max_error=max(errors,default=0),
                      program_selected=int(program.selected.sum()))
    pd.DataFrame(rosters).to_csv(TASK/'programs_rerun-20261109T2230Z/rosters.tsv',sep='\t',index=False)
    pd.DataFrame(families).to_csv(TASK/'programs_rerun-20261109T2230Z/families.tsv',sep='\t',index=False)
    (TASK/'programs_rerun-20261109T2230Z/independent_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
