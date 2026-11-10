"""Replace repeated-run cells with joint molecules; retain frozen annotations/weights.

No cell or program discovery is performed. A joint-called barcode is retained
only when it was in the original atlas and all retained run annotations agree.
The original observed-cell selection is therefore inherited, not reoptimized.
Fixed-weight log-CP10k projections are explicitly distinct from native Hotspot
neighbor-smoothed PCA scores, whose graph was learned with repeated cells.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from anndata.io import read_elem, write_elem

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
TASK = ROOT / 'Analysis/Hypotheses/NUCLEUS-FIX-20261107'
OLD = ROOT / 'Analysis/SingleCell/candidates/source-diagnosis-correction-20261026'
MEMBER = ROOT / 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv'
sys.path.insert(0, str(ROOT / 'Analysis/ATAC/Regulatory_Methods/scripts'))
from final_rna_pseudobulk import read_matrix
from final_rna_molecule_qc import check_molecule_accounting


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
    pilot = os.environ.get('NUCLEUS_FIX_PILOT_CAPTURE')
    out = TASK / ('pilot_' + pilot if pilot else 'unique_counts'); out.mkdir(exist_ok=False)
    captures = pd.read_csv(TASK / 'capture_manifest.tsv', sep='\t')
    if pilot:
        captures = captures[captures.capture.eq(pilot)]
        assert len(captures)==1 and captures.iloc[0].n_runs>1
    repeated = captures[captures.n_runs > 1]
    for path in repeated.matrix:
        assert Path(path).is_file(), path
    old = h5py.File(OLD / 'integrated_atlas_source_qualified_v2.h5ad', 'r')
    obs = read_elem(old['obs']); var = read_elem(old['raw/var'])
    assert var.gene_ids.is_unique
    membership = pd.read_csv(MEMBER, sep='\t')
    programs = membership[['program_uid','cell_type','module']].drop_duplicates().sort_values(['cell_type','module'])
    assert len(programs) == 117 and len(membership) == 7093
    assert membership.source_gene.isin(var.index).all()
    pos = pd.Series(np.arange(len(var)), index=var.index)
    pp = pd.Series(np.arange(len(programs)), index=programs.program_uid)
    weights = sparse.csr_matrix((membership.original_l1_weight,
        (membership.source_gene.map(pos), membership.program_uid.map(pp))), shape=(len(var),117))
    assert np.allclose(np.asarray(weights.sum(0)).ravel(),1)
    old_hash = hashlib.sha256(MEMBER.read_bytes()).hexdigest()
    outfile = h5py.File(out / 'unique_counts.h5ad', 'x')
    outfile.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
    raw = outfile.create_group('raw'); raw.attrs.update({'encoding-type':'raw','encoding-version':'0.1.0'})
    group = raw.create_group('X'); group.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0'})
    data = group.create_dataset('data',shape=(0,),maxshape=(None,),dtype='int64',chunks=True,compression='lzf')
    indices = group.create_dataset('indices',shape=(0,),maxshape=(None,),dtype='int32',chunks=True,compression='lzf')
    pointers = [0]; metas=[]; census=[]; sums={}; scores={}; sizes={}; library_checks=[]
    baseline_scores={}; baseline_sizes={}

    def project(matrix, frame, score_totals, cell_totals):
        totals=np.asarray(matrix.sum(1)).ravel()
        norm=sparse.diags(np.divide(1e4,totals,out=np.zeros(len(totals)),where=totals>0)) @ matrix
        norm.data=np.log1p(norm.data)
        projected=(norm @ weights).toarray()
        for (sample,ct), positions in frame.groupby(['sample','cell_type'],observed=True,sort=False).indices.items():
            key=(str(sample),str(ct)); j=np.asarray(positions)
            score_totals[key]=score_totals.get(key,0)+projected[j].sum(0)
            cell_totals[key]=cell_totals.get(key,0)+len(j)

    def consume(matrix, frame):
        matrix=matrix.tocsr(); matrix.eliminate_zeros()
        assert matrix.shape==(len(frame),len(var)) and np.isfinite(matrix.data).all()
        assert (matrix.data>=0).all() and np.equal(matrix.data,np.floor(matrix.data)).all()
        start=len(data); end=start+matrix.nnz
        data.resize((end,)); indices.resize((end,))
        data[start:end]=matrix.data.astype('int64'); indices[start:end]=matrix.indices
        pointers.extend((matrix.indptr[1:]+start).tolist()); metas.append(frame)
        project(matrix,frame,scores,sizes)
        for (sample,ct), positions in frame.groupby(['sample','cell_type'],observed=True,sort=False).indices.items():
            key=(str(sample),str(ct)); j=np.asarray(positions)
            values=np.asarray(matrix[j].sum(0)).ravel().astype('int64')
            sums[key]=sums.get(key,0)+values

    repeat_runs=set(';'.join(repeated.runs).split(';'))
    excluded_runs=set(';'.join(captures.excluded_runs.fillna('')).split(';')) - {''}
    excluded_rows=obs[obs['sample'].astype(str).isin(excluded_runs)]
    excluded_rows.groupby('sample',observed=True).size().rename('excluded_ADT_rows').to_csv(out/'excluded_assay_rows.tsv',sep='\t')
    keep=~obs['sample'].astype(str).isin(repeat_runs | excluded_runs).to_numpy()
    active=obs['sample'].astype(str).isin(repeat_runs).to_numpy() if pilot else np.ones(len(obs),dtype=bool)
    keep &= active
    source=old['raw/X']; ptr=source['indptr'][:]
    for start in range(0,len(obs),5000):
        stop=min(start+5000,len(obs)); mask=keep[start:stop]
        if not active[start:stop].any(): continue
        lo,hi=ptr[start],ptr[stop]
        matrix=sparse.csr_matrix((source['data'][lo:hi],source['indices'][lo:hi],ptr[start:stop+1]-lo),shape=(stop-start,len(var)))
        project(matrix[active[start:stop]],obs.iloc[start:stop].loc[active[start:stop]],baseline_scores,baseline_sizes)
        if not mask.any(): continue
        frame=obs.iloc[start:stop].loc[mask].copy()
        frame['capture']=frame['sample'].astype(str).map(dict(zip(captures.representative_run,captures.capture)))
        assert frame.capture.notna().all()
        consume(matrix[mask],frame)
    for r in repeated.itertuples():
        previous=obs[obs['sample'].astype(str).isin(r.runs.split(';'))].copy()
        previous['barcode']=previous.index.str.split('_',n=1).str[1]
        labels=previous.groupby('barcode',observed=True).cell_type.agg(['first','nunique'])
        matrix,barcodes,features=read_matrix(Path(r.matrix))
        molecule_check=check_molecule_accounting(Path(r.matrix).parents[2],barcodes,features.gene_id.to_numpy(),matrix)
        j=pd.Index(features.gene_id).get_indexer(var.gene_ids)
        assert (j>=0).all(), 'Gene-axis mismatch'
        selected=labels.reindex(barcodes)
        # Cell calling is joint. The annotation/QC inclusion is inherited from
        # the old atlas; conflicting and new barcodes remain explicitly excluded.
        detected=np.diff(matrix[:,j].tocsr().indptr)
        qc=(detected>=200)&(detected<=10000)
        retained=selected['nunique'].eq(1).to_numpy() & qc
        counts=matrix[retained][:,j].tocsr()
        template=previous.iloc[0]
        frame=pd.DataFrame({c:[template[c]]*int(retained.sum()) for c in obs.columns})
        frame.index=pd.Index(r.capture+'_'+barcodes[retained],name=obs.index.name)
        frame['sample']=r.representative_run
        frame['cell_type']=selected.loc[retained,'first'].to_numpy()
        frame['cell_type_raw']=frame.cell_type
        frame['cell_type_conf']=previous.groupby('barcode',observed=True).cell_type_conf.min().reindex(barcodes[retained]).to_numpy()
        frame['leiden']='historical_cluster_unassigned'
        frame['capture']=r.capture
        assert frame.index.is_unique
        census.append(dict(dataset=r.dataset,capture=r.capture,old_rows=len(previous),
            old_unique_barcodes=len(labels),joint_called=len(barcodes),retained=int(retained.sum()),
            ambiguous=int(selected['nunique'].gt(1).sum()),not_in_old_atlas=int(selected['nunique'].isna().sum()),
            joint_gene_QC_excluded=int((selected['nunique'].eq(1).to_numpy() & ~qc).sum()),
            old_not_joint_called=int((~labels.index.isin(barcodes)).sum())))
        for start in range(0,len(frame),5000):consume(counts[start:start+5000],frame.iloc[start:start+5000].copy())
        library_checks.append(dict(capture=r.capture,retained_umis=int(counts.sum()),matrix=str(r.matrix),
            called_cell_umis=molecule_check['called_cell_umis'],molecule_accounting_exact=True))
        print('replaced',r.dataset,r.capture,len(previous),'->',len(frame),flush=True)
    final=pd.concat(metas)
    assert final.index.is_unique and not final[['capture']].isna().any().any()
    identity=final.capture.astype(str)+'|'+final.index.str.rsplit('_',n=1).str[-1]
    assert identity.is_unique, 'Repeated physical nucleus remains'
    group.create_dataset('indptr',data=np.asarray(pointers,dtype='int64'))
    group.attrs['shape']=(len(final),len(var))
    write_elem(raw,'var',var); write_elem(raw,'varm',{})
    outfile['X']=group  # Explicit raw-count atlas; no old normalized X/embedding.
    write_elem(outfile,'var',var)
    for column in final:
        if final[column].dtype==object: final[column]=final[column].fillna('').astype(str).astype('category')
    write_elem(outfile,'obs',final)
    write_elem(outfile,'uns',{'count_state':'unique_capture_integer_counts',
        'embedding_state':'not_computed; historical neighborhoods must not be reused',
        'program_membership_sha256':old_hash})
    outfile.close(); old.close()
    pd.DataFrame(census).to_csv(out/'capture_census.tsv',sep='\t',index=False)
    pd.DataFrame(library_checks).to_csv(out/'count_checks.tsv',sep='\t',index=False)
    metadata=pd.read_csv(OLD/'metadata/library_metadata_source_F.tsv',sep='\t').set_index('sample')
    used=sorted(final['sample'].astype(str).unique()); metadata=metadata.loc[used].copy()
    metadata['n_cells']=final.groupby('sample',observed=True).size().reindex(used).to_numpy()
    for col in metadata.columns:
        if col.startswith('n_') and col not in ('n_cells','n_donor_cells','n_source_records'): metadata[col]=0
    pb=out/'pseudobulk'; pb.mkdir()
    for ct in sorted({k[1] for k in sums}):
        keys=sorted(k for k in sums if k[1]==ct)
        for sample,_ in keys: metadata.loc[sample,'n_'+ct.replace(' ','_')]=sizes[(sample,ct)]
        keys=[k for k in keys if sizes[k]>=10]
        if keys:
            table=pd.DataFrame({k[0]:sums[k] for k in keys},index=var.index)
            table.index.name='gene'; table.to_csv(pb/(ct.replace(' ','_').replace('/','_')+'_pseudobulk.csv'))
    md=out/'metadata'; md.mkdir()
    metadata.reset_index().to_csv(md/'library_metadata_source_F.tsv',sep='\t',index=False)
    # Both standard filenames select the same corrected library census.
    metadata.reset_index().to_csv(md/'library_metadata.tsv',sep='\t',index=False)
    donors=pd.read_csv(OLD/'metadata/donor_metadata.tsv',sep='\t')
    donors.to_csv(md/'donor_metadata.tsv',sep='\t',index=False)
    output=[]
    for (sample,ct), values in scores.items():
        row=metadata.loc[sample]
        for k,program in enumerate(programs.itertuples()):
            output.append(dict(sample=sample,donor=row.donor,dataset=row.dataset,lineage=ct,
                cell_type=program.cell_type,module=program.module,program_uid=program.program_uid,
                n_cells=sizes[(sample,ct)],score_sum=values[k],score=values[k]/sizes[(sample,ct)]))
    pd.DataFrame(output).to_csv(out/'fixed_weight_projection_by_capture.tsv.gz',sep='\t',index=False)
    original_metadata=pd.read_csv(OLD/'metadata/library_metadata_source_F.tsv',sep='\t').set_index('sample')
    baseline=[]
    for (sample,ct),values in baseline_scores.items():
        row=original_metadata.loc[sample]
        for k,program in enumerate(programs.itertuples()):
            baseline.append(dict(sample=sample,donor=row.donor,dataset=row.dataset,lineage=ct,
                cell_type=program.cell_type,module=program.module,program_uid=program.program_uid,
                n_cells=baseline_sizes[(sample,ct)],score_sum=values[k],score=values[k]/baseline_sizes[(sample,ct)]))
    pd.DataFrame(baseline).to_csv(out/'historical_fixed_weight_projection_by_run.tsv.gz',sep='\t',index=False)
    assert hashlib.sha256(MEMBER.read_bytes()).hexdigest()==old_hash
    (out/'completion.json').write_text(json.dumps(dict(n_cells=len(final),n_features=len(var),
        n_captures=final.capture.nunique(),n_donors=final.donor.nunique(),programs=117,
        source_membership_rows=len(membership),membership_sha256=old_hash,raw_integer_counts=True,
        unique_capture_barcode=True,new_embedding=False,
        projection_definition='fixed L1 weights applied to log1p(counts/total*10000); not native Hotspot PCA',
        inherited_QC_and_unanimous_labels=True,technical_pilot=bool(pilot)),indent=2))


if __name__=='__main__': main()
