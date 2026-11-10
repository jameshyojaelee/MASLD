#!/usr/bin/env python3
"""Recover unique library/barcode RNA counts and conservative existing labels."""
import argparse
import gzip
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits

from input_counts import obs_column
from final_rna_identity import source_rna_libraries
from final_rna_molecule_qc import check_molecule_accounting

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
A = Path(__file__).resolve().parents[1]


def read_matrix(path):
    with h5py.File(path,'r') as h:
        m = h['matrix']
        matrix = sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
        barcodes = m['barcodes'].asstr()[:]
        features = pd.DataFrame(dict(gene_id=m['features/id'].asstr()[:],gene_symbol=m['features/name'].asstr()[:],
                                     feature_type=m['features/feature_type'].asstr()[:]))
    keep = features.feature_type.eq('Gene Expression').to_numpy()
    return matrix[:,keep],barcodes,features[keep].reset_index(drop=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--allow-partial',action='store_true',help='Technical pilot only; no disease fit qualification')
    p.add_argument('--verified-pairing',type=Path,help='Documented source cross-assay key: donor_id,library_id,source_evidence')
    args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(exist_ok=False)
    pairing=pd.read_csv(ROOT/'data/GSE244832/metadata/donor_pairing.csv',dtype=str)
    native=pd.read_csv(ROOT/'Analysis/SingleCell/reviews/control-label-review-20261021/GSE244832_native_diagnosis_comparison.tsv',sep='\t',dtype=str)
    sources=source_rna_libraries(args.candidate)
    assert native.rna_id.is_unique and set(native.rna_id)==set(sources)
    diagnosis=[]
    for row in native.to_dict('records'):
        source=sources[row['rna_id']]
        assert (row['public_donor_id'],row['deposited_diagnosis'],row['gsm']) == (source['source_donor_id'],source['condition'],source['gsm'])
        diagnosis.append(dict(**source,local_crosswalk_condition=row['local_diagnosis'],
            agrees_with_local_crosswalk=row['local_diagnosis']==source['condition'],
            diagnosis_source='deposited_RNA_workbook_checked_against_native_diagnosis_review'))
    pd.DataFrame(diagnosis).to_csv(args.output/'source_diagnosis_decisions.tsv',sep='\t',index=False)
    documented=set()
    if args.verified_pairing:
        key=pd.read_csv(args.verified_pairing,sep='\t',dtype=str,keep_default_na=False)
        assert {'donor_id','library_id','source_evidence'}.issubset(key)
        assert key.source_evidence.str.len().gt(0).all() and not key.duplicated(['donor_id','library_id']).any()
        documented=set(zip(key.donor_id,key.library_id))
    identities=pd.read_csv(A/'results/20261023/inputs/audit/GSE244832.rna_identity.tsv.gz',sep='\t',dtype=str)
    reference=ROOT/'Analysis/ATAC/Human_Multiome/results/label_transfer/rna_reference_GSE244832.h5ad'
    with h5py.File(reference,'r') as h:
        obs=h['obs']; index=obs.attrs.get('_index','_index')
        if isinstance(index,bytes): index=index.decode()
        labels=pd.DataFrame(dict(source_cell_id=obs_column(obs,index),source_label=obs_column(obs,'cell_type')))
    identities=identities.merge(labels,on='source_cell_id',validate='one_to_one')
    from input_counts import LABEL_MAP
    identities['lineage']=identities.source_label.map(LABEL_MAP['GSE244832']).fillna('other')
    # Run labels are annotations of the same nucleus, never extra observations.
    votes=identities.groupby(['rna_library_gsm','raw_barcode']).agg(
        n_source_rows=('lineage','size'),n_broad_labels=('lineage','nunique'),lineage=('lineage','first')).reset_index()
    votes.loc[votes.n_broad_labels != 1,'lineage']='ambiguous'
    votes.to_csv(args.output/'collapsed_label_evidence.tsv.gz',sep='\t',index=False)
    rosters=[]; donors=[]; arrays=[]; genes=None; states=[]
    completed={}
    for marker in (args.candidate/'rna_recovery').glob('*/completion.json'):
        record=json.loads(marker.read_text())
        if record.get('status')=='completed':
            assert record['library'] not in completed,'Multiple successful recounts for one library'
            completed[record['library']]=marker.parent
    if not args.allow_partial:
        assert set(completed)==set(pairing.rna_sample),'All library decisions/counts required before disease inference'
    for row in pairing.to_dict('records'):
        library=row['rna_sample']
        source=sources[library]
        assert row['rna_gsm']==source['gsm']
        if library not in completed:
            states.append(dict(library=library,status='not_completed',reason='no_successful_joint_recount'))
            continue
        directory=completed[library]
        matrix,barcodes,axis=read_matrix(directory/'count/outs/filtered_feature_bc_matrix.h5')
        check_molecule_accounting(directory,barcodes,axis.gene_id.to_numpy(),matrix)
        assert len(set(barcodes))==len(barcodes) and axis.gene_id.is_unique
        assert np.equal(matrix.data,np.floor(matrix.data)).all() and (matrix.data>=0).all()
        if genes is None: genes=axis
        else: assert genes.equals(axis),'Recount gene axes differ'
        umi=np.asarray(matrix.sum(axis=1)).ravel(); detected=np.diff(matrix.indptr)
        mito=genes.gene_symbol.str.startswith('MT-').to_numpy()
        mito_pct=100*np.asarray(matrix[:,mito].sum(axis=1)).ravel()/np.maximum(umi,1)
        cells=pd.DataFrame(dict(raw_barcode=barcodes,library_id=library,donor_id=row['donor_id'],
                               source_donor_id=source['source_donor_id'],specimen=source['specimen'],cohort='GSE244832',condition=source['condition'],
                               n_counts=umi,n_genes=detected,pct_mito=mito_pct))
        cells=cells.merge(votes[votes.rna_library_gsm.eq(row['rna_gsm'])].drop(columns='rna_library_gsm'),on='raw_barcode',how='left',validate='one_to_one',sort=False)
        cells['cell_id']=library+'_'+cells.raw_barcode
        cross_eligible=(row['donor_id'],library) in documented and row['confidence']=='HIGH'
        cells['cross_assay_match_eligible']=cross_eligible
        cells['cross_assay_confidence']=row['confidence']
        cells['same_nucleus_pairing']=False
        cells['local_crosswalk_condition']=row['condition']
        cells['diagnosis_source']='deposited_RNA_workbook_checked_against_native_diagnosis_review'
        cells['qc_pass']=(cells.n_counts>=500)&(cells.n_genes>=200)&(cells.pct_mito<=20)
        cells['lineage']=cells.lineage.fillna('unannotated')
        cells['eligible']=cells.qc_pass&cells.lineage.isin(['hepatocyte','stellate','macrophage'])
        cells['exclusion_reason']=np.where(~cells.qc_pass,'fails_existing_harmonized_count_gene_mito_QC',
            np.where(~cells.lineage.isin(['hepatocyte','stellate','macrophage']),'no_unanimous_existing_target_lineage_annotation',''))
        cells['runs']=';'.join(row['rna_srrs'].split(';'))
        cells['gsm']=row['rna_gsm']; cells['batch_jb']=row['batch_jb']
        for lineage in ['hepatocyte','stellate','macrophage']:
            idx=np.flatnonzero(cells.eligible&cells.lineage.eq(lineage))
            values=np.asarray(matrix[idx].sum(axis=0)).ravel().astype(np.int64)
            donors.append(dict(cohort='GSE244832',donor_id=row['donor_id'],library_id=library,lineage=lineage,
                               condition=source['condition'],local_crosswalk_condition=row['condition'],batch_jb=row['batch_jb'],n_cells=len(idx),
                               total_umis=int(values.sum()),eligible=len(idx)>=100,
                               cross_assay_match_eligible=cross_eligible))
            arrays.append(values)
        rosters.append(cells)
        states.append(dict(library=library,status='qualified_unique_integer_counts',n_called_cells=len(cells),
                           n_eligible_target_cells=int(cells.eligible.sum()),n_source_runs=len(row['rna_srrs'].split(';')),
                           total_umis=int(umi.sum()),count_directory=str(directory)))
    assert genes is not None
    donors=pd.DataFrame(donors); counts=sparse.csr_matrix(np.vstack(arrays))
    assert np.array_equal(np.asarray(counts.sum(axis=1)).ravel(),donors.total_umis)
    with threadpool_limits(limits=1):
        with gzip.open(args.output/'pseudobulk.mtx.gz','wb') as handle: mmwrite(handle,counts.T)
    genes.to_csv(args.output/'genes.tsv',sep='\t',index=False)
    donors.to_csv(args.output/'donors.tsv',sep='\t',index=False)
    pd.concat(rosters,ignore_index=True).to_csv(args.output/'cell_roster.tsv.gz',sep='\t',index=False)
    pd.DataFrame(states).to_csv(args.output/'library_decisions.tsv',sep='\t',index=False)
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed',technical_pilot=args.allow_partial,
        n_libraries=len(completed),n_unique_cells=sum(len(x) for x in rosters),n_genes=len(genes),
        annotation='unanimous_broad_label_across_existing_run_annotations; unannotated_and_discordant_excluded',
        QC='existing_harmonized_min_counts500_min_genes200_mito_le20; prior_atlas_labels_inherit_prior_QC',
        new_doublet_model=False,new_ambient_RNA_correction=False,same_nucleus_pairing=False,scientific_adoption=False),indent=2))


if __name__=='__main__': main()
