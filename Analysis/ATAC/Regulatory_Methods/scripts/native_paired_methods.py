#!/usr/bin/env python3
"""Bounded native TRIPOD stages for genuine GSE296875 paired raw counts.

The shared disjoint donor-confined 50-cell aggregation is an adaptation, not the
native Seurat aggregation. Training-stage links need separate donor-held-out
prediction evaluation. Native matched-test p-values are not donor significance.
RNA-only and prior-only native TRIPOD/LINGER interfaces are unsupported here.
LINGER atlas files must be separately qualified; scNN is not substituted.

Official sources: https://github.com/yuchaojiang/TRIPOD/tree/
c8421358f6bfddeb6b219cbc836a40aae9d20016/package/R
LINGER source inspected: PyPI LingerGRN1.110 wheel SHA256
2db89d2953b12c79c5c9076fda53def08f956d70c1fe63eb31f9494177c7f2ec.
"""
import argparse
import ast
import copy
import hashlib
import inspect
import json
import os
import shutil
from pathlib import Path
import subprocess
import textwrap
import numpy as np
import pandas as pd
from scipy import sparse
import network_compare as shared


LINGER_GLOBAL = ['bulk_gene_all.txt','TFName.txt','Match2.txt','all_hg19.txt',
                 'motifWeight.txt','hg38_hg19_pair.bed','RE_gene_corr_hg19.bed']
LINGER_CHROM = ['MotifTarget_Matrix_{c}.txt','MotifTarget_matrix_{c}.bed',
                'hg38_Peaks_{c}.bed','{c}_gene.txt','{c}_index.txt',
                '{c}_index_all.txt','fisher_{c}.pt','all_models_{c}.pt']


def linger_inventory(root, chromosome):
    names = LINGER_GLOBAL + [x.format(c=chromosome) for x in LINGER_CHROM]
    return {'status':'not_run', 'reason':'atlas_compatibility_and_prior_overlap_unqualified',
            'files':[{'name':n,'exists':bool(root and (root/n).is_file())} for n in names],
            'remaining_requirements':['version1.110-compatible atlas manifest and checksums',
              'prior training contexts and evaluation-reference overlap qualification',
              'bounded native target restriction and CPU scheduling qualification'],
            'controls':{'rna_only':'unsupported','prior_only':'unsupported',
                        'pairing_shuffle':'not_run_native_atlas_missing'}}


def linger_model_preflight(args):
    """Inspect genuine authorized atlas serialization, never synthesize priors."""
    shared.require_compute()
    import torch
    from LingerGRN import LINGER_tr
    torch.set_num_threads(1)
    result=linger_inventory(args.prior_root,args.chromosome)
    result.update(prior_root=str(args.prior_root.resolve()),chromosome=args.chromosome)
    missing=[x['name'] for x in result['files'] if not x['exists']]
    if missing:
        result['reason']='missing_required_atlas_files'; result['missing']=missing
    else:
        prefix=args.prior_root
        tables={name:pd.read_csv(prefix/name,sep='\t',header=None)
                for name in [args.chromosome+'_gene.txt',args.chromosome+'_index.txt',args.chromosome+'_index_all.txt']}
        n=len(tables[args.chromosome+'_gene.txt'])
        if len(tables[args.chromosome+'_index.txt'])!=n or len(tables[args.chromosome+'_index_all.txt'])!=n:
            raise ValueError('Atlas gene/index row axes disagree')
        models={}; state_dicts=None; fishers=None
        for stem in ['fisher_','all_models_']:
            path=prefix/(stem+args.chromosome+'.pt')
            if path.stat().st_size>args.max_model_bytes:
                raise MemoryError('Bounded pilot serialization exceeds explicit byte limit')
            # Official user-approved atlas only: pickle loading can execute code.
            obj=torch.load(path,map_location='cpu',weights_only=False)
            models[stem]={'object_type':type(obj).__name__,'n_objects':len(obj),
                         'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
            if stem=='all_models_': state_dicts=obj
            else: fishers=obj
        checked=[]
        indices=tables[args.chromosome+'_index.txt']
        chosen=list(range(min(n,3)))
        if args.target_symbols:
            names=tables[args.chromosome+'_gene.txt'].iloc[:,0].astype(str)
            chosen=[]
            for symbol in args.target_symbols.split(','):
                found=np.flatnonzero(names.eq(symbol))
                if len(found)!=1: raise ValueError('Target not uniquely observable in atlas: '+symbol)
                chosen.append(int(found[0]))
        for index in chosen:
            model=LINGER_tr.Net(int(indices.iloc[index,3]),'ReLU')
            model.load_state_dict(state_dicts[index],strict=True)
            fisher=fishers[index][0]
            if tuple(fisher.shape)!=tuple(model.fc1.weight.shape):
                raise ValueError('Fisher/first-layer parameter axes disagree')
            checked.append({'gene_row0':index,'symbol':str(tables[args.chromosome+'_gene.txt'].iloc[index,0]),
                            'input_size':int(indices.iloc[index,3])})
        result.update(status='structurally_compatible_not_fitted',
          reason='native_import_and_serialization_passed_not_biological_qualification',
          n_genes=n,models=models,checked_gene_rows=checked,torch_version=torch.__version__,
          native_module=LINGER_tr.__file__)
        result['controls']['pairing_shuffle']='supported_prepared_stage_native_fit_not_run'
    with args.output.open('x') as h: json.dump(result,h,indent=2)


def stages(donors, seed):
    outer = shared.donor_folds(donors, 5, seed)
    for fold in range(5):
        train = np.flatnonzero(outer != fold)
        yield f'outer{fold}', train, set(np.asarray(donors)[outer == fold])
        inner = shared.donor_folds(np.asarray(donors)[train],3,shared.stable_seed(seed,('inner',fold)))
        for j in range(3):
            yield f'outer{fold}_inner{j}', train[inner != j], set(np.asarray(donors)[train[inner == j]]) | set(np.asarray(donors)[outer == fold])


def qualified_context(cells):
    if set(cells.cohort) != {'GSE296875'} or len(cells[['cohort','lineage']].drop_duplicates())!=1:
        raise ValueError('Require exactly one genuine GSE296875 cohort/lineage per invocation')
    return {field:';'.join(sorted(set(cells[field])))
            for field in ['cohort','lineage','lineage_state','source_lineage'] if field in cells}


def native_tripod_cpm(matrix):
    """TRIPOD prep.R152-173: aggregate counts, full-assay CPM1e6, no log."""
    totals=np.asarray(matrix.sum(axis=1)).ravel()
    if (totals<=0).any(): raise ValueError('Cannot form native CPM from zero-depth aggregate')
    return (sparse.diags(1e6/totals)@matrix).tocsr()


def linger_mean_inputs(raw_rna,raw_atac,assignment):
    normalized=shared.normalize_sparse(raw_rna); normalized.data=np.log1p(normalized.data)
    transformed_atac=raw_atac.copy(); transformed_atac.data=np.log1p(transformed_atac.data)
    return (assignment@normalized).tocsr(),(assignment@transformed_atac).tocsr()


def prepare_linger_stage(args):
    """Native per-nucleus units, with disjoint donor aggregation substituted.

    Source pseudo_bulk.py39-51,79-83: RNA library-normalize/log per nucleus;
    ATAC log raw counts without depth normalization; average transformed values.
    Its min3-nucleus filter is learned on training nuclei only. Full original
    RNA axes define library totals before any gene-symbol/target restriction.
    """
    shared.require_compute()
    counts,cells,peaks=shared.load_input(args.paired_root,args.seed,100)
    context=qualified_context(cells)
    if any(not (args.paired_root/name).is_file() for name in ['paired_export.json','RNA_export.json']):
        raise ValueError('Genuine paired raw export records required')
    aggregate,meta,groups,_=shared.build_metacells(counts,cells,args.seed,50,100)
    chosen=np.concatenate(groups)
    if len(set(chosen))!=len(chosen) or set(meta.groupby('donor_id').size())!={2}:
        raise ValueError('Require disjoint two50-cell aggregates per donor')
    annotation=shared.read_table(args.genes).set_index('gene_id')
    _,idx=shared.gene_peak_indices(peaks,annotation.loc[args.target])
    tss=int(annotation.loc[args.target,'tss0'])
    idx=sorted(idx,key=lambda i:(abs((peaks.iloc[i].start0+peaks.iloc[i].end)/2-tss),peaks.iloc[i].peak_id))[:25]
    raw_atac=counts[chosen][:,idx].tocsr()
    del counts,aggregate
    rna_genes=shared.read_table(args.paired_root/'genes.tsv')
    raw_rna=shared.read_matrix(args.paired_root,['rna_counts.npz','RNA_counts.mtx.gz'])
    if raw_rna.shape!=(len(cells),len(rna_genes)): raise ValueError('RNA axes disagree')
    raw_rna=raw_rna[chosen].tocsr()
    assignment=sparse.csr_matrix((np.full(len(chosen),1/50),
      (np.repeat(np.arange(len(groups)),50),np.arange(len(chosen)))),shape=(len(groups),len(chosen)))
    mean_rna,mean_atac=linger_mean_inputs(raw_rna,raw_atac,assignment)
    options={name:(train,excluded) for name,train,excluded in stages(meta.donor_id.to_numpy(),args.seed)}
    if args.stage!='all' and args.stage not in options: raise ValueError('Unknown donor training stage')
    args.output.mkdir(parents=True,exist_ok=False)
    for stage_name in options if args.stage=='all' else [args.stage]:
        train,excluded=options[stage_name]
        destination=args.output/stage_name if args.stage=='all' else args.output
        destination.mkdir(parents=True,exist_ok=True)
        write_linger_inputs(args,stage_name,train,excluded,context,meta,raw_rna,
                            raw_atac,mean_rna,mean_atac,rna_genes,peaks,idx,destination)


def write_linger_inputs(args,stage_name,train,excluded,context,meta,raw_rna,
                       raw_atac,mean_rna,mean_atac,rna_genes,peaks,idx,destination):
    nucleus_rows=np.concatenate([np.arange(i*50,(i+1)*50) for i in train])
    unique_symbol=~rna_genes.gene_symbol.duplicated(keep=False)&rna_genes.gene_symbol.ne('')
    gene_keep=unique_symbol.to_numpy()&(np.asarray((raw_rna[nucleus_rows]>0).sum(axis=0)).ravel()>=3)
    peak_keep=np.asarray((raw_atac[nucleus_rows]>0).sum(axis=0)).ravel()>=3
    if args.target not in set(rna_genes.loc[gene_keep,'gene_id']): raise ValueError('Native target unobservable')
    y=mean_rna[train][:,gene_keep].toarray(); x=mean_atac[train][:,peak_keep].toarray()
    ids=rna_genes.loc[gene_keep,'gene_id'].tolist(); pids=peaks.iloc[idx].loc[peak_keep,'peak_id'].tolist()
    for control in ['paired','pairing_shuffle']:
        path=destination/control; path.mkdir()
        yy=y.copy()
        if control=='pairing_shuffle':
            rng=np.random.default_rng(shared.stable_seed(args.seed,args.target+stage_name))
            for donor in sorted(set(meta.iloc[train].donor_id)):
                rows=np.flatnonzero(meta.iloc[train].donor_id.to_numpy()==donor)
                yy[rows]=yy[rng.permutation(rows)]
        pd.DataFrame(yy,columns=ids).to_csv(path/'rna.tsv',sep='\t',index=False)
        pd.DataFrame(x,columns=pids).to_csv(path/'atac.tsv',sep='\t',index=False)
        mm=meta.iloc[train].copy(); mm['split']='training'; mm.to_csv(path/'metacells.tsv',sep='\t',index=False)
        info=dict(context,gene_id=args.target,stage=stage_name,control=control,
          excluded_donors=sorted(excluded),status='prepared_not_fitted',maximum_nuclei_per_donor=100,
          input_units='LINGER_per_nucleus_transform_then_mean',
          RNA_units='mean_log1p_RNA_UMI_per_10000_full_original_gene_axis',
          ATAC_units='mean_log1p_raw_native_peak_count_no_depth_normalization',
          input_genes=len(ids),input_peaks=len(pids),
          native_min_cells_filter='at_least3_training_nuclei',
          ambiguous_symbol_features_excluded=int((~unique_symbol).sum()),
          aggregation='adapted_disjoint_donor_confined50_nuclei',native_donor_weighting='unweighted_metacells')
        (path/'stage.json').write_text(json.dumps(info,indent=2)+'\n')


def prepare(args):
    shared.require_compute()
    if args.cell_cap not in {100,1000}:
        raise ValueError('Bounded native comparison requires cap100 or cap1000 per donor')
    atac,cells,peaks = shared.load_input(args.paired_root, args.seed, args.cell_cap)
    context=qualified_context(cells)
    for marker in ['paired_export.json','RNA_export.json']:
        if not (args.paired_root/marker).is_file():
            raise ValueError('Missing raw paired export marker: '+marker)
    rna = shared.read_matrix(args.paired_root,['rna_counts.npz','RNA_counts.mtx.gz'])
    genes = shared.read_table(args.paired_root/'genes.tsv')
    if rna.shape != (len(cells),len(genes)) or genes.gene_id.duplicated().any():
        raise ValueError('RNA paired-axis mismatch')
    if not np.isfinite(rna.data).all() or (rna.data<0).any() or not np.equal(rna.data,np.floor(rna.data)).all():
        raise ValueError('Expected genuine nonnegative integer RNA counts')
    aggregate,meta,groups,_ = shared.build_metacells(atac,cells,args.seed,50,args.cell_cap)
    if len(set(meta.donor_id))<5: raise ValueError('Need >=5 eligible donors')
    assignment = sparse.csr_matrix((np.ones(sum(map(len,groups))),
      (np.repeat(np.arange(len(groups)),[len(g) for g in groups]),np.concatenate(groups))),shape=(len(groups),len(cells)))
    arna = (assignment@rna).tocsr()
    # Shared native-assay normalization; only bounded columns are made dense.
    aggregate = native_tripod_cpm(aggregate); arna = native_tripod_cpm(arna)
    ann = shared.read_table(args.genes); ann.index=ann.gene_id
    if ann.index.duplicated().any(): raise ValueError('Duplicate target annotation')
    targets=args.targets.split(',')
    if len(targets)>5: raise ValueError('Pilot bounded to <=5 prespecified targets')
    motif=shared.read_table(args.motif_map)
    required={'peak_id','motif_id','tf_gene_id','mapping_state'}
    if not required.issubset(motif): raise ValueError('Exact motif map requires '+str(required))
    if not set(motif.peak_id).issubset(set(peaks.peak_id)):
        raise ValueError('Motif map contains non-native peaks')
    if not set(motif.mapping_state).issubset({'exact_symbol','complex_component'}):
        raise ValueError('Unqualified motif-to-TF mapping')
    lookup=dict(zip(genes.gene_id,range(len(genes))))
    args.output.mkdir(parents=True,exist_ok=False)
    sources={str(path.resolve()):hashlib.sha256(path.read_bytes()).hexdigest()
             for path in [args.genes,args.motif_map,args.paired_root/'cells.tsv',args.paired_root/'genes.tsv',args.paired_root/'peaks.tsv']}
    (args.output/'source_sha256.json').write_text(json.dumps(sources,indent=2)+'\n')
    records=[]
    for target in targets:
        if target not in lookup or target not in ann.index: raise ValueError('Unknown target '+target)
        _,idx=shared.gene_peak_indices(peaks,ann.loc[target])
        # Geometric ordering is prespecified; no held-out expression used.
        idx=sorted(idx,key=lambda i:(abs((peaks.iloc[i].start0+peaks.iloc[i].end)/2-int(ann.loc[target,'tss0'])),peaks.iloc[i].peak_id))[:25]
        pids=peaks.iloc[idx].peak_id.tolist()
        m=motif[motif.peak_id.isin(pids)&motif.tf_gene_id.isin(lookup)].copy()
        tfs=sorted(set(m.tf_gene_id)-{target})[:10]
        m=m[m.tf_gene_id.isin(tfs)]
        if not len(m):
            records.append({'gene_id':target,'status':'not_run','reason':'no_exact_mapped_cis_motif'}); continue
        ids=[target]+tfs
        x=aggregate[:,idx].toarray(); y=arna[:,[lookup[g] for g in ids]].toarray()
        for stage,train,excluded in stages(meta.donor_id.to_numpy(),args.seed):
            if set(meta.iloc[train].donor_id)&excluded: raise AssertionError('Donor leakage')
            for control in ['paired','pairing_shuffle']:
                path=args.output/target/stage/control; path.mkdir(parents=True)
                yy=y[train].copy()
                if control=='pairing_shuffle':
                    rng=np.random.default_rng(shared.stable_seed(args.seed,target+stage))
                    for donor in sorted(set(meta.iloc[train].donor_id)):
                        rows=np.flatnonzero(meta.iloc[train].donor_id.to_numpy()==donor)
                        yy[rows]=yy[rng.permutation(rows)]
                pd.DataFrame(x[train],columns=pids).to_csv(path/'atac.tsv',sep='\t',index=False)
                pd.DataFrame(yy,columns=ids).to_csv(path/'rna.tsv',sep='\t',index=False)
                mm=meta.iloc[train].copy(); mm['split']='training'; mm.to_csv(path/'metacells.tsv',sep='\t',index=False)
                m.to_csv(path/'motifs.tsv',sep='\t',index=False)
                pd.DataFrame({'gene_id':[target]}).to_csv(path/'target.tsv',sep='\t',index=False)
                info={'gene_id':target,'stage':stage,'control':control,'excluded_donors':sorted(excluded),
                      'aggregation':'adapted_disjoint_50_cell_within_donor','status':'prepared_not_fitted',
                      'rna_only':'unsupported_native_interface','prior_only':'unsupported_native_interface'}
                info.update(context)
                info['input_units']='TRIPOD_full_axis_CPM1e6_no_log'
                info['native_donor_weighting']='unweighted_metacells_and_matched_pairs'
                info['maximum_nuclei_per_donor']=args.cell_cap
                info['analysis_context']='donor_balanced_native_comparison_100_nuclei_per_donor' if args.cell_cap==100 else 'bounded_native_development_diagnostic'
                (path/'stage.json').write_text(json.dumps(info,indent=2)+'\n'); records.append(info)
    (args.output/'stages.json').write_text(json.dumps(records,indent=2)+'\n')


def fit_batch(args):
    """Batch native fits, publishing only complete five-by-three stage families."""
    shared.require_compute()
    records=json.loads((args.stages/'stages.json').read_text())
    records=[r for r in records if r.get('status')=='prepared_not_fitted' and r['control']==args.control]
    if not records: raise ValueError('No prepared stages')
    args.output.mkdir(parents=True,exist_ok=False)
    snapshot=args.output/'source_snapshot'; snapshot.mkdir()
    for path in [Path(__file__),Path(__file__).with_suffix('.R')]:
        shutil.copy2(path,snapshot/path.name)
    jobs=[]; training_counts=[]
    for r in records:
        if r.get('input_units')!='TRIPOD_full_axis_CPM1e6_no_log':
            raise ValueError('Native TRIPOD requires full-axis CPM1e6 without log')
        stage=args.stages/r['gene_id']/r['stage']/r['control']
        meta=shared.read_table(stage/'metacells.tsv'); qualified_context(meta)
        training_counts.extend(meta.groupby('donor_id').size().tolist())
        if set(meta.donor_id)&set(r['excluded_donors']): raise ValueError('Donor leakage')
        jobs.append({'stage':str(stage.resolve()),'output':str((args.output/'fits'/r['gene_id']/r['stage']).resolve())})
    jobfile=args.output/'batch.tsv'; pd.DataFrame(jobs).to_csv(jobfile,sep='\t',index=False)
    subprocess.run([str(args.rscript),str(Path(__file__).with_suffix('.R')),str(jobfile)],check=True,timeout=args.timeout)
    method='TRIPOD_Xt_L1' if args.control=='paired' else 'TRIPOD_Xt_L1_pairing_shuffle'
    config={'external_stage_links':{method:{}},'external_stage_roots':{method:str(args.output.resolve())}}
    failures=[]
    for stage in sorted(set(r['stage'] for r in records)):
        group=[(r,j) for r,j in zip(records,jobs) if r['stage']==stage]
        if any(not (Path(j['output'])/'native_fit.json').is_file() for r,j in group):
            failures.append({'stage':stage,'reason':'native_fit_failed'}); continue
        rows=[]; donors=None
        for r,j in group:
            meta=shared.read_table(Path(j['stage'])/'metacells.tsv')
            actual=set(meta.donor_id)
            if donors is not None and donors!=actual: raise ValueError('Gene-specific donor mismatch')
            donors=actual
            table=pd.read_csv(Path(j['output'])/'links.tsv',sep='\t')
            table['method']=method; rows.append(table)
        out=args.output/'combined'/stage; out.mkdir(parents=True)
        pd.concat(rows,ignore_index=True).to_csv(out/'links.tsv',sep='\t',index=False)
        (out/'links.training_donors.json').write_text(json.dumps({'donor_ids':sorted(donors)},indent=2)+'\n')
        (out/'completion.json').write_text(json.dumps({'status':'completed','scope':'training_links_only','stage':stage})+'\n')
        outer=stage.split('_')[0].replace('outer','')
        key=outer+(':'+stage.split('_')[1] if '_' in stage else ':outer')
        config['external_stage_links'][method][key]=str((out/'links.tsv').resolve())
    complete=len(config['external_stage_links'][method])==20 and not failures
    (args.output/'external_stage_config.json').write_text(json.dumps(config,indent=2)+'\n')
    (args.output/'completion.json').write_text(json.dumps({'status':'completed' if complete else 'incomplete',
      'scope':'native_training_links_not_donor_evaluation','failures':failures,'control':args.control,
      'native_donor_weighting':'unweighted_metacells_and_matched_pairs',
      'rna_only':'unsupported','prior_only':'unsupported'},indent=2)+'\n')
    counts=pd.Series(training_counts)
    (args.output/'native_scope.json').write_text(json.dumps({
      'input_units':'TRIPOD_full_axis_CPM1e6_no_log',
      'native_donor_weighting':'unweighted_metacells_and_matched_pairs',
      'analysis_context':'donor_balanced_native_comparison_100_nuclei_per_donor' if set(counts)=={2} else 'bounded_native_development_diagnostic',
      'metacells_per_training_donor_min':int(counts.min()),
      'metacells_per_training_donor_max':int(counts.max()),
      'disease_inference':'not_qualified','native_p_values':'not_donor_calibrated',
      'aggregation':'adapted_disjoint_seeded_50_nuclei_within_donor'},indent=2)+'\n')


def linger_background_indices(n,requested=50,replace=False):
    if replace or n<1 or requested<1: raise ValueError('Background requires distinct training rows')
    return np.random.choice(n,min(n,requested),replace=False)


def validate_linger_background(meta,indices,excluded):
    indices=np.asarray(indices,dtype=int)
    if len(set(indices))!=len(indices) or (indices<0).any() or (indices>=len(meta)).any():
        raise ValueError('Background contains duplicate or nontraining rows')
    selected=meta.iloc[indices]
    if set(selected.donor_id)&set(excluded): raise ValueError('Heldout donor in SHAP background')
    if selected.metacell_id.duplicated().any(): raise ValueError('Duplicated background aggregate')
    return {'row_indices':indices.tolist(),'metacell_ids':selected.metacell_id.tolist(),
            'donor_ids':selected.donor_id.tolist(),'distinct_training_rows':True,
            'heldout_donors_present':False}


def validate_linger_prior_features(indices,bulk_tf,bulk_re,observed_tf,observed_re,tf_names,re_names,input_size):
    axis=[('TF',int(i),tf_names[int(i)]) for i in bulk_tf]+[('RE',int(i),re_names[int(i)]) for i in bulk_re]
    expected=[('TF',int(i),tf_names[int(i)]) for i in observed_tf]+[('RE',int(i),re_names[int(i)]) for i in observed_re]
    if len(axis)!=input_size or len(indices)!=len(expected): raise ValueError('Prior feature axis length mismatch')
    rows=[]
    for position,(column,wanted) in enumerate(zip(indices,expected)):
        column=int(column)
        if not 0<=column<len(axis): raise ValueError('Prior column outside full feature axis')
        actual=axis[column]
        if actual!=wanted: raise ValueError(f'Prior identity mismatch: input{position}, expected{wanted}, actual{actual}')
        rows.append({'input_position0':position,'prior_column0':column,'feature_type':wanted[0],
          'global_prior_id1':wanted[1],'expected_identity':wanted[2],'actual_identity':actual[2],'match':True})
    return rows


def clone_linger_background_adaptation(module,meta,info,prior_root,cap_background=True):
    """Clone sc_nn with explicit prior-identity and SHAP-background repairs.

    No installed file is modified. Native seed functions are observationally
    wrapped; they still set native42 before NN initialization and background.
    """
    if module.SEED!=42 or module.GLOBAL_SEED!=42: raise ValueError('Unexpected native seed defaults')
    original=textwrap.dedent(inspect.getsource(module.sc_nn))
    tree=ast.parse(original); original_tree=copy.deepcopy(tree); found=[]; offsets=[]; audit_blocks=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='background' for t in node.targets):
            if not isinstance(node.value,ast.Subscript) or not isinstance(node.value.slice,ast.Call): continue
            call=node.value.slice
            if ast.unparse(call.func)!='np.random.choice' or ast.unparse(call.args[0])!='X_tr.shape[0]' or not isinstance(call.args[1],ast.Constant) or call.args[1].value!=50:
                raise ValueError('Native background implementation differs from qualified1.110')
            if len(call.keywords)!=1 or call.keywords[0].arg!='replace' or call.keywords[0].value.value is not False:
                raise ValueError('Native background replacement policy differs')
            call.func=ast.copy_location(ast.Name(id='_linger_background_choice',ctx=ast.Load()),call.func)
            found.append(call)
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='indexall' for t in node.targets):
            for sub in ast.walk(node.value):
                if isinstance(sub,ast.Subscript) and ast.unparse(sub)=='merge_TF.shape[0]':
                    offsets.append(sub)
            if offsets and any(sub is offsets[-1] for sub in ast.walk(node.value)):
                class RepairOffset(ast.NodeTransformer):
                    def visit_Subscript(self,sub):
                        if ast.unparse(sub)=='merge_TF.shape[0]':
                            return ast.copy_location(ast.parse('len(TFidxtemp_b)',mode='eval').body,sub)
                        return self.generic_visit(sub)
                node.value=RepairOffset().visit(node.value)
        if isinstance(node,ast.If) and ast.unparse(node.test)=='good == 1':
            hook=ast.parse("_linger_prior_alignment(indexall,TFidxtemp_b,REidxtemp_b,merge_TF['TF'].tolist(),REidxtemp_b_m,input_size)").body[0]
            node.body.insert(0,hook); audit_blocks.append(hook)
    if len(found)!=1: raise ValueError('Require exactly one native SHAP background selector')
    if len(offsets)!=1 or len(audit_blocks)!=1: raise ValueError('Require exactly one qualified prior offset and observer')
    reversed_tree=copy.deepcopy(tree)
    class ReverseObservers(ast.NodeTransformer):
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='_linger_prior_alignment': return None
            return self.generic_visit(node)
        def visit_Assign(self,node):
            if any(isinstance(t,ast.Name) and t.id=='indexall' for t in node.targets):
                class ReverseOffset(ast.NodeTransformer):
                    def visit_Call(self,call):
                        if ast.unparse(call)=='len(TFidxtemp_b)': return ast.parse('merge_TF.shape[0]',mode='eval').body
                        return self.generic_visit(call)
                node.value=ReverseOffset().visit(node.value)
            return self.generic_visit(node)
    reversed_tree=ReverseObservers().visit(reversed_tree)
    for node in ast.walk(reversed_tree):
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='_linger_background_choice':
            node.func=ast.parse('np.random.choice',mode='eval').body
    if ast.dump(reversed_tree,include_attributes=False)!=ast.dump(original_tree,include_attributes=False):
        raise ValueError('Unexpected native NN equation change outside qualified repairs')
    audit={'stage':info['stage'],'gene_id':info['gene_id'],'control':info['control'],
      'native_SEED':42,'native_GLOBAL_SEED':42,'python_hash_seed_at_process_start':os.environ.get('PYTHONHASHSEED'),
      'native_source_function_sha256':hashlib.sha256(original.encode()).hexdigest(),
      'adapted_ast_sha256':hashlib.sha256(ast.dump(tree,include_attributes=False).encode()).hexdigest(),
      'NN_function_changes':['RE prior offset uses bulk TF count, preserving feature identities',
        'SHAP background uses min(50,n_training) without replacement' if cap_background else 'native50 SHAP background'],
      'native_placeholder_accuracy_values':'hardcoded0.5/0.5; unusable; heldout readout required',
      'seed_calls':[],'backgrounds':[],'prior_alignments':[]}
    tf_names={i+1:name for i,name in enumerate(pd.read_csv(prior_root/'TFName.txt',header=None)[0])}
    re_names={}
    def alignment(indices,bulk_tf,bulk_re,observed_tf,observed_re,input_size):
        missing=set(map(int,bulk_re))-set(re_names)
        if missing:
            last_needed=max(missing)
            with (prior_root/'all_hg19.txt').open() as f:
                for i,line in enumerate(f,1):
                    if i in missing: re_names[i]=line.rstrip('\n\r')
                    if i>=last_needed: break
        rows=validate_linger_prior_features(indices,bulk_tf,bulk_re,observed_tf,observed_re,tf_names,re_names,input_size)
        audit['prior_alignments'].append({'bulk_TFs':len(bulk_tf),'bulk_REs':len(bulk_re),
          'observed_TFs':len(observed_tf),'observed_RE_inputs':len(observed_re),'features':rows})
    module._linger_prior_alignment=alignment
    def wrap_seed(fn,name):
        default=inspect.signature(fn).parameters['seed'].default
        def observed(seed=default):
            result=fn(seed)
            audit['seed_calls'].append({'function':name,'seed':int(seed),
              'numpy_state_sha256':hashlib.sha256(np.random.get_state()[1].tobytes()).hexdigest(),
              'torch_initial_seed':int(module.torch.initial_seed()),
              'torch_state_sha256':hashlib.sha256(module.torch.get_rng_state().numpy().tobytes()).hexdigest()})
            return result
        return observed
    module.set_all_seeds=wrap_seed(module.set_all_seeds,'set_all_seeds')
    module.set_global_seed=wrap_seed(module.set_global_seed,'set_global_seed')
    def choice(n,requested,replace=False):
        if n!=len(meta): raise ValueError('Native NN background population differs from training metadata')
        if int(np.random.get_state()[1][0])!=42 or module.torch.initial_seed()!=42:
            raise ValueError('Native NN/background seed differs from qualified42')
        indices=linger_background_indices(n,requested,replace) if cap_background else np.random.choice(n,requested,replace=replace)
        record=validate_linger_background(meta,indices,info['excluded_donors'])
        record.update(n_training=n,requested_size=requested,actual_size=len(indices),numpy_seed=42,torch_seed=42)
        audit['backgrounds'].append(record)
        return indices
    module._linger_background_choice=choice
    exec(compile(ast.fix_missing_locations(tree),'<LINGER1.110_training_background_adaptation>','exec'),module.__dict__)
    return audit


def fit_linger_stage(args):
    """Native1.110 fit with explicit bounded TF/peak/response adaptations."""
    shared.require_compute()
    qualification=json.loads(args.model_preflight.read_text())
    if qualification.get('status')!='structurally_compatible_not_fitted':
        raise ValueError('Native atlas model preflight has not passed')
    if qualification.get('prior_root')!=str(args.prior_root.resolve()) or qualification.get('chromosome')!=args.chromosome:
        raise ValueError('Preflight prior root/chromosome differs')
    info=json.loads((args.stage/'stage.json').read_text())
    if info.get('input_units')!='LINGER_per_nucleus_transform_then_mean':
        raise ValueError('LINGER requires separate per-nucleus native-unit input preparation')
    meta=shared.read_table(args.stage/'metacells.tsv'); context=qualified_context(meta)
    if set(meta.donor_id)&set(info['excluded_donors']): raise ValueError('Donor leakage')
    rna=pd.read_csv(args.stage/'rna.tsv',sep='\t'); atac=pd.read_csv(args.stage/'atac.tsv',sep='\t')
    if len(rna)>1000 or atac.shape[1]>25: raise ValueError('Pilot size exceeded')
    genes=shared.read_table(args.genes); symbols=dict(zip(genes.gene_id,genes.gene_symbol))
    selected=[symbols[g] for g in rna.columns]; target=symbols[info['gene_id']]
    if target not in {r['symbol'] for r in qualification.get('checked_gene_rows',[])}:
        raise ValueError('Actual target model/index has not passed compatibility check')
    if len(set(selected))!=len(selected) or any(not x for x in selected):
        raise ValueError('Ambiguous stable-gene-to-symbol mapping')
    symbol_counts=genes.gene_symbol.value_counts()
    if any(symbol_counts.get(symbol,0)!=1 for symbol in selected):
        raise ValueError('Gene symbols are not unique in authoritative RNA feature annotation')
    peaks=shared.read_table(args.peaks).set_index('peak_id')
    # Native source IDs encode GRanges1-based starts. Convert from declared
    # BED columns, never by parsing the native identifier.
    p=peaks.loc[atac.columns]
    native=dict(zip(p.chrom+':'+p.start0.astype(str)+'-'+p.end.astype(str),atac.columns))
    if len(native)!=atac.shape[1] or set(p.chrom)!={args.chromosome}:
        raise ValueError('Peak identity/chromosome disagreement')
    os.environ['CUDA_VISIBLE_DEVICES']=''; os.environ['SLURM_CPUS_PER_TASK']='1'
    import torch
    import joblib
    import pybedtools
    from LingerGRN import preprocess as pre, LINGER_tr as train, LL_net
    torch.set_num_threads(1)
    args.output.mkdir(parents=True,exist_ok=False)
    root=args.output.resolve(); prior=str(args.prior_root.resolve())+'/'
    for folder in ['data','temp','native']: (root/folder).mkdir()
    pd.Series(list(native)).to_csv(root/'data'/'Peaks.txt',sep='\t',index=False,header=False)
    pybedtools.set_tempdir(str(root/'temp'))
    # Native preprocessing hardcodes8workers; only scheduling is adapted.
    def serial_parallel(*a,**kw):
        kw['n_jobs']=1
        return joblib.Parallel(*a,**kw)
    pre.Parallel=serial_parallel; train.get_system_resources=lambda:(1,8.0)
    seed_audit=clone_linger_background_adaptation(train,meta,info,args.prior_root,args.shap_background=='training_cap')
    expression=rna.T.copy(); expression.index=selected
    openness=atac.T.copy(); openness.index=list(native)
    previous=Path.cwd(); os.chdir(root)
    try:
        out=str(root/'native')+'/'
        pre.preprocess(expression,openness,prior,'hg38','LINGER',out)
        # TFexp retains the native TF axes; restrict only target/index rows.
        names=pd.read_csv(out+'Symbol.txt',header=None,sep='\t')
        keep=np.flatnonzero(names[0].eq(target))
        if len(keep)!=1: raise ValueError('Target absent from native atlas')
        for name in ['Symbol.txt','Exp.txt','index.txt']:
            table=pd.read_csv(out+name,header=None,sep='\t',keep_default_na=False)
            table.iloc[keep].to_csv(out+name,header=False,index=False,sep='\t')
        train.training(prior,'LINGER',out,'ReLU','human')
        if seed_audit is not None:
            if len(seed_audit['backgrounds'])!=1 or len(seed_audit['prior_alignments'])!=1:
                raise ValueError('Missing or repeated native feature/background audit')
            (root/'seed_and_background_audit.json').write_text(json.dumps(seed_audit,indent=2)+'\n')
        links=LL_net.cis_shap(args.chromosome,out)
        if not set(links.RE).issubset(native) or not set(links.TG).issubset({target}):
            raise ValueError('Native result outside exact supplied features')
        links['peak_id']=links.RE.map(native); links['gene_id']=info['gene_id']
        links=links.sort_values(['score','peak_id'],ascending=[False,True])
        links['rank']=np.arange(len(links))+1
        links['method']='LINGER_native_bounded' if info['control']=='paired' else 'LINGER_native_bounded_pairing_shuffle'
        links['score_unit']='native_mean_absolute_SHAP'
        links['method_variant']='LINGER1.110_corrected_prior_RE_offset_and_audited_training_background'
        links['evidence_state']='prior_assisted_training_association_not_causal'
        for field,value in context.items(): links[field]=value
        links.to_csv(root/'links.tsv',sep='\t',index=False)
        (root/'links.training_donors.json').write_text(json.dumps({'donor_ids':sorted(set(meta.donor_id))})+'\n')
        info.update(context)
        info.update(status='native_training_fit_completed_not_evaluated',
          native_donor_weighting='unweighted_metacells',
          adaptations=['donor_confined_aggregation','bounded_peak_universe','single_target_response','serial_CPU','ambiguous_gene_symbols_excluded'],
          prior_training_context='unverified',prior_evaluation_reference_overlap='unverified',
          shap_background=args.shap_background,native_accuracy_placeholders='unusable_hardcoded0.5_values',
          prior_alignment='RE_prior_columns_offset_by_bulk_TF_count_local_repair',
          rna_only='unsupported_native_interface',prior_only='unsupported_native_interface')
        (root/'completion.json').write_text(json.dumps(info,indent=2)+'\n')
    finally: os.chdir(previous)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('linger-preflight'); q.add_argument('--prior-root',type=Path); q.add_argument('--chromosome',default='chr1'); q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('linger-model-preflight'); q.add_argument('--prior-root',type=Path,required=True)
    q.add_argument('--chromosome',required=True); q.add_argument('--output',type=Path,required=True)
    q.add_argument('--max-model-bytes',type=int,default=1024**3)
    q.add_argument('--target-symbols',default='')
    q=sub.add_parser('fit-linger-stage')
    for name in ['stage','genes','peaks','prior-root','model-preflight','output']: q.add_argument('--'+name,type=Path,required=True)
    q.add_argument('--chromosome',required=True)
    q.add_argument('--shap-background',choices=['native50','training_cap'],default='native50')
    q=sub.add_parser('prepare-tripod')
    for name in ['paired-root','genes','motif-map','output']: q.add_argument('--'+name,type=Path,required=True)
    q.add_argument('--targets',required=True); q.add_argument('--seed',type=int,default=42); q.add_argument('--cell-cap',type=int,default=1000)
    q=sub.add_parser('prepare-linger-stage')
    for name in ['paired-root','genes','output']: q.add_argument('--'+name,type=Path,required=True)
    q.add_argument('--target',required=True); q.add_argument('--stage',default='outer0'); q.add_argument('--seed',type=int,default=42)
    q=sub.add_parser('fit-tripod'); q.add_argument('--stage',type=Path,required=True); q.add_argument('--rscript',type=Path,required=True); q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('fit-tripod-batch'); q.add_argument('--stages',type=Path,required=True)
    q.add_argument('--rscript',type=Path,required=True); q.add_argument('--output',type=Path,required=True)
    q.add_argument('--control',choices=['paired','pairing_shuffle'],default='paired'); q.add_argument('--timeout',type=int,default=3600)
    a=p.parse_args()
    if a.command=='prepare-tripod': prepare(a)
    elif a.command=='prepare-linger-stage': prepare_linger_stage(a)
    elif a.command=='fit-tripod-batch': fit_batch(a)
    elif a.command=='linger-model-preflight': linger_model_preflight(a)
    elif a.command=='fit-linger-stage': fit_linger_stage(a)
    elif a.command=='fit-tripod':
        shared.require_compute()
        info=json.loads((a.stage/'stage.json').read_text())
        if info.get('input_units')!='TRIPOD_full_axis_CPM1e6_no_log':
            raise ValueError('Native TRIPOD requires full-axis CPM1e6 without log')
        meta=shared.read_table(a.stage/'metacells.tsv')
        if set(meta.donor_id)&set(info['excluded_donors']): raise ValueError('Donor leakage in stage')
        subprocess.run([str(a.rscript),str(Path(__file__).with_suffix('.R')),str(a.stage),str(a.output)],check=True,timeout=600)
        for name in ['links.tsv','trios.tsv']:
            table=pd.read_csv(a.output/name,sep='\t')
            table['outer_fold']=int(info['stage'].split('_')[0].replace('outer',''))
            table['training_stage']=info['stage']; table['control']=info['control']
            for field in ['cohort','lineage','lineage_state','source_lineage']:
                if field in meta: table[field]=';'.join(sorted(set(meta[field])))
            table.to_csv(a.output/name,sep='\t',index=False)
        info['status']='native_training_fit_completed_not_evaluated'
        (a.output/'completion.json').write_text(json.dumps(info,indent=2)+'\n')
        (a.output/'links.training_donors.json').write_text(json.dumps(
            {'donor_ids':sorted(set(meta.donor_id))},indent=2)+'\n')
    else:
        with a.output.open('x') as h: json.dump(linger_inventory(a.prior_root,a.chromosome),h,indent=2)


if __name__=='__main__': main()
