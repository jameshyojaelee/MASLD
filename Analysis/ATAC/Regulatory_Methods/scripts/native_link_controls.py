#!/usr/bin/env python3
"""RNA-only and prior-only link-learning REFERENCE controls, not native ablations.

These proxy scores differ from native TRIPOD coefficients and LINGER SHAP.
They are neither probabilities nor TF activity. Held-donor prediction and modality
controls are still required; reference-control comparisons alone cannot prove
multimodal benefit. Inputs must be corrected, complete, immutable training stages.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def read(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def digest(path):
    sha=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for piece in iter(lambda:handle.read(1024*1024),b''):sha.update(piece)
    return sha.hexdigest()


def require_valid(root):
    root=Path(root)
    for path in root.iterdir():
        if any(word in path.name.lower() for word in ['invalidat','obsolete','superseded']):
            raise ValueError('Invalidated source: '+str(path))
    marker=root/'completion.json'
    if not marker.exists() or json.loads(marker.read_text()).get('status')!='completed':
        raise ValueError('Incomplete source: '+str(root))


def weighted_correlations(target, expression, donors):
    """Equal total weight per biological donor; constant columns return zero."""
    y=np.asarray(target,dtype=float);x=np.asarray(expression,dtype=float)
    if x.ndim==1:x=x[:,None]
    donors=np.asarray(donors,dtype=str)
    if len(y)!=len(donors) or len(x)!=len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Invalid RNA/donor axes')
    _,inverse,counts=np.unique(donors,return_inverse=True,return_counts=True)
    weights=1./counts[inverse];weights/=weights.sum()
    yc=y-np.sum(weights*y);xc=x-np.sum(weights[:,None]*x,axis=0)
    numerator=np.sum(weights[:,None]*yc[:,None]*xc,axis=0)
    denominator=np.sqrt(np.sum(weights*yc*yc)*np.sum(weights[:,None]*xc*xc,axis=0))
    return np.divide(numerator,denominator,out=np.zeros_like(numerator),where=denominator>0)


def tripod_controls(rna, motifs, peaks, tss0, target, donors):
    # Dedup shared/repeated motif annotations; self-TF cannot predict itself.
    eligible=motifs[motifs.tf_gene_id.ne(target)&motifs.tf_gene_id.isin(rna.columns)].drop_duplicates(['peak_id','tf_gene_id'])
    tf_names=sorted(set(eligible.tf_gene_id))
    correlations=dict(zip(tf_names,np.abs(weighted_correlations(rna[target],rna[tf_names],donors))))
    geometry=peaks.set_index('peak_id')
    rows=[]
    for peak,group in eligible.groupby('peak_id',sort=True):
        if peak not in geometry.index:raise ValueError('Motif outside exact peak geometry')
        p=geometry.loc[peak]
        distance=abs((float(p.start0)+float(p.end))/2-float(tss0))
        rows.append({'gene_id':target,'peak_id':peak,'RNA_only':max(correlations[t] for t in set(group.tf_gene_id)), 'prior_only':-distance,'n_mapped_TFs':len(set(group.tf_gene_id))})
    return pd.DataFrame(rows,columns=['gene_id','peak_id','RNA_only','prior_only','n_mapped_TFs'])


def linger_controls(target, target_symbol, expression, tf_expression, tf_names, binding, donors, identity, first_layer):
    """Explicit verified feature map prevents observed-TF-count offset corruption."""
    needed={'peak_id','RE_input_position0','prior_column0','expected_identity','actual_identity','match','bulk_TFs'}
    if not needed.issubset(identity):raise ValueError('Missing verified LINGER feature identity fields')
    if not identity.match.astype(str).str.lower().eq('true').all() or not identity.expected_identity.eq(identity.actual_identity).all():
        raise ValueError('Unverified prior identity')
    if not identity.bulk_TFs.astype(int).eq(638).all():raise ValueError('DPM1 bulkTF axis must638, never observed486')
    if first_layer.shape[1]!=1180:raise ValueError('Unexpected DPM1 bulkNN first-layer axis')
    if not np.isfinite(first_layer).all():raise ValueError('Nonfinite static first-layer weights')
    columns=identity.prior_column0.astype(int).to_numpy()
    if (columns<638).any() or (columns>=first_layer.shape[1]).any():raise ValueError('RE column outside correctly offset prior axis')
    tf_names=np.asarray(tf_names,dtype=str)
    correlations=np.abs(weighted_correlations(expression,tf_expression,donors))
    correlations[tf_names==target_symbol]=0
    binding=np.asarray(binding,dtype=float)
    if binding.shape!=(len(tf_names),int(identity.RE_input_position0.astype(int).max())+1):
        raise ValueError('TF-binding axes disagree with explicit RE input identities')
    if not np.isfinite(binding).all() or (binding<0).any():raise ValueError('Expected nonnegative static TF binding')
    rows=[]
    for peak,group in identity.groupby('peak_id',sort=True):
        positions=group.RE_input_position0.astype(int).unique()
        priors=group.prior_column0.astype(int).unique()
        rna_score=float(np.max(correlations[:,None]*binding[:,positions]))
        prior_score=float(np.max(np.abs(first_layer[:,priors]).sum(axis=0)))
        rows.append({'gene_id':target,'peak_id':peak,'RNA_only':rna_score,'prior_only':prior_score,'n_mapped_prior_RE_columns':len(priors)})
    return pd.DataFrame(rows)


def publish(output, results, donors, sources, family, metadata, raw_inputs):
    # Publication requires exact20 stage families, after source hashes rechecked.
    expected={str(f)+':'+s for f in range(5) for s in ['outer','inner0','inner1','inner2']}
    if set(results)!=expected:raise ValueError('Require exactly20 complete outer/inner stage maps')
    for path,sha in sources.items():
        if digest(path)!=sha:raise ValueError('Source changed during control computation: '+path)
    for path,sha in raw_inputs.items():
        if digest(path)!=sha:raise ValueError('Raw input changed before publication: '+path)
    config={'external_stage_links':{},'external_stage_roots':{}}
    for control in ['RNA_only','prior_only']:
        method=family+'_reference_'+control
        config['external_stage_links'][method]={};config['external_stage_roots'][method]=str(output.resolve())
        for key,table in sorted(results.items()):
            frame=table.copy();frame['score']=frame[control];frame=frame.sort_values(['gene_id','score','peak_id'],ascending=[True,False,True])
            frame['rank']=frame.groupby('gene_id').cumcount()+1
            frame['method']=method;frame['score_unit']='max_absolute_donor_weighted_RNA_correlation' if control=='RNA_only' and family=='TRIPOD' else 'negative_TSS_distance_bp' if control=='prior_only' and family=='TRIPOD' else 'max_absolute_RNA_correlation_times_static_binding' if control=='RNA_only' else 'maximum_L1_norm_of_bulkNN_first_layer_RE_column'
            frame['method_variant']='reference_link_learning_proxy_not_native_ablation'
            frame['evidence_state']='reference_score_not_probability_TF_activity_or_causality'
            frame['outer_fold']=int(key.split(':')[0]);frame['training_stage']=key
            folder=output/'combined'/key.replace(':','_')/method;folder.mkdir(parents=True)
            frame.to_csv(folder/'links.tsv',sep='\t',index=False)
            (folder/'links.training_donors.json').write_text(json.dumps({'donor_ids':sorted(donors[key])})+'\n')
            (folder/'completion.json').write_text(json.dumps({'status':'completed','scope':'training_reference_links_only'})+'\n')
            config['external_stage_links'][method][key]=str((folder/'links.tsv').resolve())
    (output/'input_sha256.json').write_text(json.dumps(raw_inputs,indent=2)+'\n')
    (output/'external_stage_config.json').write_text(json.dumps(config,indent=2)+'\n')
    (output/'controls_source_snapshot.py').write_bytes(Path(__file__).read_bytes())
    (output/'source_sha256.json').write_text(json.dumps(sources,indent=2)+'\n')
    (output/'completion.json').write_text(json.dumps({'status':'completed','scope':'reference_link_learning_controls_not_native_ablations','native_ablation_claim':False,'multimodal_benefit_claim':False,'scientific_support_calls':False,**metadata},indent=2)+'\n')


def verify_tripod_derivation(args):
    """Reconstruct all stage values with full-assay denominators, one assay at once."""
    import gc
    import importlib.util
    from scipy import sparse
    spec=importlib.util.spec_from_file_location('frozen_control_aggregation',args.core_source)
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    cells=read(args.paired_root/'cells.tsv')
    _,meta,groups,_=core.build_metacells(sparse.csr_matrix((len(cells),1)),cells,42,50,100)
    if len(meta)!=78 or meta.donor_id.nunique()!=39:raise ValueError('Expected full39-donor balanced78-aggregate source')
    records=[r for r in json.loads((args.stages/'stages.json').read_text()) if r.get('control')=='paired']
    rna_ids=set();peak_ids=set()
    for r in records:
        path=args.stages/r['gene_id']/r['stage']/'paired'
        rna_ids.update(pd.read_csv(path/'rna.tsv',sep='\t',nrows=0).columns)
        peak_ids.update(pd.read_csv(path/'atac.tsv',sep='\t',nrows=0).columns)
    projections={};denominators={}
    for assay,ids,names,feature_table,feature_column in [('ATAC',sorted(peak_ids),['cell_counts.npz','ATAC_counts.mtx.gz'],'peaks.tsv','peak_id'),('RNA',sorted(rna_ids),['rna_counts.npz'],'genes.tsv','gene_id')]:
        features=read(args.paired_root/feature_table)
        positions=features.set_index(feature_column).assign(position=np.arange(len(features))).loc[ids,'position'].to_numpy()
        matrix=core.read_matrix(args.paired_root,names)
        if matrix.shape!=(len(cells),len(features)):raise ValueError('Raw source axes differ')
        totals=np.asarray(matrix.sum(axis=1)).ravel()
        small=matrix[:,positions].tocsr();del matrix;gc.collect()
        aggregate=np.vstack([np.asarray(small[group].sum(axis=0)).ravel() for group in groups])
        libraries=np.array([totals[group].sum() for group in groups])
        if (libraries<=0).any():raise ValueError('Zero source assay depth')
        projections[assay]=pd.DataFrame(aggregate*(1e6/libraries[:,None]),columns=ids,index=meta.metacell_id)
        denominators[assay]=libraries.tolist();del small;gc.collect()
    expected={};outer=core.donor_folds(meta.donor_id.to_numpy(),5,42)
    for fold in range(5):
        train=np.flatnonzero(outer!=fold);expected[f'outer{fold}']=train
        inner=core.donor_folds(meta.donor_id.to_numpy()[train],3,core.stable_seed(42,('inner',fold)))
        for j in range(3):expected[f'outer{fold}_inner{j}']=train[inner!=j]
    checks=[]
    for r in records:
        path=args.stages/r['gene_id']/r['stage']/'paired';actual_meta=read(path/'metacells.tsv');wanted_meta=meta.iloc[expected[r['stage']]]
        fields=['cohort','lineage','donor_id','metacell_id','n_cells']
        if not actual_meta[fields].astype(str).reset_index(drop=True).equals(wanted_meta[fields].astype(str).reset_index(drop=True)):
            raise ValueError('Prepared stage differs from reconstructed aggregate/fold membership')
        errors={}
        for assay,name in [('RNA','rna.tsv'),('ATAC','atac.tsv')]:
            emitted=pd.read_csv(path/name,sep='\t');want=projections[assay].loc[wanted_meta.metacell_id,emitted.columns].to_numpy()
            if not np.allclose(emitted.to_numpy(),want,rtol=1e-12,atol=1e-9):raise ValueError('Prepared values do not reconstruct from raw '+assay)
            errors[assay]=float(np.max(np.abs(emitted.to_numpy()-want)))
        emitted=pd.read_csv(path.parent/'pairing_shuffle'/'rna.tsv',sep='\t')
        want=projections['RNA'].loc[wanted_meta.metacell_id,emitted.columns].to_numpy().copy();rng=np.random.default_rng(core.stable_seed(42,r['gene_id']+r['stage']))
        for donor in sorted(set(wanted_meta.donor_id)):
            rows=np.flatnonzero(wanted_meta.donor_id.to_numpy()==donor);want[rows]=want[rng.permutation(rows)]
        if not np.allclose(emitted.to_numpy(),want,rtol=1e-12,atol=1e-9):raise ValueError('Prepared RNA shuffle does not reconstruct')
        checks.append({'gene_id':r['gene_id'],'stage':r['stage'],'RNA_max_absolute_error':errors['RNA'],'ATAC_max_absolute_error':errors['ATAC'],'all_training_memberships_and_shuffle_values_verified':True})
    membership=[{'metacell_id':meta.iloc[i].metacell_id,'cell_id':cells.iloc[j].cell_id,'donor_id':meta.iloc[i].donor_id} for i,group in enumerate(groups) for j in group]
    pd.DataFrame(membership).to_csv(args.output/'reconstructed_aggregation_membership.tsv',sep='\t',index=False)
    receipt={'status':'verified_direct_raw_reconstruction','frozen_core':str(args.core_source.resolve()),'frozen_core_sha256':digest(args.core_source),'seed':42,'cell_cap':100,'aggregate_size':50,'n_biological_donors':39,'n_aggregates':78,'library_denominators_full_assay_axis':denominators,'n_gene_stage_checks':len(checks),'n_distinct_training_stages':len(expected),'checks':checks,'tolerance':{'rtol':1e-12,'atol':1e-9},'historical_binding':'derivation_directly_reconstructed_from_verified_raw_files_not_assumed_from_manifest_copy'}
    (args.output/'raw_stage_derivation.json').write_text(json.dumps(receipt,indent=2)+'\n')


def tripod(args):
    require_valid(args.native_root)
    verify_tripod_derivation(args)
    sources={str(p.resolve()):digest(p) for p in [Path(__file__),args.raw_input_manifest,args.core_source,args.output/'raw_stage_derivation.json',args.stages/'stages.json',args.native_root/'completion.json',args.genes,args.peaks]}
    prepared_sources=json.loads((args.stages/'source_sha256.json').read_text())
    for path,sha in prepared_sources.items():
        if digest(path)!=sha:raise ValueError('Prepared source identity changed')
        if path in args.raw_inputs and args.raw_inputs[path]!=sha:raise ValueError('Raw/prepared source identity disagreement')
    sources.update(prepared_sources)
    records=json.loads((args.stages/'stages.json').read_text());records=[r for r in records if r.get('control')=='paired']
    annotation=read(args.genes).set_index('gene_id');peaks=read(args.peaks)
    results={};donors={}
    for record in records:
        folder=args.stages/record['gene_id']/record['stage']/'paired'
        for path in folder.iterdir():
            if path.is_file():sources[str(path.resolve())]=digest(path)
        info=json.loads((folder/'stage.json').read_text())
        if info['input_units']!='TRIPOD_full_axis_CPM1e6_no_log':raise ValueError('Wrong TRIPOD RNA units')
        key=record['stage'].replace('outer','',1).replace('_',':')
        if ':' not in key:key+=':outer'
        native_stage=args.native_root/'combined'/record['stage'];require_valid(native_stage)
        for name in ['completion.json','links.training_donors.json']:
            path=native_stage/name;sources[str(path.resolve())]=digest(path)
        meta=read(folder/'metacells.tsv');actual=set(meta.donor_id)
        if actual!=set(json.loads((native_stage/'links.training_donors.json').read_text())['donor_ids']):raise ValueError('Completed native source training donor mismatch')
        if actual&set(info['excluded_donors']):raise ValueError('Held-out donor leakage')
        if key in donors and donors[key]!=actual:raise ValueError('Gene-specific donor coverage mismatch')
        donors[key]=actual
        rna=pd.read_csv(folder/'rna.tsv',sep='\t');motifs=read(folder/'motifs.tsv')
        if len(rna)!=len(meta):raise ValueError('RNA/metacell axis mismatch')
        table=tripod_controls(rna,motifs,peaks,float(annotation.loc[record['gene_id'],'tss0']),record['gene_id'],meta.donor_id)
        shuffle_path=folder.parent/'pairing_shuffle'/'rna.tsv'
        sources[str(shuffle_path.resolve())]=digest(shuffle_path)
        shuffled=tripod_controls(pd.read_csv(shuffle_path,sep='\t'),motifs,peaks,float(annotation.loc[record['gene_id'],'tss0']),record['gene_id'],meta.donor_id)
        if not table.peak_id.equals(shuffled.peak_id) or not np.allclose(table[['RNA_only','prior_only']],shuffled[['RNA_only','prior_only']],rtol=1e-12,atol=1e-12):
            raise ValueError('Reference changed after within-donor whole-RNA shuffle')
        table['cohort']=info['cohort'];table['lineage']=info['lineage'];table['lineage_state']=info['lineage_state']
        results.setdefault(key,[]).append(table)
    publish(args.output,{k:pd.concat(v,ignore_index=True) for k,v in results.items()},donors,sources,'TRIPOD',{'input_units':'TRIPOD_full_axis_CPM1e6_no_log','weighting':'equal_total_weight_per_donor','controls_do_not_use_current_ATAC_values':True,'actual_whole_RNA_pairing_shuffle_invariance_checked':True,'native_objective_difference':'motif/distance or RNA correlation proxy, not native interaction coefficient'},args.raw_inputs)


def parse_ids(value):
    return [int(x) for x in str(value).split('_')]


def linger_identity(fit, stage, prior, peaks, global_re_names):
    audit=json.loads((fit/'seed_and_background_audit.json').read_text())
    alignment=audit['prior_alignments']
    if len(alignment)!=1 or alignment[0]['bulk_TFs']!=638 or alignment[0]['bulk_REs']!=542:
        raise ValueError('Unqualified bulk DPM1 prior alignment')
    alignment=alignment[0]
    index=pd.read_csv(fit/'native/index.txt',sep='\t',header=None,keep_default_na=False)
    if len(index)!=1 or index.iloc[0,0]!='DPM1':raise ValueError('Expected target-only DPM1 index')
    re_indices=parse_ids(index.iloc[0,1]);tf_indices=parse_ids(index.iloc[0,2]);re_global=[x+1 for x in parse_ids(index.iloc[0,3])]
    tf_names=pd.read_csv(fit/'native/TFName.txt',header=None,sep='\t')[0].astype(str).to_numpy()
    prior_tf=pd.read_csv(prior/'TFName.txt',header=None,sep='\t')[0].astype(str).to_numpy()
    physical=pd.read_csv(stage/'atac.tsv',sep='\t',nrows=0).columns.tolist()
    geometry=peaks.set_index('peak_id').loc[physical]
    coordinates=(geometry.chrom+':'+geometry.start0.astype(str)+'-'+geometry.end.astype(str)).tolist()
    declared=pd.read_csv(fit/'data/Peaks.txt',header=None,sep='\t')[0].tolist()
    if coordinates!=declared:raise ValueError('Physical peak BED identity/order mismatch')
    bulk_index=pd.read_csv(prior/'chr20_index.txt',header=None,sep='\t',keep_default_na=False).iloc[134]
    bulk_re=parse_ids(pd.read_csv(prior/'chr20_index_all.txt',header=None,sep='\t',keep_default_na=False).iloc[134,0])
    bulk_tf=parse_ids(bulk_index[2])
    if bulk_index[0]!='DPM1' or len(bulk_tf)!=638 or len(bulk_re)!=542:raise ValueError('Actual bulk feature axis mismatch')
    needed=set(re_global)-set(global_re_names)
    if needed:
        with (prior/'all_hg19.txt').open() as handle:
            for i,line in enumerate(handle,1):
                if i in needed:global_re_names[i]=line.rstrip('\n\r')
                if i>=max(needed):break
    features=alignment['features']
    if len(features)!=len(tf_indices)+len(re_indices):raise ValueError('Observed feature axis mismatch')
    records=[]
    for position,tf_index in enumerate(tf_indices):
        item=features[position];global_id=int(item['global_prior_id1']);column=int(item['prior_column0'])
        if item['feature_type']!='TF' or not item['match'] or item['input_position0']!=position or tf_names[tf_index]!=item['expected_identity'] or prior_tf[global_id-1]!=item['actual_identity'] or bulk_tf[column]!=global_id:
            raise ValueError('Observed TF identity does not match actual prior axis')
    for j,(physical_index,global_id) in enumerate(zip(re_indices,re_global)):
        position=len(tf_indices)+j;item=features[position];column=int(item['prior_column0'])
        if item['feature_type']!='RE' or not item['match'] or item['input_position0']!=position or int(item['global_prior_id1'])!=global_id or global_re_names[global_id]!=item['expected_identity'] or item['expected_identity']!=item['actual_identity'] or column<638 or bulk_re[column-638]!=global_id:
            raise ValueError('Observed RE identity does not match bulkTF638-offset prior axis')
        records.append(dict(item,peak_id=physical[physical_index],RE_input_position0=j,physical_peak_index0=physical_index,bulk_TFs=638))
    return pd.DataFrame(records),tf_indices,re_indices,tf_names


def linger(args):
    require_valid(args.native_root)
    qualification=json.loads(args.model_preflight.read_text())
    checked=[x for x in qualification['checked_gene_rows'] if x['symbol']=='DPM1']
    if len(checked)!=1 or checked[0]['gene_row0']!=134 or checked[0]['input_size']!=1180:
        raise ValueError('Missing actual DPM1 model preflight')
    model=args.prior_root/'all_models_chr20.pt'
    if digest(model)!=qualification['models']['all_models_']['sha256']:raise ValueError('Actual model SHA disagrees with preflight')
    import torch
    torch.set_num_threads(1)
    states=torch.load(model,map_location='cpu',weights_only=False)
    first_layer=states[134]['fc1.weight'].detach().cpu().numpy().copy();del states
    peaks=read(args.peaks)
    tracked=[Path(__file__),args.raw_input_manifest,args.model_preflight,args.peaks,model,args.native_root/'completion.json']
    sources={str(p.resolve()):digest(p) for p in tracked}
    for name in ['TFName.txt','all_hg19.txt','chr20_index.txt','chr20_index_all.txt']:
        path=args.prior_root/name;sources[str(path.resolve())]=digest(path)
    results={};donors={};global_re_names={};identities=args.output/'verified_feature_identity';identities.mkdir()
    for outer in range(5):
        for suffix in ['', '_inner0','_inner1','_inner2']:
            stage_name=f'outer{outer}{suffix}';key=str(outer)+(':'+suffix[1:] if suffix else ':outer')
            fit=args.native_root/'fits'/stage_name/'paired';stage=args.native_root/'inputs'/stage_name/'paired'
            stage_marker=args.native_root/'stages'/stage_name/'paired';require_valid(stage_marker)
            info=json.loads((stage/'stage.json').read_text());meta=read(stage/'metacells.tsv')
            if info['input_units']!='LINGER_per_nucleus_transform_then_mean':raise ValueError('Incorrect LINGER source units')
            if set(meta.donor_id)&set(info['excluded_donors']):raise ValueError('Heldout donor leakage')
            binding_receipt=json.loads((fit/'source_and_input_sha256.json').read_text())
            for name,sha in binding_receipt['inputs'].items():
                if digest(stage/name)!=sha:raise ValueError('Fit does not bind actual immutable stage inputs')
            identity,tf_indices,re_indices,tf_names=linger_identity(fit,stage,args.prior_root,peaks,global_re_names)
            units=json.loads((fit/'native_units_check.json').read_text())
            fit_marker=json.loads((fit/'completion.json').read_text())
            if units.get('status')!='PASS' or set(units['donor_ids'])!=set(meta.donor_id) or units['training_aggregates']!=len(meta) or units['native_TFs']!=len(tf_indices):raise ValueError('Unqualified native unit/donor check')
            if fit_marker.get('input_units')!='LINGER_per_nucleus_transform_then_mean' or fit_marker.get('prior_alignment')!='RE_prior_columns_offset_by_bulk_TF_count_local_repair':raise ValueError('Wrong native fit units/prior repair variant')
            identity.to_csv(identities/(stage_name+'.tsv'),sep='\t',index=False)
            for folder in [stage,fit/'native']:
                for path in folder.iterdir():
                    if path.is_file() and path.suffix in ['.tsv','.txt','.json','.bed']:sources[str(path.resolve())]=digest(path)
            for name in ['source_and_input_sha256.json','seed_and_background_audit.json','completion.json','native_units_check.json']:
                path=fit/name;sources[str(path.resolve())]=digest(path)
            expression=pd.read_csv(fit/'native/Exp.txt',sep='\t',header=None).iloc[0].to_numpy()
            tf=pd.read_csv(fit/'native/TFexp.txt',sep='\t',header=None).to_numpy()[tf_indices].T
            binding=pd.read_csv(fit/'native/TF_binding.txt',sep='\t',header=None).to_numpy()[re_indices][:,tf_indices].T
            table=linger_controls(info['gene_id'],'DPM1',expression,tf,tf_names[tf_indices],binding,meta.donor_id,identity,first_layer)
            for field in ['cohort','lineage','lineage_state']:table[field]=info[field]
            results[key]=table;donors[key]=set(meta.donor_id)
    publish(args.output,results,donors,sources,'LINGER',{'input_units':'LINGER_per_nucleus_transform_then_mean_then_native_preprocessing','weighting':'equal_total_weight_per_donor','bulk_TFs':638,'bulk_model_gene_row0':134,'prior_only_uses_current_RNA_or_ATAC_values':False,'static_prior_conditioning':'fixed_assay_observability_and_verified_prior_feature_identities','native_objective_difference':'RNA-correlation/binding or bulk first-layer L1 proxy, not native SHAP','prior_evaluation_reference_overlap':'unverified'},args.raw_inputs)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('tripod')
    for name in ['stages','native-root','raw-input-manifest','paired-root','core-source','genes','peaks','output']:p.add_argument('--'+name,required=True,type=Path)
    p=sub.add_parser('linger')
    for name in ['native-root','raw-input-manifest','prior-root','model-preflight','peaks','output']:p.add_argument('--'+name,required=True,type=Path)
    args=parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Run in existing compute allocation')
    if args.output.exists():raise FileExistsError(args.output)
    require_valid(args.native_root)
    args.raw_inputs=json.loads(args.raw_input_manifest.read_text())
    if not args.raw_inputs:raise ValueError('Empty raw source manifest')
    for path,sha in args.raw_inputs.items():
        if str(Path(path).resolve())!=path or digest(path)!=sha:raise ValueError('Unverified raw native source input: '+path)
    if str(args.peaks.resolve()) not in args.raw_inputs:raise ValueError('Raw source manifest does not cover supplied peaks')
    args.output.mkdir(parents=True)
    (args.output/'raw_input_binding.json').write_text(json.dumps({'source_manifest':str(args.raw_input_manifest.resolve()),'source_manifest_sha256':digest(args.raw_input_manifest),'all_hashes_verified_before_computation':True,'n_raw_files':len(args.raw_inputs)},indent=2)+'\n')
    (tripod if args.command=='tripod' else linger)(args)


if __name__=='__main__':main()
