"""Registered seventeen-donor native ChIP transport comparison; never refit.

Requires successful H3 job21998882 and its RNA hashes bound before H3 reads.
Source-person disjointness remains unverified. Counts define qualified
single-read ChIP labels, not reconstruction of the source CUT&RUN assay.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback
import types

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
H3 = REC / 'codex_receiving_h3_reads_21998882'
SOURCES = {
    'export_training_h3_target_transform':'3499fb927373a04c7fc7121bd89598709774deb38f652b26d5506275e09ee907',
    'profile_error_decomposition':'f75a8ebacd9383962d11d16849efb153b621e330f6947b9abc22876fe2601d29',
    'codex_paired_profile_evaluation':'74eabdea369e15159c58651cad5999b4417cd5881231d203a3b988915328ae70',
    'apply_fixed_h3_residual_target':'0f508019599333f02f0c7164da00dd4bd99610f3c84f03366ac2492d2f3f05f3',
    'assemble_fixed_h3_counts':'0fb8317a9ffa4ccebd3e621915d9efbe72bdf93a8753fd5c4e983631b114d81f',
    'codex_receiving_h3_reads':'f641ea6485515b7f986f3ce279f449c5e6b85e99f676c3a5ac959fc93bf8b509',
}
CHECKS = {
    REC/'fixed_h3_assembly_checks_21998841/summary.json':'473234d4e51bfc2502e5618397096be1c74ec36ba962896a6215251168d7e65e',
    REC/'fixed_h3_target_checks_21998710/summary.json':'810cbf988603330c1b2a64bb5c59d64c51a629327b50cee80513a1d621ebeb86',
    REC/'codex_paired_profile_checks_21998524/summary.json':'71046ad2b7261aa78f573eb2efd9df52f6200fe997ec7fd202055bf4b1556b25',
    REC/'codex_paired_profile_checks_21998526/summary.json':'4a3e2854c5ccb296f8dd19ae681e5606f4927460fbd2c840b306efe64b19ceda',
}
OUTPUT_CAP = 200_000_000


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def guard(paths):
    for p,h in paths.items():
        require(sha(p)==h,'Frozen input changed: '+str(p))


def read_json(path):
    require(path.is_file() and path.stat().st_size<=10_000_000,'Missing/bounded receipt: '+str(path))
    return json.loads(path.read_text())


def write_json(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2);f.write('\n')


def load_sources(output):
    """Execute guarded snapshots while retaining original fixed-root __file__."""
    archive=output/'executed_sources';archive.mkdir()
    for name,expected in SOURCES.items():
        original=HERE/(name+'.py');data=original.read_bytes()
        require(hashlib.sha256(data).hexdigest()==expected,'Source snapshot changed: '+name)
        (archive/original.name).write_bytes(data)
        require(name not in sys.modules,'Numerical helper already imported: '+name)
        module=types.ModuleType(name);module.__file__=str(original);sys.modules[name]=module
        exec(compile(data,str(original),'exec'),module.__dict__)
    for p in [Path(__file__),HERE/'run_codex_receiving_profile_evaluation.sbatch']:
        (archive/p.name).write_bytes(p.read_bytes())
    return {name:sys.modules[name] for name in SOURCES}


def admitted_tables(modules):
    """Metadata/hash admission precedes count or prediction deserialization."""
    helper=modules['codex_receiving_h3_reads'];assembler=modules['assemble_fixed_h3_counts']
    complete=read_json(H3/'summary.json');protocol=read_json(H3/'protocol.json')
    require(complete.get('all17_reads_processed') and tuple(complete['completed'])==helper.DONORS
            and complete.get('biological_n')==17 and not complete.get('predictive_accuracy_evaluated')
            and set(complete['donor_summary_sha256'])==set(helper.DONORS),'H3 cohort not complete/exact')
    require(sha(H3/'protocol.json')==complete['protocol_sha256']
            and tuple(protocol['fixed_donors'])==helper.DONORS
            and protocol['h3_only_files']==202 and protocol['compressed_bytes']==56783483922
            and protocol['RNA_predictions_frozen_before_H3'] and not protocol['model_fit']
            and not protocol['predictive_accuracy_evaluated'],'H3 protocol/measurement boundary differs')
    # Requalify the exact B1/four-batch native RNA receipts and paired hashes.
    prerequisite,_,_=helper.prerequisites()
    bound={Path(p):h for p,h in protocol['input_sha256'].items()}
    require(all(bound.get(p)==h for p,h in prerequisite.items()),
            'RNA/index prerequisite was not hash-bound in H3 protocol before H3 processing')
    require(bound.get(HERE/'codex_receiving_h3_reads.py')==SOURCES['codex_receiving_h3_reads'],
            'H3 executed helper lineage differs')
    guard(bound)
    frozen=dict(bound);frozen.update(CHECKS)
    frozen[H3/'summary.json']=sha(H3/'summary.json');frozen[H3/'protocol.json']=complete['protocol_sha256']
    tables=[]
    for donor in helper.DONORS:
        base=H3/donor;p=base/'summary.json';s=read_json(p)
        require(sha(p)==complete['donor_summary_sha256'][donor] and s['donor']==donor and s['h3_reads_processed']
                and not s['predictions_refitted'] and not s['count_axis_assembly_or_target_evaluation_run'],
                'H3 donor completion/identity differs')
        count=Path(s['count_table'])
        require(count==base/'native_featurecounts.tsv','Count path escaped fixed H3 donor')
        frozen[p]=complete['donor_summary_sha256'][donor]
        frozen[count]=s['count_sha256'];frozen[Path(str(count)+'.summary')]=s['count_summary_sha256']
        require(all(Path(column).parent==base for column in s['columns']) and
                all(binding['donor']==donor for binding in s['columns'].values()),'BAM/run donor ownership differs')
        tables.append(assembler.CountTable(path=count,columns=s['columns'],count_command=tuple(s['count_command']),
                                          prefilter_commands=s['prefilter_commands']))
    frozen.update({HERE/(name+'.py'):h for name,h in SOURCES.items()})
    frozen[Path(__file__)]=sha(Path(__file__))
    frozen[HERE/'run_codex_receiving_profile_evaluation.sbatch']=sha(HERE/'run_codex_receiving_profile_evaluation.sbatch')
    guard(frozen)
    return tables,frozen,protocol


def prediction_paths():
    directories={'B1':REC/'codex_native_rna_pilot_21998663'}
    for job,group in zip((21998814,21998815,21998816,21998817),
                         (('B7','B8','B9','B10'),('B12','B15','B19','B21'),
                          ('B22','B24','B26','B36'),('B38','B41','B46','B47'))):
        directories.update({d:REC/('codex_native_rna_batch_'+str(job))/d for d in group})
    return directories


def evaluate_completed_transport(modules,output):
    """Callable for reviewed compute execution; no retraining or target tuning."""
    require(os.environ.get('SLURM_JOB_ID'),'Receiving evaluation only under SLURM')
    import numpy as np
    tables,frozen,h3_protocol=admitted_tables(modules)
    assembler=modules['assemble_fixed_h3_counts'];target_helper=modules['apply_fixed_h3_residual_target']
    evaluator=modules['codex_paired_profile_evaluation']
    require(evaluator.SEED==20260930 and evaluator.BOOTSTRAPS==10000 and evaluator.REGIONS==96460,
            'Registered inference protocol differs')
    guard(target_helper.GUARDS);frozen.update(target_helper.GUARDS)
    write_json(output/'protocol.json',dict(input_sha256={str(p):h for p,h in frozen.items()},
        biological_n=17,regions=96460,donors=evaluator.DONORS,bootstrap_draws=10000,seed=20260930,
        primary='mean donor MSE(global rrr) minus mean donor MSE(local rrr_offset_cis)',
        hypothesis_family='one registered paired mean-error contrast; no regional tests or p-values',
        target='log2 modeled-region CPM minus fixed [1, native concentration descriptors] @ B99',
        baseline='unchanged fixed training residual mean; zero secondary; no receiver centering or calibration',
        qualification='approximate native RNA quantification and single-read ChIP transport of CUT&RUN-trained frozen models',
        source_person_disjointness='unverified',bootstrap_limits='conditional on frozen models, B99, reference/counting and measured labels; no refitting, technical error or assay-shift uncertainty',
        H3_protocol_sha256=sha(H3/'protocol.json'),models_refitted=False,receiving_target_selected=False,
        donors_filtered_by_diagnostics=False))
    assembled=assembler.assemble_fixed_h3_counts(tables)
    require(tuple(assembled['donor_ids'])==evaluator.DONORS,'Assembler changed fixed donor order')
    definition=target_helper.load_fixed_definition()
    target,diagnostics=target_helper.apply_fixed_target(assembled['counts'],assembled['donor_ids'],assembled['region_keys'])
    prediction={'global':[],'local':[]};regions=assembled['region_keys'];directories=prediction_paths()
    for donor in evaluator.DONORS:
        for arm,form in [('global','rrr'),('local','rrr_offset_cis')]:
            p=directories[donor]/(arm+'_prediction.npz')
            require(p in frozen and sha(p)==frozen[p],'Prediction not frozen before receiving H3')
            # Fixed successful producer stores string axes without object pickle.
            with np.load(p,allow_pickle=False) as archive:
                require(set(archive.files)=={'profile','region_key','sample_id','form','transport','training_mean'},
                        'Prediction archive schema differs')
                require(tuple(archive['sample_id'].astype(str))==(donor,)
                        and tuple(archive['region_key'].astype(str))==tuple(regions)
                        and archive['form'].item()==form and archive['transport'].item()=='raw',
                        'Prediction native donor/region/form/transport differs')
                require(np.array_equal(archive['training_mean'],definition['training_mean']),
                        'Shipped fixed training-mean baseline differs from frozen export')
                values=archive['profile'].copy()
                require(values.shape==(1,96460) and np.isfinite(values).all(),'Nonfinite/partial prediction')
                prediction[arm].append(values[0])
    global_values=np.stack(prediction['global']);local_values=np.stack(prediction['local'])
    result=evaluator.evaluate_fixed_arrays(target,global_values,local_values,definition['training_mean'],
        evaluator.DONORS,evaluator.DONORS,evaluator.DONORS,regions,regions,regions,regions)
    guard(frozen)
    np.save(output/'complete_native_single_read_counts.npy',assembled['counts'],allow_pickle=False)
    np.save(output/'fixed_B99_residual_target.npy',target,allow_pickle=False)
    write_json(output/'count_assembly_evidence.json',assembled['evidence'])
    write_json(output/'descriptor_diagnostics.json',diagnostics)
    write_json(output/'paired_profile_results.json',result)
    output_bytes=sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    require(output_bytes<=OUTPUT_CAP,'Evaluation output cap exceeded; outputs preserved')
    guard(frozen)
    write_json(output/'summary.json',dict(status='complete_fixed_transport_evaluation',biological_n=17,regions=96460,
        inputs_unchanged_after_evaluation=True,protocol_sha256=sha(output/'protocol.json'),
        output_sha256={str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()},
        saved_output_bytes_before_summary=output_bytes,numpy=np.__version__,python=sys.version,
        source_person_disjointness_verified=False,independently_validated=False,model_fitted=False,
        primary_family='one registered paired mean-error contrast',
        interpretation='Qualified cross-assay transport; conditional donor bootstrap does not include refitting, measurement or assay-shift uncertainty'))
    print(json.dumps(dict(status='complete_fixed_transport_evaluation',biological_n=17,regions=96460,
                         mean_paired_error_reduction=result['mean_paired_error_reduction'],
                         percentile95=result['paired_error_reduction_percentile95'])),flush=True)


def main():
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Receiving evaluation only under SLURM')
    guard(CHECKS);guard({HERE/(n+'.py'):h for n,h in SOURCES.items()})
    # Successful synthetic receipt21998526 covers the exact current full-axis helper;
    # 21998524 records the earlier arithmetic checks, not current source lineage.
    require(read_json(REC/'codex_paired_profile_checks_21998526/summary.json')['script_sha256']==SOURCES['codex_paired_profile_evaluation'],
            'Current comparison helper lacks its successful full-axis receipt')
    output=REC/('codex_receiving_profile_evaluation_'+os.environ['SLURM_JOB_ID']);output.mkdir(exist_ok=False)
    try:
        modules=load_sources(output);evaluate_completed_transport(modules,output)
    except BaseException:
        write_json(output/'failure.json',dict(outputs_preserved=True,traceback=traceback.format_exc()))
        raise


if __name__=='__main__':
    main()
