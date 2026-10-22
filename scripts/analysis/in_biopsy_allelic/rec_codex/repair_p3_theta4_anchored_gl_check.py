#!/usr/bin/env python3
"""Three fixed synthetic theta4 points with isolated anchored BB and original GL."""
import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = Path(__file__).resolve().parent
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
CAMPAIGN = REC/'p3_parallel_20261003T130323Z'
MODE_RECEIPT = 'repair_p3_bb_mode_anchor_check_22042221'
DECIMAL_RECEIPTS = {1:'repair_p3_decimal_positive_gene1_22042086',
                    2:'repair_p3_decimal_positive_gene2_22042087',
                    3:'repair_p3_decimal_positive_gene3_22042088'}
ORIGINAL_RECEIPTS = {1:'repair_p3_bfix_vectors_22040217',
                     2:'repair_p3_bfix_vectors_22040218',
                     3:'repair_p3_bfix_vectors_22040219'}
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
GL_SHA = '45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c'
MODE_SHA = 'a31258aaac0fb39e91dda7aa634f0e810a6d0da144511b34d9c1b73a7028c266'
MODULES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py')
CAP = 32*1024**2


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def metadata(args):
    pins = {Path(__file__):args.expected_checker_sha,
        HERE/'repair_p3_integral.py':KERNEL_SHA,HERE/'repair_p3_composite_gl.py':GL_SHA,
        HERE/'repair_p3_bb_mode_anchor.py':MODE_SHA,
        FIXTURE/'SHA256SUMS':'0d3d06fc32a8e588d5e1c54c943a582e5054b6d2ab565290a4d57ced545d2925',
        CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
        args.receipt_pins:args.expected_receipt_pins_sha}
    for path,digest in pins.items():
        require(sha(path)==digest,'Reviewed source or metadata differs: '+path.name)
    campaign = json.loads((CAMPAIGN/'campaign.json').read_text())
    for name in MODULES:
        pins[HERE/name] = campaign['files'][str(CAMPAIGN/'executed_sources'/name)]
    published = {name.removeprefix('./'):digest for digest,name in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','genes.tsv','tags.tsv',
                 'params/theta4_perturbed.json'):
        pins[FIXTURE/name] = published[name]
    receipts = json.loads(args.receipt_pins.read_text())['receipt_sha256']
    expected = {directory+'/'+file for directory in (MODE_RECEIPT,*DECIMAL_RECEIPTS.values())
                for file in ('protocol.json','summary.json')}
    expected.update(directory+'/'+file for gene,directory in ORIGINAL_RECEIPTS.items()
                    for file in ('protocol.json','theta4_perturbed_gene%d.json' % gene))
    require(set(receipts)==expected,'Exact mode, all three Decimal and original refusal receipt pins required')
    for relative,digest in receipts.items():
        pins[REC/relative] = digest
    for path,digest in pins.items():
        require(sha(path)==digest,'Required completed receipt/source/input missing or differs: '+path.name)
    mode = json.loads((REC/MODE_RECEIPT/'summary.json').read_text())
    require(mode['required_independent_points']==18 and mode['completed_independent_points']==18 and
            mode['selected_decimal_family_criterion_met'] is True and mode['no_adoption'] is True,
            'Mode BB selected 18/18 Decimal family has not passed; broader accuracy not implied')
    mode_protocol = json.loads((REC/MODE_RECEIPT/'protocol.json').read_text())
    protocols,summaries,original_protocols,original_records = {},{},{},{}
    for gene,directory in DECIMAL_RECEIPTS.items():
        protocols[gene] = json.loads((REC/directory/'protocol.json').read_text())
        summaries[gene] = json.loads((REC/directory/'summary.json').read_text())
        summary,protocol = summaries[gene],protocols[gene]
        require(summary['completed'] is True and summary['gene']==gene and
                summary['case']=='theta4_perturbed' and len(summary['records'])==6 and
                {(r['decimal_digits'],r['order']) for r in summary['records']}==
                {(d,n) for d in (80,110) for n in (20,40,80)} and len(summary['refinements'])==7,
                'Completed six-record Decimal family required for every target gene')
        require(protocol['theta']=='theta4_perturbed' and protocol['gene']==gene and
                protocol['fixed_v']==.0004 and protocol['fixture_seed']==20260930,
                'Decimal fixed target differs')
        original_directory = REC/ORIGINAL_RECEIPTS[gene]
        original_protocols[gene] = json.loads((original_directory/'protocol.json').read_text())
        original_records[gene] = json.loads((original_directory/('theta4_perturbed_gene%d.json' % gene)).read_text())
        old = original_records[gene]
        require('candidate_refusal' in old and 'candidate' not in old and old['theta']=='theta4_perturbed' and
                old['gene']==gene and old['p']==protocol['fixed_p'] and old['v']==protocol['fixed_v'] and
                old['rows']==protocol['fixed_gene_rows'] and
                original_protocols[gene]['requested_points']==[['theta4_perturbed',gene]],
                'Original refusal or fixed-target identity differs')
    for protocol,key in [(mode_protocol,'source_input_pins')]+[
            (p,'source_and_input_sha256') for p in (*protocols.values(),*original_protocols.values())]:
        for absolute,digest in protocol[key].items():
            path = Path(absolute)
            require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve()),
                    'Receipt source/input path escaped project')
            require(path not in pins or pins[path]==digest,'Conflicting pinned source/input')
            pins[path] = digest
    for path,digest in pins.items():
        require(sha(path)==digest,'Completed receipt source/input changed: '+path.name)
        if path.suffix=='.py':
            ast.parse(path.read_text())
    return pins,mode_protocol,protocols,summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene',type=int,choices=(1,2,3),required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--receipt-pins',type=Path,required=True)
    parser.add_argument('--expected-receipt-pins-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    pins,mode_protocol,decimal_protocols,decimal_summaries = metadata(args)
    protocol = dict(theta='theta4_perturbed',gene=args.gene,fixture_seed=20260930,
        selected_BB_family_pass_is_not_broader_adoption=True,
        source_and_input_sha256={str(p):d for p,d in pins.items()},
        original_GL_orders=[16,32,64],original_GL_grid_spacings=[.25,.125],
        maximum_evaluations=100000,output_cap_bytes=CAP,
        no_model_fits=True,no_calibration=True,full_Hessian_accuracy_verified=False)
    if args.inspect:
        print(json.dumps(protocol,sort_keys=True))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    states = {}
    for job in (22042221,22042086,22042087,22042088):
        command = ['scontrol','show','job','-o',str(job)]
        terminal = subprocess.run(command,check=True,capture_output=True,text=True,timeout=15)
        fields = dict(token.split('=',1) for token in terminal.stdout.split() if '=' in token)
        require(fields.get('JobId')==str(job) and fields.get('JobState')=='COMPLETED' and
                fields.get('ExitCode')=='0:0',
                'Admission requires available terminal COMPLETED exit-zero mode and Decimal jobs; no implicit wait')
        states[str(job)] = dict(command=command,state=fields['JobState'],exit_code=fields['ExitCode'])
    protocol['admission_job_states'] = states
    output = REC/('repair_p3_theta4_anchored_gl_gene%d_%s' % (args.gene,os.environ['SLURM_JOB_ID']))
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    preserved = output/'original_refusal_receipts'
    preserved.mkdir()
    old_directory = REC/ORIGINAL_RECEIPTS[args.gene]
    for name in ('protocol.json','theta4_perturbed_gene%d.json' % args.gene):
        source = old_directory/name
        target = preserved/name
        target.write_bytes(source.read_bytes())
        require(sha(target)==pins[source],'Preserved original refusal receipt differs')
    archived = {}
    for source,digest in pins.items():
        if source.suffix=='.py':
            target = archive/source.name
            require(target.name not in archived or archived[target.name]==digest,
                    'Different source bytes share an archive basename')
            if not target.exists():
                target.write_bytes(source.read_bytes())
            require(sha(target)==digest,'Archive source differs')
            archived[target.name] = digest
    began = time.monotonic()

    def serial(value):
        if isinstance(value,np.ndarray):
            return value.tolist()
        if isinstance(value,np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    def diagnostic(value):
        if isinstance(value,np.ndarray):
            return diagnostic(value.tolist())
        if isinstance(value,np.generic):
            return diagnostic(value.item())
        if isinstance(value,dict):
            return {k:diagnostic(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)):
            return [diagnostic(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):
            return {'nonfinite_float':'nan' if np.isnan(value) else 'positive_infinity' if value>0 else 'negative_infinity'}
        return value

    def save(name,value):
        body = json.dumps(value,default=serial,indent=2,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())<=CAP-65536,
                'Diagnostic output cap exhausted; failure reserve retained')
        partial,final = output/(name+'.json.partial'),output/(name+'.json')
        with partial.open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
        require(not final.exists(),'Refusing result overwrite')
        partial.rename(final)
        print(json.dumps(dict(saved=name,elapsed_seconds=time.monotonic()-began)),flush=True)

    try:
        sys.path.insert(0,str(archive))
        require(not any(n[:-3] in sys.modules for n in MODULES),'Original numerical module already imported')
        import numpy as np
        import scipy
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES:
            module = sys.modules[name[:-3]]
            require(Path(module.__file__).resolve().parent==archive.resolve() and
                    sha(module.__file__)==pins[HERE/name],'Original import escaped pinned archive')
        private = load('anchored_gl_private_original',archive/'likelihood_v1.py')
        anchor = load('anchored_gl_private_BB',archive/'repair_p3_bb_mode_anchor.py')
        gl = load('anchored_gl_private_integrator',archive/'repair_p3_composite_gl.py')
        private.bb_terms = anchor.stable_bb_terms
        versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(versions==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),'Reviewed environment differs')
        design,settings = loader.load()
        model = design.model()
        require(len(design.participants)==60 and len(model.rows)==173 and model.G==8 and
                model.Z==7 and model.size==28 and settings['seed']==20260930,'Shared synthetic design differs')
        p,v = loader.parameters(FIXTURE/'params/theta4_perturbed.json')
        rows = tuple(r for r in model.rows if r.gene==args.gene)
        count_sha = hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest()
        row_body = json.dumps([vars(r) for r in rows],sort_keys=True,allow_nan=False)
        row_sha = hashlib.sha256(row_body.encode()).hexdigest()
        parent = decimal_protocols[args.gene]
        require(p.shape==(28,) and np.isfinite(p).all() and v==.0004 and p.tolist()==parent['fixed_p'] and
                count_sha==parent['full_fixture_count_vector_sha256'] and row_sha==parent['selected_row_metadata_sha256'] and
                len(rows)==parent['fixed_gene_rows'] and p.tolist()==mode_protocol['p'] and
                v==mode_protocol['v'] and count_sha==mode_protocol['counts_sha'],'Original target p/v/rows/counts differ')
        protocol.update(environment=versions,fixed_p=p.tolist(),fixed_v=v,fixed_gene_rows=len(rows),
            full_fixture_count_vector_sha256=count_sha,selected_row_metadata_sha256=row_sha)
        save('protocol',protocol)
        candidate = private.Model(rows,gene_count=8)
        refused = False
        try:
            ll,g,h,report = gl.integrate_gene(candidate,p,v,rows,private.conditional_terms,with_hessian=True,
                expected_kernel_sha=KERNEL_SHA,kernel_path=archive/'repair_p3_integral.py')
            require(g.shape==(29,) and h.shape==(29,29) and np.isfinite(ll) and np.isfinite(g).all() and
                    np.isfinite(h).all(),'Anchored GL result nonfinite or wrong shape')
            save('candidate',dict(logmass=ll,gradient=g,hessian=h,integration=diagnostic(report)))
            comparisons = []
            with localcontext() as context:
                context.prec = 110
                for record in decimal_summaries[args.gene]['records']:
                    comparisons.append(dict(order=record['order'],decimal_digits=record['decimal_digits'],
                        logmass_difference=str(Decimal.from_float(float(ll))-Decimal(record['logmass'])),
                        variance_score_difference=str(Decimal.from_float(float(g[-1]))-Decimal(record['variance_score'])),
                        Hvv_difference=str(Decimal.from_float(float(h[-1,-1]))-Decimal(record['Hvv'])),
                        comparison_is_diagnostic_not_full_Hessian_accuracy=True))
            save('Decimal_comparisons',comparisons)
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):
            refused = True
            failure = sys.exception() if hasattr(sys,'exception') else sys.exc_info()[1]
            save('candidate_refusal',dict(traceback=traceback.format_exc(),
                partial_reports=diagnostic(getattr(failure,'integration_reports',None)),
                evaluations=getattr(failure,'likelihood_evaluations',None),original_refusals_preserved=True))
        for path,digest in pins.items():
            require(sha(path)==digest,'Pinned source/input/receipt changed during diagnostic')
        save('summary',dict(completed=True,candidate_refused=refused,theta='theta4_perturbed',gene=args.gene,
            full_Hessian_accuracy_verified=False,no_model_fit=True,no_calibration=True,
            elapsed_seconds=time.monotonic()-began,all_original_accuracy_requirements_remain_required=True))
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())<=CAP,'Final output cap exceeded')
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(completed=False,traceback=traceback.format_exc(),partials_preserved=True,
                no_model_fit=True,no_calibration=True),handle,indent=2)
            handle.write('\n')
        raise


if __name__=='__main__':
    main()
