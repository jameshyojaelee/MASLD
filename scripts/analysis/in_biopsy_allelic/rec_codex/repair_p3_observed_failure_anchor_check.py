#!/usr/bin/env python3
"""Isolated arithmetic diagnostic at saved synthetic observed-fit call 12."""
import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import math
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
FAILED = REC/'repair_p3_parallel_observed_fit_22042114'
MODE = REC/'repair_p3_bb_mode_anchor_check_22042221'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
COUNT_SHA = '7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
MODULES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py')
CAP = 32*1024**2
HELPERS = {
    'repair_p3_integral.py':KERNEL_SHA,
    'repair_p3_composite_gl.py':'45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c',
    'repair_p3_bb_mode_anchor.py':'a31258aaac0fb39e91dda7aa634f0e810a6d0da144511b34d9c1b73a7028c266',
    'repair_p3_decimal_positive.py':'6ff2857ce42e408478f9c5dc7cc13ed6c70bcfdd77bd8f2169fff8de9656ac40',
    'repair_p3_boundary_decimal_fast.py':'d192f1aa311b54d3fe94af3744cf14a8c02856a95ffc5e3e590c7efa1432c017',
    'repair_p3_bb_mode_anchor_check.py':'559eb29393f1165dc77608b4d54845ca3bad624ae4697f48dedcd75c2f12f684'}


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


def metadata(expected_checker_sha):
    pins = {Path(__file__):expected_checker_sha,
        FAILED/'failure.json':'d3e7946dc977906052a34a36a12fc9898ba0d82fd77497a52df41f83f8a7f305',
        FAILED/'protocol.json':'05c69cdb3aa293a7b1150bfc16c0a28487d76d731307236341e4e2c0de5b0604',
        MODE/'summary.json':'070d7fd02720ce736bd65aaaf7e25fbd163da28bfb26a655fe517fbda5283846',
        MODE/'protocol.json':'74a3e719068757a45dd776bf08f6decf0f450d0c025bfcffc34064f6b8014e88',
        FIXTURE/'SHA256SUMS':'0d3d06fc32a8e588d5e1c54c943a582e5054b6d2ab565290a4d57ced545d2925',
        CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94'}
    pins.update({HERE/name:digest for name,digest in HELPERS.items()})
    for path,digest in pins.items():
        require(sha(path)==digest,'Reviewed source/input/receipt differs: '+path.name)
    campaign = json.loads((CAMPAIGN/'campaign.json').read_text())
    for name in MODULES:
        pins[HERE/name] = campaign['files'][str(CAMPAIGN/'executed_sources'/name)]
    published = {name.removeprefix('./'):digest for digest,name in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','genes.tsv','tags.tsv'):
        pins[FIXTURE/name] = published[name]
    failure = json.loads((FAILED/'failure.json').read_text())
    parent = json.loads((FAILED/'protocol.json').read_text())
    mode = json.loads((MODE/'summary.json').read_text())
    mode_protocol = json.loads((MODE/'protocol.json').read_text())
    require(mode['required_independent_points']==18 and mode['completed_independent_points']==18 and
            mode['selected_decimal_family_criterion_met'] is True and mode['no_adoption'] is True,
            'Previously selected BB family has not passed; no broader acceptance implied')
    point = failure['last_integration_failure']
    require(failure['phase']=='observed_marginal_fit' and failure['events_saved']==12 and
            failure['numeric_failure_details_omitted_due_to_output_cap'] is False and
            point['call']==12 and point['gene']==2 and point['variance']==.0009636462009520182 and
            len(point['parameters'])==28 and len(point['refused_gene_ladders'])==6 and
            point['no_full_sum_formed'] is True,
            'Exact saved observed failure identity or retained ladder differs')
    require(parent['observed_count_vector_sha256']==COUNT_SHA and parent['fixture_seed']==20260930,
            'Original observed synthetic fixture identity differs')
    for protocol,key in ((parent,'source_guards'),(mode_protocol,'source_input_pins')):
        for absolute,digest in protocol[key].items():
            path = Path(absolute)
            require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve()),
                    'Receipt source/input escaped project')
            require(path not in pins or pins[path]==digest,'Conflicting source/input pins')
            pins[path] = digest
    for path,digest in pins.items():
        require(sha(path)==digest,'Pinned source/input changed: '+path.name)
        if path.suffix=='.py':
            ast.parse(path.read_text())
    return pins,parent,point


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    pins,parent,point = metadata(args.expected_checker_sha)
    protocol = dict(original_job=22042114,call=12,gene=2,fixed_p=point['parameters'],fixed_v=point['variance'],
        fixture_seed=20260930,count_vector_sha256=COUNT_SHA,participants=60,rows=173,target_rows=24,
        source_and_input_sha256={str(p):d for p,d in pins.items()},
        representative_points=['u=0','u=+sqrt(v)','u=-sqrt(v)'],maximum_decimal_n=5000,
        maximum_decimal_prefix_log_evaluations=20000000,
        selected_BB_criteria=dict(logpmf='1e-9',score='1e-8',Hessian='1e-5',Decimal_stability='1e-25'),
        full_GL_with_hessian=False,original_GL_orders=[16,32,64],maximum_GL_evaluations=100000,
        Decimal_orders=[20,40,80],Decimal_digits=[80,110],output_cap_bytes=CAP,
        all_comparisons_diagnostic=True,no_model_fits=True,no_calibration=True,
        full_Hessian_accuracy_verified=False,no_mathematical_repair_claim=True)
    if args.inspect:
        print(json.dumps(protocol,sort_keys=True))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    command = ['scontrol','show','job','-o','22042221']
    terminal = subprocess.run(command,check=True,capture_output=True,text=True,timeout=15)
    fields = dict(token.split('=',1) for token in terminal.stdout.split() if '=' in token)
    require(fields.get('JobId')=='22042221' and fields.get('JobState')=='COMPLETED' and
            fields.get('ExitCode')=='0:0','Selected BB admission job must be available and completed exit zero')
    protocol['selected_BB_admission'] = dict(command=command,state=fields['JobState'],exit_code=fields['ExitCode'])
    output = REC/('repair_p3_observed_failure_anchor_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    # Only the named original mathematical modules and six reviewed helpers are
    # executable dependencies. Parent launchers remain pinned metadata, avoiding
    # basename collisions between different historical launcher revisions.
    for source in (Path(__file__),*(HERE/name for name in (*MODULES,*HELPERS))):
        target = archive/source.name
        target.write_bytes(source.read_bytes())
        require(sha(target)==pins[source],'Executable archive differs')
    preserved = output/'original_failure_receipts'
    preserved.mkdir()
    for name in ('failure.json','protocol.json'):
        target = preserved/name
        target.write_bytes((FAILED/name).read_bytes())
        require(sha(target)==pins[FAILED/name],'Original full failure packet differs')
    began = time.monotonic()

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
        body = json.dumps(diagnostic(value),indent=2,allow_nan=False)+'\n'
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
        require(not any(n[:-3] in sys.modules for n in MODULES),'Original numerical modules already imported')
        import numpy as np
        import scipy
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES:
            module = sys.modules[name[:-3]]
            require(Path(module.__file__).resolve().parent==archive.resolve() and
                    sha(module.__file__)==pins[HERE/name],'Original import escaped reviewed archive')
        private = load('observed_anchor_private_original',archive/'likelihood_v1.py')
        anchor = load('observed_anchor_private_BB',archive/'repair_p3_bb_mode_anchor.py')
        gl = load('observed_anchor_private_GL',archive/'repair_p3_composite_gl.py')
        decimal_positive = load('observed_anchor_Decimal_positive',archive/'repair_p3_decimal_positive.py')
        decimal_bb = load('observed_anchor_Decimal_BB',archive/'repair_p3_boundary_decimal_fast.py')
        bb_checker = load('observed_anchor_BB_regression',archive/'repair_p3_bb_mode_anchor_check.py')
        private.bb_terms = anchor.stable_bb_terms
        versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(versions==parent['environment']==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),
                'Original observed numerical environment differs')
        design,settings = loader.load()
        model = design.model()
        rows = tuple(r for r in model.rows if r.gene==2)
        count_sha = hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest()
        require(len(design.participants)==60 and len(model.rows)==173 and model.G==8 and model.Z==7 and
                model.size==28 and len(rows)==24 and count_sha==COUNT_SHA and settings['seed']==20260930,
                'Exact observed synthetic rows/design differ')
        p,v = np.asarray(point['parameters']),float(point['variance'])
        require(p.shape==(28,) and np.isfinite(p).all() and v>0,'Saved parameter vector invalid')
        protocol.update(environment=versions,selected_row_metadata_sha256=hashlib.sha256(
            json.dumps([vars(r) for r in rows],sort_keys=True,allow_nan=False).encode()).hexdigest())
        representative = max(rows,key=lambda row:row.n)
        require(representative.n<=5000 and 3*representative.n*3*36<=20000000,
                'Representative Decimal work cap exceeded')
        protocol['representative_depth'] = representative.n
        protocol['decimal_prefix_log_evaluation_upper_bound'] = 3*representative.n*3*36
        save('protocol',protocol)
        independent = []
        for index,u in enumerate((0.,math.sqrt(v),-math.sqrt(v))):
            eta,r = model.predictors(representative,p,u)[:2]
            # Native predictors omit cohort omega; BB conditional applies this
            # orientation correction once. This pure BB comparison does likewise.
            eta -= representative.orientation*representative.omega
            record = dict(n=representative.n,k=representative.a,u=u,eta=eta,r=r,
                individual=representative.individual,selected_family_criterion_met=False)
            try:
                lp,g,H = anchor.stable_bb_terms([representative.a],representative.n,eta,r)
                record['candidate'] = dict(logpmf=lp,score=g,Hessian=H)
                ladder = bb_checker.decimal_fd(representative.n,representative.a,eta,r,decimal_bb.bb_vector)
                record['Decimal_reference'] = ladder
                final = ladder[-1]
                with localcontext() as context:
                    context.prec = 110
                    def flatten(a):
                        return [Decimal(a['logpmf']),*map(Decimal,a['score']),
                            *(Decimal(x) for row in a['hessian'] for x in row)]
                    differences = [[str(x-y) for x,y in zip(flatten(a),flatten(final))] for a in ladder]
                    stability = max(abs(Decimal(x)) for row in differences for x in row)
                    raw = [Decimal.from_float(float(lp[0]))-Decimal(final['logpmf']),
                        *(Decimal.from_float(float(x))-Decimal(y) for x,y in zip(g[0],final['score'])),
                        *(Decimal.from_float(float(x))-Decimal(y) for x,y in zip(H[0].flat,
                            [x for row in final['hessian'] for x in row]))]
                    record.update(raw_candidate_minus_Decimal=list(map(str,raw)),
                        Decimal_ladder_minus_final=differences,max_Decimal_stability_difference=str(stability))
                    record['selected_family_criterion_met'] = bool(np.isfinite(lp).all() and np.isfinite(g).all() and
                        np.isfinite(H).all() and abs(raw[0])<=Decimal('1e-9') and
                        max(map(abs,raw[1:3]))<=Decimal('1e-8') and max(map(abs,raw[3:]))<=Decimal('1e-5') and
                        stability<=Decimal('1e-25'))
            except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                record['numerical_refusal'] = dict(traceback=traceback.format_exc())
            independent.append(record)
            save('representative_BB_'+str(index),record)
        selected_pass = all(r['selected_family_criterion_met'] for r in independent)
        refused = None
        comparisons = []
        if selected_pass:
            candidate = private.Model(rows,gene_count=8)
            gl_result = None
            try:
                ll,g,h,report = gl.integrate_gene(candidate,p,v,rows,private.conditional_terms,with_hessian=False,
                    expected_kernel_sha=KERNEL_SHA,kernel_path=archive/'repair_p3_integral.py')
                require(g.shape==(29,) and h is None and np.isfinite(ll) and np.isfinite(g).all(),
                        'GL score-only result invalid')
                gl_result = (ll,g)
                refused = False
                save('GL_candidate',dict(logmass=ll,gradient=g,Hessian=None,integration=report))
            except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                refused = True
                error = sys.exc_info()[1]
                save('GL_candidate_refusal',dict(traceback=traceback.format_exc(),
                    partial_reports=getattr(error,'integration_reports',None),
                    likelihood_evaluations=getattr(error,'likelihood_evaluations',None)))
            def seeds(order):
                return np.sqrt(2)*np.polynomial.hermite.hermgauss(order)[0]
            reference = decimal_positive.check(rows,p,v,8,7,decimal_bb.bb_vector,seeds,save)
            save('Decimal_reference_summary',reference)
            if gl_result is not None:
                with localcontext() as context:
                    context.prec = 110
                    for record in reference['records']:
                        comparisons.append(dict(order=record['order'],decimal_digits=record['decimal_digits'],
                            logmass_difference=str(Decimal.from_float(float(gl_result[0]))-Decimal(record['logmass'])),
                            variance_score_difference=str(Decimal.from_float(float(gl_result[1][-1]))-Decimal(record['variance_score'])),
                            diagnostic_only=True))
                save('GL_minus_Decimal',comparisons)
        for path,digest in pins.items():
            require(sha(path)==digest,'Pinned source/input/receipt changed during diagnostic')
        save('summary',dict(completed=True,selected_three_BB_points_pass=selected_pass,
            GL_and_Decimal_skipped_unless_selected_BB_points_pass=not selected_pass,
            GL_candidate_refused=refused,original_failed_fit_preserved=True,
            no_full_Hessian_accuracy=True,no_model_fit=True,no_calibration=True,no_mathematical_repair_claim=True,
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
