#!/usr/bin/env python3
"""Independent seven-by-seven Hessian stencil reference; compute only, no fits."""
import argparse
import ast
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback
import warnings

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
CAMPAIGN = REC / 'p3_parallel_20261003T130323Z'
PROTOCOL = REC / 'repair_p3_case_protocol_20261005.json'
HERE = Path(__file__).resolve().parent
SNAPSHOT = REC / 'repair_p3_complex_launch_20261005T141156Z'
REFERENCE = SNAPSHOT / 'repair_p3_complex_reference.py'
CHECKER_SOURCE = SNAPSHOT / 'repair_p3_complex_check.py'
CHECKER_SHA = 'ca9630e4f939ed1684eb81a67ad337b56bab289679f66466adca1e366a87c1fd'
COMPOSITES = {'failed0gene2':REC/'repair_p3_fixed_point_failed0gene2_22039200',
              'positive18':REC/'repair_p3_fixed_point_positive18_22039201'}
REFERENCE_SHA = 'f722eb76a2fb288f4f4e697dcd149d714b89221ec791ba6a43db8cfebb2f7ec3'
SCALARS = {'positive18':REC/'repair_p3_scalar_prefix_positive18_22039119',
           'failed0gene2':REC/'repair_p3_scalar_prefix_failed0gene2_22039118'}
BASE_KERNEL = HERE / 'repair_p3_integral.py'
BASE_KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
ANCESTOR_SHA = '89434fd7e32b06cb93d9d31f18b646eaa646736ec329ab48cc00fa4409561f72'
ANCESTOR = REC / 'repair_p3_launch_20261005T121740Z/repair_p3_fixed_point_check.py'
CAP = 32 * 1024**2
CASES = {'failed0gene2':0,'positive18':18}
G = 2
STEPS = (1e-3, 3e-4, 1e-4)
HESSIAN_STEPS = (1e-3,3e-4,1e-4,3e-5,1e-5)
COMPLEX_STEPS = (1e-10,1e-15,1e-20)
FIXED_AXES = (2,10,16,17,24,25,28)
GUARDS = {
    PROTOCOL: '8e7ddf8b0a7bb18afd9692b53c3560f4258fdc373c4b822aee36aa94395eabe9',
    CAMPAIGN / 'campaign.json': 'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
    REFERENCE: REFERENCE_SHA,
    CHECKER_SOURCE: CHECKER_SHA,
    COMPOSITES['failed0gene2']/'candidate_baseline.json':'edec063451383f9a1386619186c2089d3ef8e6ee34d854e80a2894fdd28e0f9a',
    COMPOSITES['failed0gene2']/'protocol.json':'7b5d793be5485a25e5fc7f5b3d97c4ddf30a1a497d697728cd64468ed61f8d40',
    COMPOSITES['failed0gene2']/'summary.json':'5ad1ea571c71525e6ae1eb18bf625baf1e6779934d64a149e0d16ca5f3137a0b',
    COMPOSITES['positive18']/'candidate_baseline.json':'4a576370be2f33ac28bf23052fe283b4a5f614bc30e983bd6a84224da857c722',
    COMPOSITES['positive18']/'protocol.json':'ed6c190424fa29501e1377b5b47b12762917e53f46338c3c3f1a408c3f9e4774',
    COMPOSITES['positive18']/'summary.json':'1daa3275c8d524857c90369b7686fc2234229df765e41f139a4b232c9e323613',
    SCALARS['positive18'] / 'independent_baseline.json': '0184daa0afae8e392272868530d31ad064089c3de06483d6cc8b87b76b109257',
    SCALARS['positive18'] / 'summary.json': 'fba2f38821c59e85ea4973e5ebe8cfe2509bc35d9910fcb7d9638078d6adb43d',
    SCALARS['positive18'] / 'protocol.json': '26f061f8e481a1b9d1e11f2a5dfc0984c82010d3b1dd9c8f270eb1418d701001',
    SCALARS['positive18'] / 'finite_differences.json': 'ebce73a49bbca829ff182e4cb0fc710c18fea794edbe8f8e476a7d145c91ac62',
    SCALARS['failed0gene2'] / 'independent_baseline.json': 'b6f8422966c6db5b10ccbe498a30e1e899f92ec16c66736f71a6c08a7e685d1a',
    SCALARS['failed0gene2'] / 'summary.json': '74f375cf748940e4bf91137e55ab6e2e8d564b7fedb48032c19d60c19e688d53',
    SCALARS['failed0gene2'] / 'protocol.json': 'b27a85f6c11342406e07e99aacafa70839d09fcd870e0449dac40e59aab97a32',
    SCALARS['failed0gene2'] / 'finite_differences.json': '2ae8d822a8e456ad6e247ffaf91bcc21f5a90e06b63c66c8d2846a0cda96451a',
    ANCESTOR: ANCESTOR_SHA,
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def metadata():
    for path, expected in GUARDS.items():
        require(sha(path) == expected, 'Pinned small protocol/source differs')
    protocol = json.loads(PROTOCOL.read_text())
    guards = dict(GUARDS)
    guards.update({ROOT / p: digest for p, digest in protocol['case_input_sha256'].items()})
    for path, expected in guards.items():
        require(sha(path) == expected, 'Pinned case metadata differs')
    campaign = json.loads((CAMPAIGN / 'campaign.json').read_text())
    # Inspect hashes only for source code here; numerical arrays are compute-only.
    for path, expected in campaign['files'].items():
        p = Path(path)
        if p.suffix == '.py':
            require(p.parent == CAMPAIGN / 'executed_sources' and sha(p) == expected, 'Immutable synthetic source differs')
            guards[p] = expected
    ast.parse(REFERENCE.read_text())
    return protocol, campaign, guards


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=tuple(CASES), required=True)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-helper-sha')
    parser.add_argument('--expected-reference-sha')
    parser.add_argument('--axes',nargs='+',type=int,choices=FIXED_AXES)
    args = parser.parse_args()
    require(args.axes is None or len(set(args.axes))==len(args.axes),'Axes must be a unique subset')
    columns=tuple(args.axes) if args.axes is not None else FIXED_AXES
    score_axes=FIXED_AXES
    SCALAR = SCALARS[args.case]
    protocol, campaign, guards = metadata()
    if args.inspect:
        print(json.dumps(dict(case=args.case,draw=CASES[args.case],fixed_gene=G,
            numerical_execution=False,no_model_fits=True,no_calibration=True,
            independent_complex_gradient_stencils=True,requested_columns=columns,score_axes=score_axes,
            all_required_columns=FIXED_AXES,real_steps=HESSIAN_STEPS,original_scalar_FD_steps=STEPS,
            complex_steps=COMPLEX_STEPS,bounds=(24.,26.),
            per_realpoint_scoreaxis_evaluation_cap=100000,
            aggregate_cap_formula='distinct_realpoints * 7 * 100000',
            maximum_central_realpoints_per_column=20,maximum_one_sided_realpoints_per_column=21,
            expected_hessian_records_per_column=len(score_axes)*len(HESSIAN_STEPS)*len(COMPLEX_STEPS)*2,
            estimate_cpu_hours_per_column=[4,18],estimate_peak_memory_gib=.5,
            reference_sha256=sha(REFERENCE),checker_ancestry_sha256=CHECKER_SHA,
            unverified_scope='All29moment rows/columns outside selected7x7, full fits and calibration.',
            metadata_and_source_checks=True,hessian_accuracy_acceptance=False)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute job required')
    require(args.expected_helper_sha==sha(Path(__file__)) and args.expected_reference_sha==sha(REFERENCE),'Reviewed complex source drift')
    output = REC / ('repair_p3_independent_hessian_' + args.case + '_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    archive = output / 'executed_sources'
    archive.mkdir()
    began, phase = time.monotonic(), 'archive_before_import'

    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    def save(name, value):
        body = json.dumps(value, indent=2, default=serial, allow_nan=False) + '\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file()) + len(body.encode()) <= CAP-65536,
                'Check output cap reached; failure reserve retained')
        final = output / (name + '.json')
        partial = output / (name + '.json.partial')
        with partial.open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
        require(not final.exists(), 'Refusing result overwrite')
        partial.rename(final)

    def event(name, **fields):
        record = dict(phase=phase, event=name, elapsed_seconds=time.monotonic()-began, **fields)
        with (output / 'progress.jsonl').open('a') as handle:
            handle.write(json.dumps(record, default=serial, allow_nan=False)+'\n')
        print(json.dumps(record, default=serial), flush=True)

    try:
        own = {Path(__file__):args.expected_helper_sha, REFERENCE:args.expected_reference_sha}
        for source in (*own, *(p for p in guards if p.suffix == '.py')):
            (archive / source.name).write_bytes(source.read_bytes())
            require(sha(archive / source.name) == (own.get(source) or guards[source]), 'Archive copy differs')
        fixture = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
        fixture_names = ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','draws/p3_uniforms.npy')
        for name in fixture_names:
            path = fixture / name
            require(path.resolve().is_relative_to(fixture.resolve()), 'Fixture input escaped synthetic root')
            expected = campaign['files'][str(path)]
            require(sha(path) == expected, 'Synthetic fixture bytes differ')
            guards[path] = expected
        phase = 'archived_imports_and_synthetic_draw'
        sys.path.insert(0, str(archive))
        require(not any(n in sys.modules for n in ('likelihood_v1','bfix_v1','design_v1','tests_v1')), 'Numerical original already imported')
        global np
        import numpy as np
        import scipy
        loader = importlib.import_module('bfix_v1')
        original = importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')  # explicit before the archived-module assertion
        for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve().parent == archive.resolve() and
                    sha(Path(sys.modules[name].__file__)) == campaign['files'][str(CAMPAIGN/'executed_sources'/ (name+'.py'))],
                    'Original import escaped guarded archive')
        reference = load_module('repair_p3_archived_complex_reference',archive/REFERENCE.name)
        require(Path(reference.__file__).resolve().parent==archive.resolve(),'Complex reference import escaped archive')
        require(reference.STEPS==COMPLEX_STEPS,'Complex step family differs')
        versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        require(versions == protocol['environment'], 'Original numerical environment versions differ')
        design, settings = loader.load()
        model = design.model()
        require(settings['B'] == 50 and len(model.rows) == 173, 'Synthetic design shape differs')
        observed = json.loads((CAMPAIGN/'baseline/observed_null.json').read_text())
        require(observed['converged'], 'Saved generating null not converged')
        null_p = np.asarray(observed['parameters'])
        uniforms = np.load(fixture/'draws/p3_uniforms.npy', mmap_mode='r', allow_pickle=False)
        require(uniforms.shape == (50,173,2), 'Original supplied uniforms differ')
        require([(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')] ==
                [(r.individual,settings['gene_order'][r.gene]) for r in model.rows], 'Original uniform-to-row order differs')
        draw = CASES[args.case]
        rows = []
        for row, (uc, uk) in zip(model.rows, uniforms[draw]):
            eta, r, _, _ = model.predictors(row, null_p)
            rows.append(original.conditional_draw(row, eta, r, uc, uk))
        drawn = original.Model(rows, gene_count=model.G)
        count_sha = hashlib.sha256(np.asarray([r.a for r in rows], dtype='<i8').tobytes()).hexdigest()
        if draw == 0:
            require(count_sha == '5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691', 'Known draw0 reconstruction differs')
            failed = json.loads((REC/'p3_first_draw_diagnostic_21998696/summary.json').read_text())
            p, v = np.asarray(failed['latest_marginal']['parameters']), float(failed['latest_marginal']['v'])
            require(v == 0.9234818604998617, 'Wrong original failed point')
        else:
            record = json.loads((CAMPAIGN/f'baseline/p3_{draw:03d}.json').read_text())
            require(record['b'] == draw and record['fit']['converged'], 'Wrong original successful point')
            p, v = np.asarray(record['fit']['parameters']), float(record['fit']['v'])
            require(v == (0.00711885901713953 if draw == 18 else 0.), 'Wrong original successful variance')
        gene_rows = [r for r in drawn.rows if r.gene == G]
        require(len(gene_rows) == 24 and len(p) == drawn.size, 'Fixed gene/coordinate shape differs')
        require(v>0,'Fixed complex reference requires saved positive variance')
        scalar_protocol=json.loads((SCALAR/'protocol.json').read_text())
        scalar_summary=json.loads((SCALAR/'summary.json').read_text())
        scalar_baseline=json.loads((SCALAR/'independent_baseline.json').read_text())
        scalar_fd=json.loads((SCALAR/'finite_differences.json').read_text())
        require(len(scalar_fd)==21 and {(r['axis'],r['step']) for r in scalar_fd}==
                {(a,h) for a in FIXED_AXES for h in STEPS}, 'Full fixed scalar FD family differs')
        require(scalar_summary['scalar_fixed_point_check_completed'] and scalar_summary['scalar_fd_record_count']==21 and
                scalar_summary['case']==args.case and scalar_summary['draw']==draw, 'Complete scalar reference case differs')
        require(scalar_protocol['fixed_p']==p.tolist() and scalar_protocol['fixed_v']==v and
                scalar_protocol['count_vector_sha256']==count_sha, 'Scalar reference point or row counts differ')
        vector=np.r_[p,v]
        bounds=drawn.bounds()+[(0.,None)]
        require(len(bounds)==len(vector),'Bound axis differs')
        def valid_vector(q):
            require(np.isfinite(q).all() and q[-1]>0,'Invalid positive real-offset point')
            for value,(lo,hi) in zip(q,bounds):
                require((lo is None or value>=lo) and (hi is None or value<=hi),'Offset violates original parameter bounds')
        valid_vector(vector)
        points=[]
        point_index={}
        def register(q):
            valid_vector(q)
            body=np.asarray(q,dtype='<f8').tobytes()
            digest=hashlib.sha256(body).hexdigest()
            if digest not in point_index:
                point_index[digest]=len(points)
                points.append(dict(point_id=len(points),vector=q.copy(),vector_sha256=digest))
            return point_index[digest]
        stencils=[]
        stencil_refusals=[]
        for column in columns:
            for real_step in HESSIAN_STEPS:
                x=vector[column]
                lo,hi=bounds[column]
                central=(lo is None or x-2*real_step>=lo) and (hi is None or x+2*real_step<=hi)
                if central:
                    method='fourth_order_central'
                    offsets=(-2,-1,1,2)
                    coefficients=(1,-8,8,-1)
                elif hi is None or x+4*real_step<=hi:
                    method='fourth_order_forward_one_sided'
                    offsets=(0,1,2,3,4)
                    coefficients=(-25,48,-36,16,-3)
                elif lo is None or x-4*real_step>=lo:
                    method='fourth_order_backward_one_sided'
                    offsets=(0,-1,-2,-3,-4)
                    coefficients=(25,-48,36,-16,3)
                else:
                    stencil_refusals.append(dict(column=column,real_step=real_step,
                        reason='No declared fourth-order stencil fits unchanged bounds'))
                    continue
                # Integer polynomial moments independently verify the stencil:
                # d/dx x=1 and derivatives at0 of constants/powers2..4 are0.
                for degree in range(5):
                    require(sum(c*m**degree for c,m in zip(coefficients,offsets))==(12 if degree==1 else 0),
                            'Declared stencil fails exact polynomial identity')
                ids=[]
                for offset in offsets:
                    q=vector.copy()
                    q[column]+=offset*real_step
                    ids.append(register(q))
                stencils.append(dict(column=column,real_step=real_step,method=method,
                    offsets=offsets,coefficients=coefficients,denominator=12*real_step,point_ids=ids))
        aggregate_cap=len(points)*len(score_axes)*100000
        save('protocol',dict(case=args.case,draw=draw,fixed_gene=G,fixed_gene_rows=24,fixed_p=p,fixed_v=v,
            count_vector_sha256=count_sha,environment=versions,
            source_guards={str(k):d for k,d in guards.items()},own_source_hashes={str(k):d for k,d in own.items()},
            requested_columns=columns,all_required_columns=FIXED_AXES,score_axes=score_axes,
            real_steps=HESSIAN_STEPS,original_scalar_FD_steps=STEPS,
            complex_steps=COMPLEX_STEPS,bounds=(24.,26.),grid_spacings=(.5,.25),
            original_parameter_bounds=bounds,stencils=stencils,stencil_refusals=stencil_refusals,
            point_roster=points,distinct_realpoints=len(points),
            per_realpoint_scoreaxis_evaluation_cap=100000,aggregate_successful_evaluation_cap=aggregate_cap,
            shared_geometry_mass_physically_attributed_once=True,
            each_scoreaxis_budget_counter_includes_shared_geometry_mass=True,
            no_original_conditional_evaluator_or_scores_used_for_reference=True,
            hessian_reference_scope='Selected seven score axes by requested selected columns only.',
            all29moment_coordinates_outside7x7_independently_unverified=True,
            no_model_fit=True,no_calibration=True,output_cap_bytes=CAP,
            estimate_cpu_hours_per_column=[4,18],estimate_peak_memory_gib=.5,
            maximum_planned_realpoints_per_column=21,
            supplemental_real_steps_prespecified_before_independent_hessian_results=True,
            supplemental_step_rationale='Small positive variance may need finer fourth-order stencils; all five steps retained.',
            estimates_not_measured_for_realoffset_family=True,hessian_accuracy_acceptance=False))
        point_results={}
        total_successful_evaluations=0
        refusals=list(stencil_refusals)
        phase='realoffset_independent_reference'
        for point in points:
            point_id=point['point_id']
            q=point['vector']
            shared=reference.Reference(gene_rows,model.G,model.Z,maximum_evaluations=100000)
            geometries,masses={},{}
            point_refusals=[]
            for bound,spacing in ((24.,.5),(26.,.25)):
                try:
                    require(shared.evaluations<100000,'Shared realpoint evaluation cap exhausted')
                    geometry=shared.geometry(q[:-1],q[-1],bound,spacing)
                    mass=shared.mass(q[:-1],q[-1],geometry)
                    geometries[bound],masses[bound]=geometry,mass
                except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
                    record=dict(point_id=point_id,bound=bound,completed=False,phase='real_geometry_mass',
                        exception_type=type(error).__name__,traceback=traceback.format_exc(),
                        reference_record=getattr(error,'reference_record',None),
                        budget_counter_evaluations=shared.evaluations)
                    point_refusals.append(record)
                    refusals.append(record)
            shared_count=min(shared.evaluations,100000)
            total_successful_evaluations+=shared_count
            save('point_%03d_metadata' % point_id,dict(point=point,masses=masses,
                shared_geometry_mass_successful_probability_evaluations=shared_count,
                shared_counter_including_refused_attempt=shared.evaluations,refusals=point_refusals))
            values={}
            axis_receipts=[]
            for score_axis in score_axes:
                engine=reference.Reference(gene_rows,model.G,model.Z,maximum_evaluations=100000)
                # Match executed single-axis check scope: account shared baseline
                # in each100k axis budget, but add its physical calls only once.
                engine.evaluations=shared_count
                records=[]
                for complex_step in COMPLEX_STEPS:
                    for bound in (24.,26.):
                        identity=dict(point_id=point_id,score_axis=score_axis,complex_step=complex_step,bound=bound)
                        try:
                            require(bound in masses,'Derivative blocked by realpoint mass refusal')
                            require(engine.evaluations<100000,'Realpoint/scoreaxis100k budget exhausted')
                            result=engine.derivative(q[:-1],q[-1],score_axis,complex_step,geometries[bound],masses[bound])
                            require(np.isfinite(result['gradient']) and all(np.isfinite(result[key]) and result[key]>=0
                                for key in ('estimated_combined_gradient_uncertainty','estimated_numerator_denominator_gradient_error',
                                            'relative_derivative_tail_allowance','finite_domain_variance_endpoint_allowance')),
                                'Independent gradient or uncertainty is nonfinite/invalid')
                            values[(score_axis,complex_step,bound)]=result
                            record=dict(**identity,completed=True,result=result,
                                axis_budget_counter=engine.evaluations,shared_baseline_budget_count=shared_count)
                        except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
                            record=dict(**identity,completed=False,exception_type=type(error).__name__,
                                traceback=traceback.format_exc(),reference_record=getattr(error,'reference_record',None),
                                axis_budget_counter=engine.evaluations,shared_baseline_budget_count=shared_count)
                            refusals.append(record)
                        records.append(record)
                successful_axis_count=max(0,min(engine.evaluations,100000)-shared_count)
                total_successful_evaluations+=successful_axis_count
                axis_receipt=dict(point_id=point_id,score_axis=score_axis,records=records,
                    maximum_axis_budget_evaluations=100000,
                    shared_baseline_budget_evaluations=shared_count,
                    derivative_successful_probability_evaluations=successful_axis_count,
                    axis_budget_counter_including_refused_attempt=engine.evaluations,
                    physical_shared_evaluations_not_added_again=True)
                axis_receipts.append(axis_receipt)
                save('point_%03d_score_%02d' % (point_id,score_axis),axis_receipt)
                event('realpoint_scoreaxis_saved',point_id=point_id,score_axis=score_axis,
                    successful_records=sum(r['completed'] for r in records),
                    total_successful_probability_evaluations=total_successful_evaluations)
                require(total_successful_evaluations<=aggregate_cap,'Explicit aggregate evaluation cap exceeded')
            point_results[point_id]=dict(values=values,axis_receipts=axis_receipts)
        phase='independent_stencils_before_candidate_comparison'
        hessian_records=[]
        for stencil in stencils:
            for score_axis in score_axes:
                for complex_step in COMPLEX_STEPS:
                    for bound in (24.,26.):
                        identity=dict(column=stencil['column'],score_axis=score_axis,real_step=stencil['real_step'],
                                      complex_step=complex_step,bound=bound)
                        try:
                            selected=[point_results[pid]['values'][(score_axis,complex_step,bound)] for pid in stencil['point_ids']]
                            coefficients=stencil['coefficients']
                            denominator=stencil['denominator']
                            estimate=sum(c*r['gradient'] for c,r in zip(coefficients,selected))/denominator
                            error=sum(abs(c)*r['estimated_combined_gradient_uncertainty'] for c,r in zip(coefficients,selected))/denominator
                            quad_error=sum(abs(c)*r['estimated_numerator_denominator_gradient_error'] for c,r in zip(coefficients,selected))/denominator
                            tail_error=sum(abs(c)*r['relative_derivative_tail_allowance'] for c,r in zip(coefficients,selected))/denominator
                            endpoint_error=sum(abs(c)*r['finite_domain_variance_endpoint_allowance'] for c,r in zip(coefficients,selected))/denominator
                            require(all(np.isfinite(value) for value in (estimate,error,quad_error,tail_error,endpoint_error)),
                                    'Independent Hessian stencil is nonfinite')
                            row=dict(**identity,completed=True,method=stencil['method'],stencil=stencil,
                                independent_hessian=float(estimate),abs_stencil_weighted_uncertified_error=float(error),
                                propagated_quad_estimate=float(quad_error),propagated_derivative_tail_allowance=float(tail_error),
                                propagated_variance_endpoint_allowance=float(endpoint_error),
                                input_gradients=[r['gradient'] for r in selected],
                                real_stencil_truncation_complex_truncation_arithmetic_unbounded=True)
                        except (KeyError,ValueError,FloatingPointError) as error:
                            row=dict(**identity,completed=False,method=stencil['method'],stencil=stencil,
                                     reason='At least one required independent-gradient record unavailable',
                                     traceback=traceback.format_exc())
                            refusals.append(row)
                        hessian_records.append(row)
        for refusal in stencil_refusals:
            for score_axis in score_axes:
                for complex_step in COMPLEX_STEPS:
                    for bound in (24.,26.):
                        row=dict(column=refusal['column'],real_step=refusal['real_step'],score_axis=score_axis,
                                 complex_step=complex_step,bound=bound,completed=False,reason=refusal['reason'])
                        hessian_records.append(row)
                        refusals.append(row)
        expected_record_count=len(columns)*len(score_axes)*len(HESSIAN_STEPS)*len(COMPLEX_STEPS)*2
        save('independent_hessian_before_candidate_comparison',hessian_records)
        # Candidate values are read only after all independent reference arithmetic.
        phase='compare_saved_candidate_only'
        composite_protocol=json.loads((COMPOSITES[args.case]/'protocol.json').read_text())
        composite_summary=json.loads((COMPOSITES[args.case]/'summary.json').read_text())
        for key in ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256','environment'):
            require(composite_protocol[key]==scalar_protocol[key],'Saved candidate point differs: '+key)
        require(composite_summary['candidate_integral_returned'],'Saved composite baseline not complete')
        candidate_record=json.loads((COMPOSITES[args.case]/'candidate_baseline.json').read_text())
        candidate_h=candidate_record['hessian']
        candidate_error=candidate_record['integration']['empirical_refinement_plus_tail_allowances'][2]
        require(len(candidate_h)==len(vector) and all(len(row)==len(vector) for row in candidate_h),'Candidate Hessian shape differs')
        comparisons=[]
        for row in hessian_records:
            if not row['completed']:
                continue
            i,j=row['score_axis'],row['column']
            difference=row['independent_hessian']-candidate_h[i][j]
            combined=row['abs_stencil_weighted_uncertified_error']+candidate_error[i][j]
            comparisons.append(dict(**row,candidate_hessian=candidate_h[i][j],
                independent_minus_candidate=difference,
                unchanged_numerical_hessian_budget=1e-5,
                nominal_budget_within_descriptive_only=abs(difference)<=1e-5,
                candidate_plus_reference_uncertified_error=combined,
                uncertainty_overlap_descriptive_only=abs(difference)<=combined,
                hessian_accuracy_not_accepted=True))
        save('candidate_comparisons',comparisons)
        refinements=[]
        for column in columns:
            for score_axis in score_axes:
                rows=[r for r in hessian_records if r['completed'] and r['column']==column and r['score_axis']==score_axis]
                pairs=[]
                for i,left in enumerate(rows):
                    for right in rows[i+1:]:
                        differing=sum(left[k]!=right[k] for k in ('real_step','complex_step','bound'))
                        if differing!=1:
                            continue
                        difference=right['independent_hessian']-left['independent_hessian']
                        error=right['abs_stencil_weighted_uncertified_error']+left['abs_stencil_weighted_uncertified_error']
                        pairs.append(dict(first={k:left[k] for k in ('real_step','complex_step','bound')},
                            second={k:right[k] for k in ('real_step','complex_step','bound')},
                            hessian_difference=difference,summed_uncertified_error=error,
                            disagreement_outside_estimated_error=abs(difference)>error,
                            agreement_not_an_accuracy_pass=True))
                refinements.append(dict(column=column,score_axis=score_axis,pairs=pairs))
        save('real_complex_step_domain_refinement',refinements)
        save('original_scalar_FD_records',scalar_fd)
        phase='final_source_hashes'
        for path,expected in {**guards,**own}.items():
            require(sha(path)==expected,'Source/fixture changed during instrumented check')
        save('summary',dict(independent_hessian_records_completed=True,case=args.case,draw=draw,fixed_gene=G,
            requested_columns=columns,score_axes=score_axes,required_combined_seven_coordinate_family=FIXED_AXES,
            real_steps=HESSIAN_STEPS,original_scalar_FD_steps=STEPS,
            complex_steps=COMPLEX_STEPS,bounds=(24.,26.),
            expected_hessian_record_count=expected_record_count,saved_hessian_record_count=len(hessian_records),
            successful_hessian_record_count=sum(r['completed'] for r in hessian_records),refusal_count=len(refusals),
            stencil_refusals=stencil_refusals,
            all_requested_records_present=len(hessian_records)==expected_record_count,
            requested_columns_numerically_complete=not refusals and len(hessian_records)==expected_record_count,
            full7x7_numerically_complete=not refusals and len(hessian_records)==expected_record_count and set(columns)==set(FIXED_AXES),
            all29moment_entries_outside7x7_independently_unverified=True,
            hessian_accuracy_accepted=False,no_original_conditional_scores_used_for_reference=True,
            candidate_hessian_used_only_after_reference_finished=True,
            distinct_realpoints=len(points),per_realpoint_scoreaxis_evaluation_cap=100000,
            aggregate_successful_probability_evaluation_cap=aggregate_cap,
            actual_successful_probability_evaluations=total_successful_evaluations,
            shared_geometry_mass_physically_counted_once=True,
            numerical_hessian_budget_unchanged=1e-5,original_scalar_FD_records_retained=len(scalar_fd),
            no_model_fits=True,no_statistical_calibration=True,original_records_preserved=True,
            uncertainty_estimates_not_certified=True,real_complex_stencil_truncation_arithmetic_not_bounded=True,
            elapsed_seconds=time.monotonic()-began,output_cap_bytes=CAP,
            saved_bytes_before_summary=sum(q.stat().st_size for q in output.rglob('*') if q.is_file())))
        event('completed')
        require(sum(q.stat().st_size for q in output.rglob('*') if q.is_file())<=CAP,'Final inclusive output cap exceeded')
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(completed=False,phase=phase,traceback=traceback.format_exc(),partials_preserved=True,
                           no_fit_or_calibration=True),handle,indent=2)
            handle.write('\n')
        raise


if __name__=='__main__':
    main()


