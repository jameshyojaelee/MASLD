#!/usr/bin/env python3
"""Extended-precision native conditional arithmetic at two saved points.

Variance scores are compared with the separately checked heat identities.
Full Hessians receive domain checks only, not independent validation.
No optimizer, resampling, calibration or human-data access.
"""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
POINT = REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420'
HEAT = {2: REC/'heat_variance_score_gene2_22084571', 5: REC/'heat_variance_score_gene5_22084572'}
GL = {2: REC/'GL_initial64_refusal_gene2_22084480/summary.json',
      5: REC/'GL_initial64_refusal_gene5_22084481/summary.json'}
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
           'repair_p3_bb_gamma_difference_v3.py', 'repair_p3_local_moments.py',
           'repair_p3_pathwise_moments_v1.py', 'repair_p3_longdouble_conditional_v1.py',
           'repair_p3_decimal_mixed.py', 'repair_p3_boundary_decimal_fast.py')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gene', type=int, choices=(2, 5), required=True)
    ap.add_argument('--pins', type=Path, required=True)
    ap.add_argument('--expected-pins-sha', required=True)
    ap.add_argument('--expected-runner-sha', required=True)
    ap.add_argument('--inspect', action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha, 'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    for path, digest in pins.items():
        p = Path(path)
        require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest, 'Input differs')
        if p.suffix == '.py':
            ast.parse(p.read_text())
    required = [*(HERE/n for n in SOURCES), POINT/'draw_counts.json', POINT/'recipe_summary.json',
                *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')),
                *(HEAT[g]/n for g in (2, 5) for n in ('summary.json', 'protocol.json')), *GL.values()]
    require(all(str(p) in pins for p in required), 'Missing pin')
    read = lambda p: json.loads(p.read_text())
    heat, hp = read(HEAT[a.gene]/'summary.json'), read(HEAT[a.gene]/'protocol.json')
    gl = read(GL[a.gene])
    failure = read(POINT/'recipe_summary.json')['fitting_integration_failure']
    require(heat['reference_check_passed'] and heat['gene'] == a.gene == hp['gene']
            and hp['parameters'] == failure['parameters'] == gl['parameters']
            and hp['variance'] == failure['variance'] == gl['variance'], 'Saved point differs')
    for name in SOURCES[:5]:
        require(hp['source_and_input_sha256'][str(HERE/name)] == pins[str(HERE/name)], 'Heat numerical source differs')
    require(hp['source_and_input_sha256'][str(POINT/'draw_counts.json')] == pins[str(POINT/'draw_counts.json')], 'Heat counts differ')
    if a.inspect:
        print(json.dumps(dict(gene=a.gene, source_inputs_checked=True, no_fit=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('pathwise_longdouble_component_budget_gene'+str(a.gene)+'_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    archive = out/'executed_sources'
    archive.mkdir()
    for p in (Path(__file__), *(HERE/n for n in SOURCES)):
        q = archive/p.name
        q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]), 'Archive differs')
    import numpy as np
    import scipy
    from scipy.integrate import quad_vec
    sys.path.insert(0, str(archive))
    original = importlib.import_module('likelihood_v1')
    bb = original.bb_terms
    loader = importlib.import_module('bfix_v1')
    gamma = importlib.import_module('repair_p3_bb_gamma_difference_v3')
    lm = importlib.import_module('repair_p3_local_moments')
    pm = importlib.import_module('repair_p3_pathwise_moments_v1')
    ld_module = importlib.import_module('repair_p3_longdouble_conditional_v1')
    mixed = importlib.import_module('repair_p3_decimal_mixed')
    probability = importlib.import_module('repair_p3_boundary_decimal_fast')
    require(loader.FIXTURE.resolve() == FIXTURE.resolve() and original.bb_terms is bb, 'Loader/arithmetic differs')
    for name in ('likelihood_v1', 'design_v1', 'bfix_v1', 'tests_v1', 'repair_p3_bb_gamma_difference_v3', 'repair_p3_local_moments', 'repair_p3_pathwise_moments_v1', 'repair_p3_longdouble_conditional_v1',
                 'repair_p3_decimal_mixed', 'repair_p3_boundary_decimal_fast'):
        require(Path(sys.modules[name].__file__).resolve().parent == archive.resolve(), 'Import escaped archive')
    design, settings = loader.load()
    template = design.model()
    counts = read(POINT/'draw_counts.json')
    require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows]
            and len(counts['counts']) == len(template.rows) == 173, 'Count order differs')
    model = original.Model([replace(r, a=c) for r, c in zip(template.rows, counts['counts'])], gene_count=8)
    rows = tuple(r for r in model.rows if r.gene == a.gene)
    p = np.asarray(failure['parameters'])
    p.setflags(write=False)
    require(model.G == 8 and model.Z == 7 and model.size == 28 and 0 < len(rows) <= 32
            and all(r.n <= 5000 for r in rows) and len(rows) == hp['rows'], 'Synthetic scope differs')
    local = lm.LocalGeneMoments(model, p, rows, None)
    moments = pm.PathwiseGeneMoments(local)
    context = ld_module.LongDoubleConditional(model, p, rows)
    v = failure['variance']
    before, row_before = p.tobytes(), repr(rows)
    def safe(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, dict): return {k: safe(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)): return [safe(v) for v in x]
        return x
    def save(name, value):
        with (out/(name+'.json')).open('x') as f:
            json.dump(safe(value), f, indent=2, allow_nan=False)
            f.write('\n'); f.flush(); os.fsync(f.fileno())
    save('protocol', dict(gene=a.gene, parameters=p, variance=v, rows=len(rows), source_and_input_sha256=pins,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__), node=os.environ.get('SLURMD_NODENAME'),
        derivative_units='d/dv; u=sqrt(v)*z; fixed standard normal measure', epsabs=1e-12, epsrel=1e-12,
        native_conditional_arithmetic='longdouble integer prefixes, direct selected support',
        longdouble_mantissa_bits=np.finfo(np.longdouble).nmant,
        separate_variance_second_moment_quadrature=True,
        scalar_second_moment_epsabs=1e-5, scalar_second_moment_epsrel=1e-13,
        no_fit=True, no_calibration=True, no_human_data=True, empirical_not_certified=True))
    start = time.perf_counter()
    records = []
    try:
        original.bb_terms = gamma.stable_bb_terms
        # Compare every saved row at three one-SD nodes before quadrature.
        maxima = np.zeros(3)
        for row in rows:
            for z in (-1., 0., 1.):
                eta, rr = model.predictors(row, p, math.sqrt(v)*z, need_hessian=False)[:2]
                ref = original.conditional_terms(row, eta, rr)
                actual = context.terms(row, eta, rr)
                maxima = np.maximum(maxima, [np.max(np.abs(np.asarray(a)-b)) for a, b in zip(actual, ref)])
        save('conditional_Gamma_probe', dict(rows=len(rows), nodes_per_row=3, maximum_differences=maxima))
        require(np.all(maxima <= [1e-9, 1e-8, 1e-7]), 'Extended conditional and Gamma disagree')
        # Independent probability-only Decimal finite differences on two rows.
        from decimal import Decimal, localcontext
        targets = [max(rows, key=lambda r: r.n)]
        opposite = [r for r in rows if r.orientation != targets[0].orientation]
        if opposite: targets.append(max(opposite, key=lambda r: r.n))
        fd_records = []
        with localcontext() as decimal_context:
            decimal_context.prec = 80
            prepared = mixed.prepare(targets, p, model.G, model.Z, probability.bb_vector)
            step = Decimal('1e-10')
            for d in prepared:
                row = d['row']
                ee, rr = d['eta0'], d['r']
                values = {(i, j): mixed.row_logmass(d, ee+i*step, rr+j*step, probability.bb_vector)
                          for i, j in ((0,0), (1,0), (-1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1))}
                f0 = values[0,0]
                gf = np.asarray([float((values[1,0]-values[-1,0])/(2*step)),
                                 float((values[0,1]-values[0,-1])/(2*step))])
                cross = float((values[1,1]-values[1,-1]-values[-1,1]+values[-1,-1])/(4*step*step))
                hf = np.asarray([[float((values[1,0]+values[-1,0]-2*f0)/(step*step)), cross],
                                 [cross, float((values[0,1]+values[0,-1]-2*f0)/(step*step))]])
                eta, rr_float = model.predictors(row, p, 0., need_hessian=False)[:2]
                actual = context.terms(row, eta, rr_float)
                diffs = [abs(actual[0]-float(f0)), float(np.max(np.abs(actual[1]-gf))), float(np.max(np.abs(actual[2]-hf)))]
                fd_records.append(dict(individual=row.individual, n=row.n, orientation=row.orientation,
                    decimal_digits=80, FD_step=str(step), differences=diffs))
        save('independent_Decimal_probe', dict(records=fd_records, probability_only_FD=True))
        require(all(all(a <= b for a, b in zip(r['differences'], (1e-9, 1e-8, 1e-7))) for r in fd_records),
                'Extended conditional and independent Decimal disagree')
        ell0 = math.fsum(context.terms(row, *model.predictors(row, p, 0., need_hessian=False)[:2],
                                      need_hessian=False)[0] for row in rows)
        _, g0, _, s0 = moments.point(p, 0., v, context.terms)
        width = len(g0)
        scales = np.maximum(1., np.abs(np.r_[1., g0, s0.ravel()]))
        # At z=0 every pathwise variance derivative is exactly zero. Include
        # both signs of a one-SD node when scaling heterogeneous units.
        for z in (-1., 1.):
            _, gs, _, ss = moments.point(p, z, v, context.terms)
            scales = np.maximum(scales, np.abs(np.r_[1., gs, ss.ravel()]))
        def integrand(z):
            ll, g, _, s = moments.point(p, float(z), v, context.terms)
            return math.exp(ll-ell0)*np.r_[1., g, s.ravel()]/scales
        for bound in (24., 26.):
            values, error, info = quad_vec(integrand, -bound, bound, epsabs=1e-12, epsrel=1e-12,
                norm='max', cache_size=4*1024**2, limit=1024, workers=1, full_output=True)
            save('quadrature'+str(int(bound)), dict(success=bool(info.success), status=int(info.status),
                message=info.message, evaluations=info.neval, error_estimate=float(error),
                scaled_result=values, scales=scales))
            require(info.success and np.isfinite(values).all() and math.isfinite(error) and error >= 0, 'Quadrature refused')
            values = values*scales
            # Integrate the variance second numerator in its own units. The
            # vector quadrature's common max-norm estimate is conservative
            # after multiplication by this component's large scale.
            component_errors = scales*error
            scalar, scalar_error, scalar_info = quad_vec(
                lambda z: integrand(z)[-1]*scales[-1], -bound, bound,
                epsabs=1e-5, epsrel=1e-13, norm='max', cache_size=4*1024**2,
                limit=1024, workers=1, full_output=True)
            save('variance_second_quadrature'+str(int(bound)), dict(
                success=bool(scalar_info.success), status=int(scalar_info.status),
                message=scalar_info.message, evaluations=scalar_info.neval,
                result=float(scalar), error_estimate=float(scalar_error),
                original_vector_result=float(values[-1]),
                original_vector_error_estimate=float(component_errors[-1])))
            require(scalar_info.success and math.isfinite(scalar) and math.isfinite(scalar_error)
                    and scalar_error >= 0, 'Variance second-moment quadrature refused')
            values[-1] = float(scalar)
            component_errors[-1] = float(scalar_error)
            mass = float(values[0])
            g = values[1:1+width]/mass
            second = values[1+width:].reshape(width, width)/mass
            h = second-np.outer(g, g)
            mt, gt, st = moments.shifted_tail_bounds(bound, v, ell0)
            lower = mass-error-mt
            require(lower > 0, 'Estimated mass allowance unresolved')
            ga = (component_errors[1:1+width]+gt+np.abs(g)*(error+mt))/lower
            sa = (component_errors[1+width:].reshape(width, width)+st+np.abs(second)*(error+mt))/lower
            ha = sa+np.outer(np.abs(g), ga)+np.outer(ga, np.abs(g))+np.outer(ga, ga)
            record = dict(z_bound=bound, logmass=ell0+math.log(mass), gradient=local.embed_vector(g),
                hessian=local.embed_matrix(h), estimated_logmass_allowance=(error+mt)/lower,
                estimated_gradient_allowance=local.embed_vector(ga), estimated_hessian_allowance=local.embed_matrix(ha),
                shifted_mass=mass, shifted_moments=values, scales=scales, log_shift=ell0,
                absolute_quadrature_component_error_estimates=component_errors,
                absolute_shifted_tail_bounds=[mt, gt, st], quadrature_error_estimate=float(error), quad_evaluations=info.neval)
            records.append(record)
            save('bound'+str(int(bound)), record)
    finally:
        original.bb_terms = bb
    left, right = records
    dl = abs(left['logmass']-right['logmass'])
    dg = np.abs(left['gradient']-right['gradient'])
    dh = np.abs(left['hessian']-right['hessian'])
    stable = (dl+left['estimated_logmass_allowance']+right['estimated_logmass_allowance'] <= 1e-9
        and np.all(dg+left['estimated_gradient_allowance']+right['estimated_gradient_allowance'] <= 1e-8)
        and np.all(dh+left['estimated_hessian_allowance']+right['estimated_hessian_allowance'] <= 1e-5))
    heat_deltas = [right['gradient'][-1]-r['domain'][-1]['variance_score'] for r in heat['records']]
    heat_pass = all(abs(d) <= 1e-8 for d in heat_deltas)
    require(p.tobytes() == before and repr(rows) == row_before and original.bb_terms is bb, 'Inputs changed')
    for path, digest in pins.items(): require(sha(Path(path)) == digest, 'Input changed during check')
    save('summary', dict(completed=True, gene=a.gene, domain_LL_difference=dl, domain_gradient_max_difference=float(dg.max()),
        domain_hessian_max_difference=float(dh.max()), domain_allowance_check_passed=bool(stable),
        heat_variance_score_differences=heat_deltas, heat_variance_score_check_passed=heat_pass,
        GL_gradient_max_difference=float(np.max(np.abs(right['gradient']-gl['gradient']))),
        fixed_point_check_passed=bool(stable and heat_pass), seconds=time.perf_counter()-start,
        no_independent_full_Hessian_validation=True, empirical_not_certified=True, no_fit=True, no_calibration=True))
    print(json.dumps(dict(gene=a.gene, fixed_point_check_passed=bool(stable and heat_pass))))
    if not (stable and heat_pass): raise SystemExit(1)


if __name__ == '__main__':
    main()
