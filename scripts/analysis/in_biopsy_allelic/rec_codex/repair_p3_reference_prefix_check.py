#!/usr/bin/env python3
"""Focused synthetic literal-versus-prefix probability checks, compute only.

No fits, quadrature, calibration, biological inputs, or candidate evaluator use.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
MODULE = Path(__file__).resolve().with_name('repair_p3_reference_prefix.py')
LITERAL_SOURCE = REC / 'repair_p3_launch_20261005T121740Z/repair_p3_fixed_point_check.py'
LITERAL_SHA = '89434fd7e32b06cb93d9d31f18b646eaa646736ec329ab48cc00fa4409561f72'
CAP = 16 * 1024**2
STEPS = (1e-3, 3e-4, 1e-4)
VALUE_TOL = 1e-8


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-helper-sha')
    parser.add_argument('--expected-module-sha')
    args = parser.parse_args()
    ast.parse(MODULE.read_text())
    require(sha(LITERAL_SOURCE) == LITERAL_SHA, 'Frozen literal-reference source differs')
    if args.inspect:
        print(json.dumps(dict(numerical_execution=False, helper_sha256=sha(Path(__file__)),
              module_sha256=sha(MODULE), literal_source_sha256=LITERAL_SHA,
              output_cap_bytes=CAP, tests='Full support, normalized selected mixtures, fixed-step FD noise; no fits.')))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    own = {Path(__file__).resolve(): args.expected_helper_sha, MODULE: args.expected_module_sha}
    require(all(d and sha(p) == d for p, d in own.items()), 'Reviewed source drift')
    output = REC / ('repair_p3_reference_prefix_check_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    archive = output / 'executed_sources'
    archive.mkdir()
    began = time.monotonic()
    records, failures = [], []
    phase = 'archive_before_import'

    def save(name, value):
        body = json.dumps(value, indent=2, allow_nan=False) + '\n'
        size = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
        require(size + len(body.encode()) <= CAP - 65536, 'Output cap reached')
        partial = output / (name + '.json.partial')
        final = output / (name + '.json')
        with partial.open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
        require(not final.exists(), 'Refusing result overwrite')
        partial.rename(final)

    def case(kind, identity, work):
        try:
            result = work()
            row = dict(kind=kind, identity=identity, completed=True, **result)
            if not row.get('value_checks_passed', True):
                failures.append(dict(kind=kind, identity=identity, reason='value_tolerance'))
        except Exception as exc:
            row = dict(kind=kind, identity=identity, completed=False,
                       exception_type=type(exc).__name__, traceback=traceback.format_exc())
            failures.append(dict(kind=kind, identity=identity, reason='exception'))
        records.append(row)

    try:
        for source, expected in own.items():
            target = archive / source.name
            with target.open('xb') as handle:
                handle.write(source.read_bytes())
            require(sha(target) == expected, 'Archive source differs')
        with (archive / LITERAL_SOURCE.name).open('xb') as handle:
            handle.write(LITERAL_SOURCE.read_bytes())
        require(sha(archive / LITERAL_SOURCE.name) == LITERAL_SHA, 'Literal archive differs')
        phase = 'numerical_import'
        import numpy as np
        import scipy
        from scipy.special import expit, gammaln, logsumexp
        from scipy.stats import binom
        spec = importlib.util.spec_from_file_location('isolated_reference_prefix', archive / MODULE.name)
        prefix = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = prefix
        spec.loader.exec_module(prefix)
        require(Path(prefix.__file__).resolve().parent == archive.resolve(), 'Reference import escaped archive')
        environment = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        require(environment == dict(python='3.10.18', numpy='2.2.6', scipy='1.15.2'), 'Expected numerical environment differs')

        # Literal implementation transcribed from the frozen checker; no candidate code.
        def rise(logshape, count):
            return 0. if count == 0 else logshape + np.logaddexp(logshape, np.log(np.arange(1,count,dtype=float))).sum()

        def literal(n, eta, rr):
            rho = 1e-6 + (1-2e-6)*expit(rr)
            logt = np.log((1-rho)/rho)
            la, lb = logt-np.logaddexp(0.,-eta), logt-np.logaddexp(0.,eta)
            return np.array([gammaln(n+1)-gammaln(a+1)-gammaln(n-a+1)+rise(la,int(a))+rise(lb,int(n-a))-rise(logt,n)
                             for a in range(n+1)])

        def fixed(n, mean, phi):
            if phi == 0:
                return binom.logpmf(np.arange(n+1), n, mean)
            t = (1-phi)/phi
            return np.array([gammaln(n+1)-gammaln(a+1)-gammaln(n-a+1)+rise(np.log(mean*t),int(a))+
                             rise(np.log((1-mean)*t),int(n-a))-rise(np.log(t),n) for a in range(n+1)])

        timing = dict(literal_seconds=0., prefix_unprepared_seconds=0., prefix_prepared_seconds=0., evaluations=0)
        phase = 'full_support'
        for n in (0, 5, 20, 200, 500):
            for eta in (-1000., -40., -6., 0., 6., 40., 1000.):
                for rr in (-15., -3., 5.):
                    def full(n=n, eta=eta, rr=rr):
                        t = time.perf_counter(); old = literal(n,eta,rr); timing['literal_seconds'] += time.perf_counter()-t
                        t = time.perf_counter(); new = prefix.logpmf_full(n,eta,rr); timing['prefix_unprepared_seconds'] += time.perf_counter()-t
                        prepared = prefix.PreparedBB(n,rr)
                        t = time.perf_counter(); reused = prepared.logpmf(eta); timing['prefix_prepared_seconds'] += time.perf_counter()-t
                        timing['evaluations'] += 1
                        reflected = prefix.logpmf_full(n,-eta,rr)[::-1]
                        error = float(np.max(np.abs(old-new)))
                        norm_old, norm_new = float(logsumexp(old)), float(logsumexp(new))
                        symmetry = float(np.max(np.abs(new-reflected)))
                        reuse_error = float(np.max(np.abs(new-reused)))
                        return dict(max_logpmf_difference=error, old_log_normalizer=norm_old, new_log_normalizer=norm_new,
                            reflection_max_error=symmetry, preparation_reuse_max_error=reuse_error,
                            value_checks_passed=bool(np.all(np.isfinite(old)) and np.all(np.isfinite(new)) and
                                max(error,abs(norm_old),abs(norm_new),symmetry,reuse_error) <= VALUE_TOL))
                    case(phase,dict(n=n,eta=eta,linear_r=rr),full)
        save('full_support', records)
        records = []
        phase = 'selected_mixtures_and_fd'
        for n in (20, 200, 500):
            lo = max(1,(10*n+99)//100)
            support = np.arange(lo,n-lo+1)
            for phi in (0., .02):
                homref, homalt = fixed(n,.01,phi)[support], fixed(n,.99,phi)[support]
                for pref,palt,px in ((.8,.15,.05),(.5,.5,0.),(1.,0.,0.)):
                    for orientation in (-1,1):
                        prior = np.array([2*pref*palt,pref**2+2*pref*px,palt**2+2*palt*px])
                        if orientation == -1:
                            prior = prior[[0,2,1]]
                        # Raw three-class weights exclude x/x; the selected-support
                        # normalizer handles their common scale exactly as the original.
                        require(abs(prior.sum()-(1-px**2)) <= 1e-14, 'Raw three-class prior differs')
                        with np.errstate(divide='ignore'):
                            lp = np.log(prior)
                        for eta in (-6.,0.,6.):
                            for rr in (-15.,-3.,5.):
                                def mixture(eta,rr,fast):
                                    effective = eta-orientation*.3  # fixed row omega; no second orientation shift
                                    het = prefix.logpmf_full(n,effective,rr) if fast else literal(n,effective,rr)
                                    logmass = logsumexp(np.stack((het[support],homref,homalt))+lp[:,None],axis=0)
                                    return logmass-logsumexp(logmass)
                                def selected(eta=eta,rr=rr):
                                    old,new = mixture(eta,rr,False),mixture(eta,rr,True)
                                    error = float(np.max(np.abs(old-new)))
                                    norm = float(abs(logsumexp(new)))
                                    fd = []
                                    # Fixed family: both coordinates, all three steps, selected endpoints/middle.
                                    # Compare noise only; no derivative pass threshold or extrapolated certification.
                                    for axis in ('eta','linear_r'):
                                        for h in STEPS:
                                            args_plus = (eta+h,rr) if axis=='eta' else (eta,rr+h)
                                            args_minus = (eta-h,rr) if axis=='eta' else (eta,rr-h)
                                            d_old = (mixture(*args_plus,False)-mixture(*args_minus,False))/(2*h)
                                            d_new = (mixture(*args_plus,True)-mixture(*args_minus,True))/(2*h)
                                            at = (0,len(support)//2,len(support)-1)
                                            fd.append(dict(axis=axis,step=h,max_fd_difference=float(np.max(np.abs(d_old-d_new))),
                                                observed_counts=[int(support[i]) for i in at],old_scores=[float(d_old[i]) for i in at],
                                                prefix_scores=[float(d_new[i]) for i in at],
                                                near_zero_score_binding_atol=1e-8,binding_rtol=1e-6,derivative_pass_not_asserted=True))
                                    return dict(selected_support_low=lo,selected_support_high=n-lo,
                                        max_selected_logpmf_difference=error,selected_log_normalizer_error=norm,
                                        finite_difference_records=fd,value_checks_passed=bool(np.all(np.isfinite(old)) and
                                            np.all(np.isfinite(new)) and max(error,norm)<=VALUE_TOL))
                                identity = dict(n=n,phi=phi,p_ref=pref,p_alt=palt,p_x=px,raw_prior_sum=float(prior.sum()),
                                                excluded_x_x_weight=px**2,orientation=orientation,eta=eta,linear_r=rr)
                                case(phase,identity,selected)
        save('selected_mixtures_and_fd',records)
        save('summary',dict(completed=True,value_checks_passed=not failures,failures=failures,
            full_support_cases=105,selected_mixture_cases=len(records),fd_records=sum(len(r.get('finite_difference_records',[])) for r in records),
            timing=timing,environment=environment,seed=20261005,random_draws_used=False,elapsed_seconds=time.monotonic()-began,
            source_sha256={str(p):d for p,d in own.items()},literal_source_sha256=LITERAL_SHA,
            output_cap_bytes=CAP,probability_tolerance=VALUE_TOL,fd_pass_gate_encoded=False,
            model_fitting=False,calibration=False,numerical_integration=False,
            limitations=['Prefix cumsum and literal sum use different floating reduction order.',
                'Probability algebra is shared rising-factorial identity; independent integration and mixture assembly remain separate.',
                'FD comparisons diagnose arithmetic noise, not analytic derivative correctness or quadrature accuracy.',
                'Timing is separate evidence and does not establish numerical accuracy.']))
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())<=CAP,'Final inclusive output cap exceeded')
        require(not failures,'Recorded probability check failures; inspect all case receipts')
        print(json.dumps(dict(completed=True,output=str(output),value_checks_passed=True,timing=timing)),flush=True)
    except Exception:
        save('failure',dict(completed=False,phase=phase,traceback=traceback.format_exc(),recorded_failures=failures,
                            elapsed_seconds=time.monotonic()-began))
        raise


if __name__ == '__main__':
    main()
