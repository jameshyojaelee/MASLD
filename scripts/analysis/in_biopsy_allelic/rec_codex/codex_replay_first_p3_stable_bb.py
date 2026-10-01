"""Compute-only replay of the first failed synthetic P3 draw, isolated.

Only BB numerical representation is replaced in an immutable private copy.
Original model, bounds, starts, fitting rules and failed-draw treatment remain.
No live module, running job, fixture or historical result is changed.
"""
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

import numpy as np
import scipy

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ORIGINAL = REC / 'bfix_p3_21998519'
SNAPSHOT = ORIGINAL / 'executed_sources'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
PRIOR = REC / 'bfix_fits_21984981/null.json'
HASHES = {
    'likelihood_v1.py': '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
    'bfix_v1.py': '324594ce783ba2464921eaee5924065f579060eded6043a48bdf4b7e392bfa56',
    'design_v1.py': 'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
    'tests_v1.py': '52cd14d43b649c20347330b64b467dbd81320f10252fe85ce3c538ceef2f092e',
    'check_bfix_tests_v13.py': 'a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d',
}
PROBE_SHA = '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488'
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv', 'draws/p3_uniforms.npy')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    job = os.environ['SLURM_JOB_ID']
    output = REC / ('codex_replay_first_p3_stable_bb_' + job)
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'
    private.mkdir()
    def save(name, value):
        (output / (name+'.json')).write_text(json.dumps(value, default=serial, indent=2)+'\n')
        print(name+' saved', flush=True)

    started = time.monotonic()
    probe = Path(__file__).with_name('codex_stable_bb_probe.py')
    require(sha(probe) == PROBE_SHA, 'Reviewed stable arithmetic source changed')
    for name, expected in HASHES.items():
        require(sha(SNAPSHOT / name) == expected, 'Immutable source hash mismatch: '+name)
        (private / name).write_bytes((SNAPSHOT / name).read_bytes())
    (private / probe.name).write_bytes(probe.read_bytes())
    (private / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    stable_spec = importlib.util.spec_from_file_location('private_stable_bb_probe', private / probe.name)
    stable_module = importlib.util.module_from_spec(stable_spec)
    sys.modules[stable_spec.name] = stable_module
    stable_spec.loader.exec_module(stable_module)
    stable_bb_terms = stable_module.stable_bb_terms
    original_receipt = json.loads((ORIGINAL / 'receipt.json').read_text())
    published = {}
    for line in (FIXTURE / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        published[name.removeprefix('./')] = digest
    inputs = {name: sha(FIXTURE / name) for name in INPUTS}
    require(all(digest == published[name] == original_receipt['fixture_sha256'][name]
                for name, digest in inputs.items()), 'Original supplied synthetic input hash mismatch')
    require(sha(PRIOR) == original_receipt['prior_null_sha256'], 'Original warm-start fit changed')
    reference = {name: json.loads((ORIGINAL / (name+'.json')).read_text())
                 for name in ('observed_null', 'observed_marginal', 'p3_000')}
    require(reference['p3_000']['fit'] is None and reference['p3_000']['exceed'] is True
            and np.isnan(reference['p3_000']['LR']), 'Historical failed draw differs')
    receipt = dict(draw_index=0, supplied_draws=50, replayed_draws=1,
                   fixture=str(FIXTURE), fixture_sha256=inputs,
                   source_sha256=HASHES, stable_probe_sha256=PROBE_SHA,
                   wrapper_sha256=sha(Path(__file__)), prior_null_sha256=sha(PRIOR),
                   original_results_sha256={name: sha(ORIGINAL / (name+'.json')) for name in reference},
                   original_failed_exceedance=reference['p3_000'],
                   modification='Only snapshot bb_terms numerical representation uses exact log-shape identity',
                   law='Original 2-start fit_marginal, shifts, bounds, section4_from_null, statistic and failure rule',
                   scientific_calibration=False, actual_genotypes_read=False,
                   receiving_outcomes_read=False, python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, job_id=job)
    save('receipt', receipt)

    # Import dependencies exclusively from copied immutable sources.
    sys.path.insert(0, str(private))
    likelihood = importlib.import_module('likelihood_v1')
    original_spec = importlib.util.spec_from_file_location('original_p3_likelihood', private / 'likelihood_v1.py')
    original_module = importlib.util.module_from_spec(original_spec)
    sys.modules[original_spec.name] = original_module
    original_spec.loader.exec_module(original_module)
    phase = 'load'
    trace = {'first_warning': None, 'first_exception': None, 'first_extreme_bb_call': None,
             'bb_calls_by_phase': {}, 'eta_min': None, 'eta_max': None}
    standard_warning = warnings.showwarning
    def record_warning(message, category, filename, lineno, file=None, line=None):
        if trace['first_warning'] is None:
            trace['first_warning'] = dict(phase=phase, message=str(message), category=category.__name__,
                                         filename=filename, lineno=lineno, stack=traceback.format_stack())
            save('numerical_trace', trace)
        standard_warning(message, category, filename, lineno, file=file, line=line)
    warnings.showwarning = record_warning
    def stable_wrapper(k, n, eta, linear_r, need_hessian=True):
        trace['bb_calls_by_phase'][phase] = trace['bb_calls_by_phase'].get(phase, 0)+1
        trace['eta_min'] = float(eta) if trace['eta_min'] is None else min(trace['eta_min'], float(eta))
        trace['eta_max'] = float(eta) if trace['eta_max'] is None else max(trace['eta_max'], float(eta))
        if trace['first_extreme_bb_call'] is None and abs(eta) > 100:
            trace['first_extreme_bb_call'] = dict(phase=phase, eta=float(eta), linear_r=float(linear_r),
                                                n=int(n), stack=traceback.format_stack())
            save('numerical_trace', trace)
        try:
            lp, grad, hess = stable_bb_terms(k, n, eta, linear_r)
        except Exception:
            if trace['first_exception'] is None:
                trace['first_exception'] = dict(phase=phase, callsite='stable_bb_terms', traceback=traceback.format_exc())
                save('numerical_trace', trace)
            raise
        return lp, grad, hess if need_hessian else None
    likelihood.bb_terms = stable_wrapper
    loader = importlib.import_module('bfix_v1')
    require(loader.FIXTURE == FIXTURE, 'Loader selected another fixture')
    try:
        design, settings = loader.load()
        model = design.model()
        require(settings['B'] == 50, 'Original draw contract changed')
        # At the historical observed parameters, compare entire LL/score/Hessian.
        original_model = original_module.Model(model.rows)
        equivalence = {}
        for name in ('observed_null', 'observed_marginal'):
            p = np.asarray(reference[name]['parameters'])
            old = original_model.evaluate(p)
            new = model.evaluate(p)
            differences = [float(np.max(np.abs(a-b))) for a, b in zip(old, new)]
            require(all(d < 1e-6 for d in differences), 'Fixed-parameter likelihood/derivative mismatch')
            equivalence[name+'_fixed_parameters'] = dict(zip(('loglik', 'gradient', 'hessian'), differences))
            if name == 'observed_marginal':
                old_marginal = original_model.marginal(p, reference[name]['v'])
                new_marginal = model.marginal(p, reference[name]['v'])
                marginal_differences = [float(np.max(np.abs(a-b))) for a, b in zip(old_marginal, new_marginal)]
                require(all(d < 1e-6 for d in marginal_differences), 'Fixed marginal likelihood/score mismatch')
                equivalence[name+'_fixed_marginal'] = dict(zip(('loglik', 'gradient'), marginal_differences))
        save('fixed_parameter_equivalence', equivalence)

        # Exactly the observed-fit prefix of the supplied runner.
        phase = 'observed_null_fit'
        old_null = json.loads(PRIOR.read_text())['fit']
        require(old_null['converged'], 'Original warm start not converged')
        null = model.fit(np.asarray(old_null['parameters']))
        save('observed_null', null)
        require(null['converged'], 'Observed null did not converge')
        null_p = np.asarray(null['parameters'])
        phase = 'observed_marginal_fit'
        marginal = model.fit_marginal(settings['tau_max_fixture'], null_p, null_fit=null)
        save('observed_marginal', marginal)
        require(marginal['converged'], 'Observed marginal did not converge')
        for name, fitted in (('observed_null', null), ('observed_marginal', marginal)):
            prior = reference[name]
            differences = dict(loglik=abs(fitted['loglik']-prior['loglik']),
                               kappa=abs(fitted['parameters'][2*model.G]-prior['parameters'][2*model.G]))
            if name == 'observed_marginal':
                differences['v'] = abs(fitted['v']-prior['v'])
            equivalence[name+'_refit'] = differences
            require(all(d <= 1e-6 for d in differences.values()), 'Observed refit differs by more than1e-6')
        save('observed_equivalence', equivalence)

        ordered = [(r['individual_id'], r['gene_id']) for r in loader.table('rows.tsv')]
        require(ordered == [(r.individual, settings['gene_order'][r.gene]) for r in model.rows], 'Row order changed')
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy', mmap_mode='r', allow_pickle=False)
        require(uniforms.shape == (50, len(model.rows), 2), 'Supplied uniform shape differs')
        observed_lr = 0. if marginal['v'] == 0. else max(0., 2*(marginal['loglik']-null['loglik']))
        fit = None
        try:
            phase = 'draw0_conditional_generation'
            rows = []
            for row, (u1, u2) in zip(model.rows, uniforms[0]):
                eta, r, _, _ = model.predictors(row, null_p)
                rows.append(likelihood.conditional_draw(row, eta, r, u1, u2))
            save('generation', dict(completed=True, synthetic_rows=len(rows),
                                    count_vector_sha256=hashlib.sha256(np.asarray([r.a for r in rows], dtype='<i8').tobytes()).hexdigest()))
            phase = 'draw0_original_kept_model_fit'
            fit = likelihood.Model(rows).fit_marginal(settings['tau_max_fixture'], null_p, section4_from_null=True)
            okay = fit['converged'] and fit['null']['converged']
            lr = 0. if fit['v'] == 0. else max(0., 2*(fit['loglik']-fit['null']['loglik']))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            okay = False
            if trace['first_exception'] is None:
                trace['first_exception'] = dict(phase=phase, callsite='original_generation_or_fit', traceback=traceback.format_exc())
        if not okay:
            lr = float('nan')
        result = dict(b=0, fit=fit, LR=lr, exceed=not okay or lr >= observed_lr-1e-6)
        save('p3_000', result)
        save('numerical_trace', trace)
        save('summary', dict(replayed_draws=1, observed_LR=observed_lr, first_draw_okay=bool(okay),
                             original_failed_exceedance_retained=receipt['original_failed_exceedance'],
                             replay_result=result, observed_equivalence=equivalence,
                             elapsed_seconds=time.monotonic()-started,
                             note='Single supplied synthetic draw diagnostic; no MC p-value, size or power claim'))
    except Exception:
        if trace['first_exception'] is None:
            trace['first_exception'] = dict(phase=phase, callsite='replay_wrapper', traceback=traceback.format_exc())
        save('numerical_trace', trace)
        raise
    finally:
        warnings.showwarning = standard_warning
    require(all(sha(SNAPSHOT / name) == expected for name, expected in HASHES.items()), 'Original snapshot changed')
    require(all(sha(FIXTURE / name) == digest for name, digest in inputs.items()), 'Synthetic input changed')
    require(sha(probe) == PROBE_SHA, 'Stable arithmetic source changed during replay')


if __name__ == '__main__':
    main()
