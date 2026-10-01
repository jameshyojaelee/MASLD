"""Compute-only local-coordinate derivative comparison on supplied B-FIX rows.

Widths64/379 add unobserved parameter slots, not participants or calibration
data. No fitting, starting-point optimization, integration or resampling.
"""
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time
import tracemalloc

import numpy as np
import scipy


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ORIGINAL = REC / 'bfix_p3_21998519'
SNAPSHOT = ORIGINAL / 'executed_sources'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HASHES = {
    'likelihood_v1.py': '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
    'bfix_v1.py': '324594ce783ba2464921eaee5924065f579060eded6043a48bdf4b7e392bfa56',
    'design_v1.py': 'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
    'tests_v1.py': '52cd14d43b649c20347330b64b467dbd81320f10252fe85ce3c538ceef2f092e',
    'check_bfix_tests_v13.py': 'a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d',
}
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sparse_model_class(likelihood):
    class SparseModel(likelihood.Model):
        def evaluate(self, p, rows=None, u=0, need_hessian=True):
            ll = 0.
            grad = np.zeros(self.size)
            hess = np.zeros((self.size, self.size)) if need_hessian else None
            shared = np.arange(2*self.G, self.size)
            for row in self.rows if rows is None else rows:
                g, G, Z = row.gene, self.G, self.Z
                k = 2*G
                factor = 1+p[k]*row.stage+np.dot(p[k+1:k+1+Z], row.z)
                eta = p[g]*factor-row.orientation*p[k+1+Z]*row.stage+u*row.stage
                linear_r = p[G+g]+np.dot(p[-3:], (row.stage, row.background, row.duplication))
                # [alpha_g, r_g, kappa, delta..., omega_S, rho_S,b,d]
                active = np.r_[g, G+g, shared]
                jac = np.zeros((2, Z+7))
                jac[0, 0] = factor
                jac[0, 2] = p[g]*row.stage
                jac[0, 3:3+Z] = p[g]*np.asarray(row.z)
                jac[0, 3+Z] = -row.orientation*row.stage
                jac[1, 1] = 1.
                jac[1, -3:] = (row.stage, row.background, row.duplication)
                value, score, hh = likelihood.conditional_terms(row, eta, linear_r, need_hessian)
                ll += value
                grad[active] += jac.T@score
                if need_hessian:
                    local_hess = jac.T@hh@jac
                    # eta is bilinear only in alpha_g with kappa/delta.
                    local_hess[0, 2] += score[0]*row.stage
                    local_hess[2, 0] += score[0]*row.stage
                    local_hess[0, 3:3+Z] += score[0]*np.asarray(row.z)
                    local_hess[3:3+Z, 0] += score[0]*np.asarray(row.z)
                    hess[np.ix_(active, active)] += local_hess
            return ll, grad, hess
    return SparseModel


def timed_evaluate(model, p, u, need_hessian, likelihood):
    # Clear the same pre-existing hom-inclusion cache for each evaluator.
    likelihood.hom_log_inclusion.cache_clear()
    tracemalloc.start()
    started = time.perf_counter()
    try:
        answer = model.evaluate(p, u=u, need_hessian=need_hessian)
        elapsed = time.perf_counter()-started
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return answer, dict(seconds=elapsed, traced_peak_bytes=peak)


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    output = REC / ('codex_sparse_gene_derivative_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'
    private.mkdir()
    def save(name, data):
        (output / (name+'.json')).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    for name, digest in HASHES.items():
        require(sha(SNAPSHOT/name) == digest, 'Immutable source changed: '+name)
        (private/name).write_bytes((SNAPSHOT/name).read_bytes())
    wrapper = Path(__file__).resolve()
    launcher = wrapper.with_name('run_codex_sparse_gene_derivative_probe.sbatch')
    (private/wrapper.name).write_bytes(wrapper.read_bytes())
    (private/launcher.name).write_bytes(launcher.read_bytes())
    original_receipt = json.loads((ORIGINAL/'receipt.json').read_text())
    inputs = {name: sha(FIXTURE/name) for name in INPUTS}
    published = {}
    for line in (FIXTURE/'SHA256SUMS').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        published[name.removeprefix('./')] = digest
    require(all(d == published[n] == original_receipt['fixture_sha256'][n]
                for n, d in inputs.items()), 'Original supplied synthetic input changed')
    observed_path = ORIGINAL/'observed_null.json'
    observed = json.loads(observed_path.read_text())
    require(observed['converged'], 'Archived observed synthetic null not converged')
    protocol = dict(widths=[8,64,379], seed=20260930, perturbation_sd=.03,
                    points=['archived_observed_null','moderate_1','moderate_2','moderate_3'],
                    u_values=[0.,-.05,.05], need_hessian_modes=[False,True],
                    derivative_atol=1e-8, derivative_rtol=1e-10,
                    source_sha256=HASHES, fixture_sha256=inputs,
                    observed_null_sha256=sha(observed_path), original_receipt_sha256=sha(ORIGINAL/'receipt.json'),
                    wrapper_sha256=sha(wrapper), launcher_sha256=sha(launcher),
                    synthetic_only=True, dimension_only_extra_slots=True,
                    fits_performed=False, integration_performed=False, law_changed=False,
                    timing='One evaluation per point and role, with tracemalloc; paired order alternates; hom cache cleared',
                    memory_limit_bytes=2*1024**3, time_limit_seconds=3600,
                    python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    blas_threads={n:os.environ.get(n) for n in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')})
    save('protocol', protocol)
    require(not any(n in sys.modules for n in ('bfix_v1','design_v1','likelihood_v1','tests_v1')),
            'Source isolation violated')
    sys.path.insert(0, str(private))
    likelihood = importlib.import_module('likelihood_v1')
    loader = importlib.import_module('bfix_v1')
    require(loader.FIXTURE == FIXTURE, 'Loader selected another fixture')
    design, _ = loader.load()
    base = design.model()
    require(base.G == 8 and len(base.rows) == 173, 'Expected supplied toy8-gene/173-row fixture')
    SparseModel = sparse_model_class(likelihood)
    base_p = np.asarray(observed['parameters'], float)
    require(base_p.shape == (base.size,), 'Archived parameter axis differs')
    rng = np.random.default_rng(protocol['seed'])
    parameters = [base_p]
    for _ in range(3):
        p = base_p+rng.normal(0., protocol['perturbation_sd'], base.size)
        for j, (lo, hi) in enumerate(base.bounds()):
            if lo is not None: p[j] = max(p[j],lo)
            if hi is not None: p[j] = min(p[j],hi)
        parameters.append(p)
    started = time.perf_counter()
    checks = []
    for width in protocol['widths']:
        dense = likelihood.Model(base.rows, gene_count=width)
        sparse = SparseModel(base.rows, gene_count=width)
        for point, small_p in enumerate(parameters):
            p = dense.initial()
            p[:base.G] = small_p[:base.G]
            p[width:width+base.G] = small_p[base.G:2*base.G]
            p[2*width:] = small_p[2*base.G:]
            for u in protocol['u_values']:
                for need_hessian in protocol['need_hessian_modes']:
                    record = dict(width=width, parameter_count=dense.size, point=protocol['points'][point],
                                  u=u, need_hessian=need_hessian,
                                  dense_row_eta_hessian_bytes=8*dense.size**2 if need_hessian else 0,
                                  local_row_hessian_bytes=8*(base.Z+7)**2 if need_hessian else 0,
                                  final_dense_hessian_bytes=8*dense.size**2 if need_hessian else 0)
                    answers = {}
                    roles = ('dense','local') if len(checks)%2 == 0 else ('local','dense')
                    record['execution_order'] = list(roles)
                    for role in roles:
                        answer, timing = timed_evaluate(dense if role=='dense' else sparse,p,u,need_hessian,likelihood)
                        answers[role] = answer
                        record[role] = timing
                    differences = {}
                    agreement = {}
                    for index, name in enumerate(('loglik','gradient','hessian')):
                        if name == 'hessian' and not need_hessian:
                            require(answers['dense'][index] is None and answers['local'][index] is None,
                                    'Score-only Hessian contract changed')
                            continue
                        before, after = answers['dense'][index], answers['local'][index]
                        require(np.isfinite(before).all() and np.isfinite(after).all(), 'Nonfinite synthetic derivatives')
                        differences[name] = float(np.max(np.abs(np.asarray(before)-after)))
                        agreement[name] = bool(np.allclose(before,after,atol=protocol['derivative_atol'],
                                                          rtol=protocol['derivative_rtol']))
                    if need_hessian:
                        agreement['hessian_symmetric'] = bool(np.allclose(
                            answers['local'][2], answers['local'][2].T,
                            atol=protocol['derivative_atol'],rtol=protocol['derivative_rtol']))
                    record['max_absolute_differences'] = differences
                    record['agreement'] = agreement
                    checks.append(record)
                    save('points', checks)
                    require(all(agreement.values()), 'Local derivative/symmetry guard failed; point saved')
                    require(time.perf_counter()-started < 3500, 'Fixed probe exceeded time buffer')
    for name, digest in HASHES.items():
        require(sha(SNAPSHOT/name) == digest, 'Original source changed during probe')
    require(all(sha(FIXTURE/name) == digest for name, digest in inputs.items()),
            'Supplied synthetic input changed during probe')
    require(sha(observed_path) == protocol['observed_null_sha256'],
            'Archived observed parameter point changed during probe')
    save('summary', dict(passed=True, comparisons=len(checks), elapsed_seconds=time.perf_counter()-started,
                         max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         limitation='Only173 supplied synthetic rows; wider axes have unobserved gene slots. No full379-gene workload, fit-runtime or calibration claim.'))


if __name__ == '__main__':
    main()
