"""REC B-FIX P1/P2/P3, family 'v1.8.4' = v1.8.3 + P3 class-order fix, 50 supplied draws.

Copy of ../v183/run_bfix_v183.py (sha256 c9dee58a2cfd6990f785e7369010958d52983605cfd182bd61074095c1e85ebf).
Only change: P3 draws use draw_v184.conditional_draw_v184, which walks genotype classes in the
order (H, R, A) for every row as B-FIX p3_uniform_use states; likelihood_v1.conditional_draw used
index order, which is (H, A, R) for orientation -1 rows. P1, P2 and the observed fit do not call the
draw and are unchanged. Each P3 record also stores which rows changed class versus the v1 order.

Based on ../run_bfix_v18.py (unchanged). Changes: ModelV183; observed section-3 fit from two
starts (v1.8.3 item 3, B-FIX start-1 entropy [20260923, 7, 0]); failed observed fit per item 4;
job header with host, CPU model, glibc, node features, thread variables, GPU and versions.

Synthetic agreement fixture, not calibration. Reads only the supplied B-FIX
inputs (globals.json, individuals.tsv, rows.tsv, cohorts.tsv, draws/*); never
make_fixture.py or rows_true_class.tsv. No random numbers are generated: every
draw is supplied by the fixture (P1 S*, P2 resample + lambda, P3 uniforms).

Roles (each a separate SLURM job):
  observed            observed section-3 and section-4 fits -> <run>/observed/
  p1  --draws a:b     P1 CRT refits for draws a..b-1
  p2k --draws a:b     P2 kappa bootstrap refits for draws a..b-1
  p2tau --source s    P2 tau profile for lambda source s (0 GSE193066, 1 PXD051911_A1)
  p3  --draws a:b     P3 bootstrap refits for draws a..b-1
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
SOURCES = {'draw_v184.py': HERE, 'run_bfix_v184.py': HERE, 'likelihood_v183.py': PARENT / 'v183', 'likelihood_v18.py': PARENT,
           'likelihood_v1.py': PARENT, 'bfix_v1.py': PARENT, 'design_v1.py': PARENT, 'tests_v1.py': PARENT}
OBS_START1_ENTROPY = (20260923, 7, 0)
REC_V18_OBSERVED = Path('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/'
                        'MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/'
                        'rec-v1.8-20261015T234851Z/observed')
FIXTURE = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/'
               'executions/model-a-bfix-20260930T115521Z')
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv', 'draws/lambda_draws.tsv',
          'draws/p2_resample.tsv', 'draws/p3_uniforms.npy', 'draws/p1_Sstar.tsv')
V17 = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/'
           'executions/codex-rec-20260929T142434Z')
V17_PRIOR_NULL = V17 / 'bfix_fits_21984981/null.json'
V17_OBSERVED = V17 / 'bfix_tests_21985999'
SOURCE_NAMES = ('GSE193066', 'PXD051911_A1')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def header():
    def run(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=30).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None
    cpu = None
    try:
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if line.startswith('model name'):
                cpu = line.split(':', 1)[1].strip()
                break
    except OSError:
        pass
    host = platform.node()
    node = run(['scontrol', 'show', 'node', host.split('.')[0]]) or ''
    features = [t.split('=', 1)[1] for t in node.split() if t.startswith(('AvailableFeatures=', 'ActiveFeatures='))]
    return {'host': host, 'cpu_model': cpu, 'glibc': os.confstr('CS_GNU_LIBC_VERSION'),
            'node_features': features, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'threads': {k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                                                       'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
            'gpu': run(['nvidia-smi', '-L']) or 'none (CPU job)', 'python_executable': sys.executable,
            'conda_prefix': sys.prefix, 'kernel': platform.release()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--role', required=True, choices=('observed', 'p1', 'p2k', 'p2tau', 'p3'))
    ap.add_argument('--run-dir', required=True)
    ap.add_argument('--draws', default=None, help='a:b half-open range')
    ap.add_argument('--source', type=int, choices=(0, 1))
    args = ap.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    if not job.isdigit():
        raise SystemExit('run under SLURM only')
    run = Path(args.run_dir)
    tag = args.role + ('' if args.draws is None else '_' + args.draws.replace(':', '-')) + \
        ('' if args.source is None else '_s%d' % args.source)
    out = run / ('observed' if args.role == 'observed' else tag + '_' + job)
    out.mkdir(parents=True, exist_ok=False)
    snap = out / 'executed_sources'
    snap.mkdir()
    hashes = {}
    for name, where in SOURCES.items():
        shutil.copy2(where / name, snap / name)
        hashes[name] = sha(snap / name)
        assert hashes[name] == sha(where / name)
    published = dict(reversed(line.split(maxsplit=1)) for line in (FIXTURE / 'SHA256SUMS').read_text().splitlines())
    published = {k.removeprefix('./'): v for k, v in published.items()}
    inputs = {n: sha(FIXTURE / n) for n in INPUTS}
    for n, d in inputs.items():
        if n in published and published[n] != d:
            raise ValueError('fixture input checksum differs: ' + n)
    (out / 'command.txt').write_text(' '.join([sys.executable] + sys.argv) + '\n')
    (out / 'pip_freeze.txt').write_text(subprocess.run([sys.executable, '-m', 'pip', 'freeze'],
                                                       capture_output=True, text=True).stdout)
    sys.path.insert(0, str(snap))
    import numpy as np
    import scipy
    from scipy.stats import norm
    from bfix_v1 import load, table
    from draw_v184 import conditional_draw_v184
    from likelihood_v183 import ModelV183, observed_failure

    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, set):
            return sorted(value)
        raise TypeError(type(value).__name__)

    def save(name, value):
        with (out / (name + '.json')).open('x') as fh:
            fh.write(json.dumps(value, default=serial, indent=1) + '\n')
        print(name, 'saved', time.strftime('%FT%T'), flush=True)

    receipt = {'role': args.role, 'spec': 'v1.8.2 + v1.8.3 (spec sha256 285c254cb65bf4bcc411c0f984c1ea13c2f23ba552d42c1c8b2e74894c883a9b) + v1.8.4 P3 class order (H, R, A) per B-FIX p3_uniform_use',
               'family': 'REC B-FIX v1.8.4', 'header': header(), 'job': job, 'draws': args.draws,
               'source': args.source, 'source_sha256': hashes, 'fixture': str(FIXTURE),
               'fixture_sha256': inputs, 'python': platform.python_version(), 'numpy': np.__version__,
               'scipy': scipy.__version__, 'host': platform.node(),
               'rng': 'none generated; all draws supplied by the fixture',
               'note': 'Synthetic agreement fixture, not calibration.'}
    save('receipt', receipt)
    start = time.monotonic()
    design, settings = load()
    assert settings['B'] == 50
    model = ModelV183.from_model(design.model())
    tau_max = settings['tau_max_fixture']

    def as_v18(m):
        return ModelV183(m.rows, model.G)

    if args.role == 'observed':
        null = model.fit_observed(OBS_START1_ENTROPY)
        se = None
        if not null['observed_fit_failed']:
            try:
                se = model.sandwich(np.asarray(null['parameters']), design.clusters)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                se = float('nan')
        failure = observed_failure(null, se)
        null['observed_failure'] = failure
        null['converged'] = failure is None
        save('observed_null', {**null, 'observed_se': se, 'stats': dict(model.stats)})
        prev = {n: json.loads((REC_V18_OBSERVED / (n + '.json')).read_text())
                for n in ('observed_null', 'observed_marginal')}
        if failure is not None:
            save('observed_summary', {'observed_fit_failed': True, 'observed_failure': failure,
                                      'elapsed_seconds': time.monotonic() - start})
            raise SystemExit('observed fit failed: ' + failure)
        marginal = model.fit_marginal(tau_max, np.asarray(null['parameters']), null_fit=null)
        marginal.pop('null')
        save('observed_marginal', {**marginal, 'stats': model.stats, 'events': model.events[:200]})
        z = float(null['parameters'][2 * model.G] / se)
        comparison = {}
        for name, fitted in (('observed_null', null), ('observed_marginal', marginal)):
            p_old = np.asarray(prev[name]['parameters'], float)
            p_new = np.asarray(fitted['parameters'], float)
            comparison[name] = {'rec_v18_loglik': prev[name]['loglik'], 'v183_loglik': fitted['loglik'],
                                'loglik_difference': fitted['loglik'] - prev[name]['loglik'],
                                'kappa_difference': float(p_new[2 * model.G] - p_old[2 * model.G]),
                                'max_abs_parameter_difference': float(np.max(np.abs(p_new - p_old))),
                                'parameters_bitwise_equal': bool(np.array_equal(p_new, p_old))}
            if name == 'observed_marginal':
                comparison[name].update(rec_v18_v=prev[name]['v'], v183_v=fitted['v'])
        save('observed_summary', {'observed_fit_failed': False, 'observed_z': z, 'observed_se': se,
                                  'observed_kappa': float(null['parameters'][2 * model.G]),
                                  'observed_loglik': null['loglik'], 'reported_start': null['reported_start'],
                                  'selection_rule': null['selection_rule'],
                                  'max_unconverged_minus_candidate': null['max_unconverged_minus_candidate'],
                                  'start_records': [{k: v for k, v in r.items() if k != 'attempt_records'}
                                                    for r in null['starts']],
                                  'observed_marginal_v': marginal['v'], 'observed_marginal_converged': marginal['converged'],
                                  'tau_max_fixture': tau_max,
                                  'observed_LR': 0. if marginal['v'] == 0. else
                                  max(0., 2 * (marginal['loglik'] - null['loglik'])),
                                  'rec_v18_comparison': comparison, 'elapsed_seconds': time.monotonic() - start})
        return

    obs = run / 'observed'
    null = json.loads((obs / 'observed_null.json').read_text())
    marginal = json.loads((obs / 'observed_marginal.json').read_text())
    summary = json.loads((obs / 'observed_summary.json').read_text())
    null_p = np.asarray(null['parameters'])
    a, b = (map(int, args.draws.split(':')) if args.draws else (0, 0))

    if args.role == 'p1':
        stage_rows = table('draws/p1_Sstar.tsv')
        for k in range(a, b):
            row = stage_rows[k]
            assert int(row['draw']) == k
            sample = design.with_stage([int(row[p.individual]) for p in design.participants])
            t0 = time.monotonic(); fit = None; err = None
            try:
                m = as_v18(sample.model())
                fit = m.fit(null_p)
                sd = m.sandwich(fit['parameters'], sample.clusters)
                z = float(fit['parameters'][2 * model.G] / sd)
                okay = bool(fit['converged'] and np.isfinite(z))
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                okay = False; z = float('nan'); err = traceback.format_exc()
            save('p1_%03d' % k, {'b': k, 'z': z, 'converged': okay,
                                 'exceed': (not okay) or abs(z) >= abs(summary['observed_z']) - 1e-6,
                                 'fit': fit, 'error': err, 'seconds': time.monotonic() - t0,
                                 'stats': m.stats if fit is not None else None})

    elif args.role == 'p2k':
        groups = design.cohort_clusters()
        gi = {c: {design.clusters[mem[0]]: j for j, mem in enumerate(v)} for c, v in groups.items()}
        idx = {k: {c: [] for c in groups} for k in range(settings['B'])}
        for r in table('draws/p2_resample.tsv'):
            k = int(r['draw']); c = r['stratum']
            assert int(r['slot']) == len(idx[k][c])
            idx[k][c].append(gi[c][r['cluster_id']])
        lam = {(int(r['draw']), r['source']): float(r['lambda']) for r in table('draws/lambda_draws.tsv')}
        for k in range(a, b):
            sample = design.resample(idx[k])
            t0 = time.monotonic(); fit = None; err = None
            try:
                m = as_v18(sample.model())
                fit = m.fit(null_p)
                okay = bool(fit['converged']); kappa = float(fit['parameters'][2 * model.G])
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                okay = False; err = traceback.format_exc()
            if not okay:
                kappa = float('nan')
            save('p2_%03d' % k, {'b': k, 'kappa': kappa, 'converged': okay,
                                 'lambda': [lam[k, s] for s in SOURCE_NAMES], 'fit': fit, 'error': err,
                                 'seconds': time.monotonic() - t0, 'stats': m.stats if fit is not None else None})

    elif args.role == 'p2tau':
        lam = {(int(r['draw']), r['source']): float(r['lambda']) for r in table('draws/lambda_draws.tsv')}
        col = np.asarray([lam[k, SOURCE_NAMES[args.source]] for k in range(settings['B'])])
        cleaned = np.where(np.isfinite(col) & (col > .05), col, 0.)
        fifth = float(np.quantile(cleaned, .05))
        record = {'source': SOURCE_NAMES[args.source], 'lambda_fifth': fifth, 'tau_max': tau_max}
        if fifth <= .05:
            record.update(p_tau=1., reason='lambda fifth percentile <= 0.05')
        else:
            tau0 = tau_max * fifth
            mp = np.asarray(marginal['parameters'])
            err = None
            try:
                fit = model.fit_marginal(tau_max, mp, fixed_v=tau0 ** 2, null_fit=null, profile_start=mp)
                fit.pop('null')
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                fit = None; err = traceback.format_exc()
            record.update(tau0=tau0, profile_fit=fit, error=err)
            if fit is None or not fit['converged'] or not marginal['converged']:
                record.update(p_tau=1., reason='profile fit failed')
            else:
                lp_tau0 = model.reference_loglik(np.asarray(fit['parameters']), tau0 ** 2)
                lp_hat = model.reference_loglik(mp, marginal['v'])
                maximum = max(lp_hat, null['loglik'])
                lr = max(0., 2 * (maximum - lp_tau0))
                root = float(np.sign(tau0 - np.sqrt(marginal['v'])) * np.sqrt(lr))
                record.update(lp_tau0_reference=lp_tau0, lp_tauhat_reference=lp_hat, LR=lr,
                              signed_root=root, p_tau=float(norm.sf(root)),
                              profile_on_v_upper_bound=bool(abs(tau0 ** 2 - 1) <= 1e-8))
        record.update(stats=model.stats, events=model.events[:200], seconds=time.monotonic() - start)
        save('p2tau_s%d' % args.source, record)

    elif args.role == 'p3':
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy')
        ordered = [(r['individual_id'], r['gene_id']) for r in table('rows.tsv')]
        assert ordered == [(r.individual, settings['gene_order'][r.gene]) for r in model.rows]
        assert uniforms.shape == (settings['B'], len(model.rows), 2)
        observed_lr = summary['observed_LR']
        for k in range(a, b):
            t0 = time.monotonic(); fit = None; err = None
            m = None; changed = None
            try:
                rows = []; changed = []
                for i, (row, (u1, u2)) in enumerate(zip(model.rows, uniforms[k])):
                    eta, r, _, _ = model.predictors(row, null_p)
                    drawn, g, g_v1 = conditional_draw_v184(row, eta, r, u1, u2)
                    rows.append(drawn)
                    if g != g_v1:
                        label = 'HRA' if row.orientation == 1 else 'HAR'
                        changed.append({'row': i, 'orientation': int(row.orientation),
                                        'class_v184': label[g], 'class_v1_order': label[g_v1]})
                m = ModelV183(rows, model.G)
                fit = m.fit_marginal(tau_max, null_p, section4_from_null=True)
                okay = bool(fit['converged'] and fit['null']['converged'])
                lr = 0. if fit['v'] == 0. else max(0., 2 * (fit['loglik'] - fit['null']['loglik']))
                lr_ref = 0. if fit['v'] == 0. else max(0., 2 * (fit['loglik_reference'] - fit['null']['loglik']))
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                okay = False; err = traceback.format_exc()
            if not okay:
                lr = lr_ref = float('nan')
            save('p3_%03d' % k, {'b': k, 'LR': lr, 'LR_reference_endpoint': lr_ref, 'converged': okay,
                                 'exceed': (not okay) or lr >= observed_lr - 1e-6,
                                 'v_hat': None if fit is None else fit['v'],
                                 'v_on_upper_bound': None if fit is None else bool(abs(fit['v'] - 1) <= 1e-8),
                                 'fit': fit, 'error': err, 'seconds': time.monotonic() - t0,
                                 'class_changed_vs_v1_order': changed,
                                 'stats': None if m is None else m.stats,
                                 'events': None if m is None else m.events[:300]})
    for name in SOURCES:
        if sha(snap / name) != hashes[name]:
            raise RuntimeError('snapshot changed: ' + name)
    print('done', args.role, time.monotonic() - start, flush=True)


if __name__ == '__main__':
    main()
