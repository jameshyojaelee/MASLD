"""REC B-FIX P1/P2/P3 under spec v1.8.4 (addendum sha256 eaa1a3b891e4234b6c4af4de6ab0578e2184a5e5daa5d84d575e8629e9670978).

Copy of ../v184-20261017T210000Z/run_bfix_v184.py (sha256 18030631795c436e428354502b90ee439bb079e97efc0dfe5aa7ea794cb05407),
which is not edited. The P3 class-order draw is imported unchanged from ../v184-20261017T210000Z/draw_v184.py
(sha256 4eb101c0...3636). Changes from that runner:
  * ModelV184 (likelihood_v184.py): every section-3 refit (P1, P2 kappa, P3) from K = 3 starts (warm, cold,
    probe j = 1) with probe entropy [20260923, 13, t, b, 1] (t = 1 P1, 2 P2 kappa, 3 P3; b = draw, 0-based);
    selection item 3.
  * Observed section-3 fit from three starts: start 1 jitter [20260923, 7, 0] (unchanged), start 2 probe
    [20260923, 13, 0, 0, 1]; selection item 5.2.
  * P3 section 4 (observed and each draw) is the item 6.1 pair from the kept section-3 theta at v = 0 and
    v = tau_max^2 (B-FIX tau_max_fixture), the K = 2 table with eps_tie,4, l3*, and LR by item 6.4.
  * P2 tau profile fits at tau0 keep v1.7 rule 5 (unchanged code path); l_p(tau_hat) uses l3, not l3* (C11).
  * Every per-draw record carries full end points and l of every start and member (item 9 step 1 export).
P2 tau fallback (stream t = 6) is not part of the B-FIX P2 runner and is not run here.

Synthetic agreement fixture, not calibration. Reads only the supplied B-FIX inputs. The only random numbers
generated are the item 4 probe starts and the REC restart jitter.

Roles (each a separate SLURM job): observed | p1 --draws a:b | p2k --draws a:b | p2tau --source s | p3 --draws a:b
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
DRAWFIX = PARENT / 'v184-20261017T210000Z'
SOURCES = {'likelihood_v184.py': HERE, 'run_bfix_v184full.py': HERE, 'draw_v184.py': DRAWFIX,
           'likelihood_v183.py': PARENT / 'v183', 'likelihood_v18.py': PARENT,
           'likelihood_v1.py': PARENT, 'bfix_v1.py': PARENT, 'design_v1.py': PARENT, 'tests_v1.py': PARENT}
OBS_START1_ENTROPY = (20260923, 7, 0)
OBS_PROBE_ENTROPY = (20260923, 13, 0, 0, 1)
PROBE_BASE = (20260923, 13)
REC_V183_OBSERVED = Path('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/'
                         'MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/'
                         'rec-v1.8.3-20261016T234346Z/observed')
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
    from likelihood_v183 import observed_failure
    from likelihood_v184 import ModelV184

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

    receipt = {'role': args.role, 'spec': 'v1.8.3b (spec sha256 285c254cb65bf4bcc411c0f984c1ea13c2f23ba552d42c1c8b2e74894c883a9b + H-14) + v1.8.4 addendum sha256 eaa1a3b891e4234b6c4af4de6ab0578e2184a5e5daa5d84d575e8629e9670978',
               'family': 'REC B-FIX v1.8.4 full', 'header': header(), 'job': job, 'draws': args.draws,
               'source': args.source, 'source_sha256': hashes, 'fixture': str(FIXTURE),
               'fixture_sha256': inputs, 'python': platform.python_version(), 'numpy': np.__version__,
               'scipy': scipy.__version__, 'host': platform.node(),
               'rng': 'fixture-supplied draws; probe starts [20260923, 13, t, b, 1] and observed [20260923, 13, 0, 0, 1]; observed start 1 [20260923, 7, 0]; REC restart jitter [20260923, 90, a]',
               'note': 'Synthetic agreement fixture, not calibration.'}
    save('receipt', receipt)
    start = time.monotonic()
    design, settings = load()
    assert settings['B'] == 50
    model = ModelV184.from_model(design.model())
    tau_max = settings['tau_max_fixture']

    def as_v18(m):
        return ModelV184(m.rows, model.G)

    if args.role == 'observed':
        null = model.fit_observed3(OBS_START1_ENTROPY, OBS_PROBE_ENTROPY)
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
        if failure is not None:
            save('observed_summary', {'observed_fit_failed': True, 'observed_failure': failure,
                                      'elapsed_seconds': time.monotonic() - start})
            raise SystemExit('observed fit failed: ' + failure)
        pair = model.section4_pair(null['parameters'], null['loglik'], tau_max ** 2)
        kept = pair['kept_member']
        marginal = {'parameters': None if kept is None else pair['members'][kept]['end_point'],
                    'v': pair['v_kept'], 'loglik_reference': pair['l4_reference'],
                    'converged': not pair['failed'], 'pair': pair}
        if kept is not None:
            marginal['loglik'] = pair['members'][kept]['loglik']
        save('observed_marginal', {**marginal, 'stats': model.stats, 'events': model.events[:200]})
        z = float(null['parameters'][2 * model.G] / se)
        prev = json.loads((REC_V183_OBSERVED / 'observed_null.json').read_text())
        p_old = np.asarray(prev['parameters'], float); p_new = np.asarray(null['parameters'], float)
        save('observed_summary', {'observed_fit_failed': False, 'observed_z': z, 'observed_se': se,
                                  'observed_kappa': float(null['parameters'][2 * model.G]),
                                  'observed_loglik': null['loglik'], 'reported_start': null['reported_start'],
                                  'selection_rule': null['selection_rule'],
                                  'max_unconverged_minus_candidate': null['max_unconverged_minus_candidate'],
                                  'start_records': [{k: v for k, v in r.items() if k != 'attempt_records'}
                                                    for r in null['starts']],
                                  'observed_section4_failed': pair['failed'],
                                  'observed_P3_failed': pair['failed'],
                                  'observed_marginal_v': pair['v_kept'], 'observed_kept_member': kept,
                                  'observed_l3_star_minus_l3': pair['l3_star_minus_l3'],
                                  'observed_missed_v0_point': pair['missed_v0_point'],
                                  'tau_max_fixture': tau_max, 'observed_LR': pair['LR'],
                                  'members': [{k: v for k, v in r.items() if k != 'attempt_records'}
                                              for r in pair['members']],
                                  'rec_v183_comparison': {
                                      'v183_loglik': prev['loglik'], 'loglik_difference': null['loglik'] - prev['loglik'],
                                      'kappa_difference': float(p_new[2 * model.G] - p_old[2 * model.G]),
                                      'max_abs_parameter_difference': float(np.max(np.abs(p_new - p_old))),
                                      'parameters_bitwise_equal': bool(np.array_equal(p_new, p_old))},
                                  'elapsed_seconds': time.monotonic() - start})
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
                fit = m.fit_refit(null_p, PROBE_BASE + (1, k, 1))
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
                fit = m.fit_refit(null_p, PROBE_BASE + (2, k, 1))
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
        elif marginal['parameters'] is None or not marginal['converged']:
            record.update(p_tau=1., reason='observed section-4 pair failed')
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
            if fit is None or not fit['converged'] or not marginal['converged'] or marginal['parameters'] is None:
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
                m = ModelV184(rows, model.G)
                fit = m.fit_refit(null_p, PROBE_BASE + (3, k, 1))
                pair = None
                if fit['converged']:
                    pair = m.section4_pair(fit['parameters'], fit['loglik'], tau_max ** 2)
                okay = bool(fit['converged'] and pair is not None and not pair['failed'])
                lr = pair['LR'] if okay else float('nan')
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                okay = False; err = traceback.format_exc(); pair = None
            if not okay:
                lr = float('nan')
            save('p3_%03d' % k, {'b': k, 'LR': lr, 'converged': okay,
                                 'exceed': (not okay) or lr >= observed_lr - 1e-6,
                                 'section3_failed': None if fit is None else bool(not fit['converged']),
                                 'section4_failed': None if pair is None else pair['failed'],
                                 'v_hat': None if pair is None else pair['v_kept'],
                                 'v_on_upper_bound': None if pair is None or pair['v_kept'] is None
                                 else bool(abs(pair['v_kept'] - 1) <= 1e-8),
                                 'l3_star_minus_l3': None if pair is None else pair['l3_star_minus_l3'],
                                 'missed_v0_point': None if pair is None else pair['missed_v0_point'],
                                 'fit': fit, 'section4_pair': pair, 'error': err, 'seconds': time.monotonic() - t0,
                                 'class_changed_vs_v1_order': changed,
                                 'stats': None if m is None else m.stats,
                                 'events': None if m is None else m.events[:300]})
    for name in SOURCES:
        if sha(snap / name) != hashes[name]:
            raise RuntimeError('snapshot changed: ' + name)
    print('done', args.role, time.monotonic() - start, flush=True)


if __name__ == '__main__':
    main()
