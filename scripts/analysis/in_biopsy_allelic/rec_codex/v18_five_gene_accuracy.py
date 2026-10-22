"""Spec §4 accuracy tests under v1.8 §1 (REC). Synthetic only; no fits.

Points:
  --gene g (0..4): the §4 five-gene fixture (5/25/75/200/500 rows, construction
      and RNG order of check_quadrature_v1.py, seed 20260923) at
      tau = 0.05 x {0.25, 1, 4}.
  --case failed0gene2 | positive18: gene 2 of B-FIX P3 draw 0 at the v1.7
      failed point (v = 0.92348...) and of draw 18 at its v1.7 fit (v = 0.00712).

At each point it records l_g from: the v1.8 production evaluation (final-point
rules, whatever integrator v1.8 selects), GH20 and GH40 at the v1.8 mode, the
v1.8 reference integrator, the original v1 GH20 (u = 0 Newton mode search) and
an independent scalar QUADPACK integral written without conditional_terms.
Pass (v1.8 §1.1): |difference| <= 1e-8 in l_g (nats).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
import warnings

HERE = Path(__file__).resolve().parent
SOURCES = ('likelihood_v18.py', 'likelihood_v1.py', 'bfix_v1.py', 'design_v1.py', 'tests_v1.py',
           'v18_five_gene_accuracy.py')
REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/'
           'executions/codex-rec-20260929T142434Z')
FIXTURE = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/'
               'executions/model-a-bfix-20260930T115521Z')
COUNTS = (5, 25, 75, 200, 500)
MULTIPLES = (.25, 1., 4.)
TOL = 1e-8


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-dir', required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--gene', type=int, choices=range(5))
    g.add_argument('--case', choices=('failed0gene2', 'positive18'))
    args = ap.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    if not job.isdigit():
        raise SystemExit('run under SLURM only')
    tag = 'gene%d' % args.gene if args.case is None else args.case
    out = Path(args.run_dir) / ('accuracy_' + tag + '_' + job)
    out.mkdir(parents=True, exist_ok=False)
    snap = out / 'executed_sources'
    snap.mkdir()
    hashes = {}
    for name in SOURCES:
        shutil.copy2(HERE / name, snap / name)
        hashes[name] = sha(snap / name)
    (out / 'command.txt').write_text(' '.join([sys.executable] + sys.argv) + '\n')
    (out / 'pip_freeze.txt').write_text(subprocess.run([sys.executable, '-m', 'pip', 'freeze'],
                                                       capture_output=True, text=True).stdout)
    sys.path.insert(0, str(snap))
    import numpy as np
    import scipy
    from scipy.integrate import quad
    from scipy.special import betaln, expit, gammaln, logsumexp
    from scipy.stats import norm
    import likelihood_v1 as L
    from likelihood_v18 import ModelV18

    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    def save(name, value):
        with (out / (name + '.json')).open('x') as fh:
            fh.write(json.dumps(value, default=serial, indent=1) + '\n')
        print(name, 'saved', time.strftime('%FT%T'), flush=True)

    save('receipt', dict(spec='v1.8', job=job, point=tag, source_sha256=hashes, tolerance_nats=TOL,
                         numpy=np.__version__, scipy=scipy.__version__, python=sys.version,
                         seed='20260923 (five-gene fixture construction, as check_quadrature_v1.py)',
                         note='Synthetic accuracy test, no fits, no calibration.'))

    # ---------------------------------------------------------------- points
    points = []
    if args.case is None:
        rows = []
        for gene, count in enumerate(COUNTS):
            for j in range(count):
                rows.append(L.Row(gene, f'{j:04d}', 'fixture', 20 + j % 17, 10, .6, .3, .1, .01, .02, .05,
                                  (-1) ** (gene + j), float(j % 3 - 1), (float(j % 5 - 2) / 10,),
                                  float(j % 7 - 3) / 10, float(j % 5 - 2) / 10, 10, 2, 10))
        model = L.Model(rows)
        p = model.initial(); p[:5] = [-1.2, -.7, .4, .9, 1.3]; p[10] = .08; p[11] = .04
        rng = np.random.default_rng(20260923)
        drawn = []
        for row in model.rows:
            eta, r, _, _ = model.predictors(row, p)
            drawn.append(L.conditional_draw(row, eta, r, *rng.random(2)))
        G, Z = 5, 1
        rs = [r for r in drawn if r.gene == args.gene]
        for multiple in MULTIPLES:
            points.append(dict(label=f'gene{args.gene}_tau{multiple}', gene=args.gene, rows=rs, p=p,
                               v=(multiple * .05) ** 2, G=G))
    else:
        from bfix_v1 import load
        design, settings = load()
        base = design.model()
        null = json.loads((REC / 'p3_parallel_20261003T130323Z/baseline/observed_null.json').read_text())
        null_p = np.asarray(null['parameters'])
        draw = 0 if args.case == 'failed0gene2' else 18
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy')
        rows = []
        for row, (u1, u2) in zip(base.rows, uniforms[draw]):
            eta, r, _, _ = base.predictors(row, null_p)
            rows.append(L.conditional_draw(row, eta, r, u1, u2))
        if draw == 0:
            saved = json.loads((REC / 'p3_first_draw_diagnostic_21998696/summary.json').read_text())['latest_marginal']
        else:
            saved = json.loads((REC / 'p3_parallel_20261003T130323Z/baseline/p3_018.json').read_text())['fit']
        p, v = np.asarray(saved['parameters']), float(saved['v'])
        G, Z = base.G, base.Z
        rs = [r for r in L.Model(rows, gene_count=G).rows if r.gene == 2]
        assert len(rs) == 24
        points.append(dict(label=args.case, gene=2, rows=rs, p=p, v=v, G=G))

    # ----------------------------------------------------- independent value
    def direct(rs, p, G, Zn, u):
        total = 0.
        for row in rs:
            k0 = 2 * G
            eta = p[row.gene] * (1 + p[k0] * row.stage + np.dot(p[k0 + 1:k0 + 1 + Zn], row.z))
            eta += u * row.stage - row.orientation * p[k0 + 1 + Zn] * row.stage - row.orientation * row.omega
            rr = p[G + row.gene] + np.dot(p[-3:], (row.stage, row.background, row.duplication))
            m = expit(eta); rho = 1e-6 + (1 - 2e-6) * expit(rr)
            lo = max(row.minor, (row.percent * row.n + 99) // 100)
            k = np.arange(lo, row.n - lo + 1)
            prior = np.array([2 * row.p_ref * row.p_alt, row.p_ref ** 2 + 2 * row.p_ref * row.p_x,
                              row.p_alt ** 2 + 2 * row.p_alt * row.p_x])
            if row.orientation == -1:
                prior = prior[[0, 2, 1]]
            parts = []
            for mean, disp in ((m, rho), (row.e, row.phi), (1 - row.e, row.phi)):
                t = (1 - disp) / disp
                parts.append(gammaln(row.n + 1) - gammaln(k + 1) - gammaln(row.n - k + 1)
                             + betaln(k + mean * t, row.n - k + (1 - mean) * t) - betaln(mean * t, (1 - mean) * t))
            with np.errstate(divide='ignore'):
                mass = logsumexp(np.asarray(parts) + np.log(prior)[:, None], axis=0)
            total += mass[row.a - lo] - logsumexp(mass)
        return float(total)

    def scalar(rs, p, G, Zn, tau, bound):
        grid = np.linspace(-bound, bound, int(np.ceil(2 * bound / .25)) + 1)
        f = lambda z: direct(rs, p, G, Zn, tau * z) - z * z / 2 - .5 * np.log(2 * np.pi)
        shift = max(f(z) for z in grid)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            mass, error, info = quad(lambda z: np.exp(f(z) - shift), -bound, bound, points=grid[1:-1],
                                     epsabs=1e-12, epsrel=1e-12, limit=6000, full_output=1)[:3]
        tail = float(np.exp(np.log(2.) + norm.logsf(bound) - (shift + np.log(mass))))
        return dict(logmass=float(shift + np.log(mass)), z_bound=bound, evaluations=info['neval'],
                    quad_error_logmass=float(error / mass), prior_tail_relative=tail,
                    warnings=[str(w.message) for w in caught])

    # ----------------------------------------------------------- evaluate
    summary = []
    for pt in points:
        t0 = time.monotonic()
        rec = dict(label=pt['label'], gene=pt['gene'], rows=len(pt['rows']), v=pt['v'], tau=pt['v'] ** .5)
        rs, p, v, G = pt['rows'], pt['p'], pt['v'], pt['G']
        Zn = len(rs[0].z)
        m = ModelV18(rs, gene_count=G)
        g = pt['gene']
        try:
            rec['production'] = m.marginal(p, v, final=True, need_gradient=False)
            rec['production_integrator'] = m.last_choices.get(g)
            pre = m._gene_pre(p, g)
            m._seen.add(g)
            search = m._mode_search(g, pre, v, final=True)
            rec['mode_search'] = dict(U=search['U'], modes=[dict(u=x['u'], h=x['h'], d2=x['d2']) for x in search['modes']],
                                      best_u=search['best']['u'], triggers=search['triggers'], grid=search['grid'],
                                      multimodal=search['multimodal'], flat=search['flat'])
            if m._sigma(search['best']) is not None:
                for nodes in (20, 40):
                    u, logw = m._gh_nodes(search['best'], nodes)
                    rec['v18_GH%d' % nodes] = m._node_integral(p, g, pre, v, u, logw, False, False)[0]
            ref = m._reference_nodes(pre, v, search)
            rec['reference'] = ref[0]; rec['reference_panels'] = ref[3]; rec['reference_W'] = ref[4]
        except Exception:
            rec['v18_error'] = traceback.format_exc()
        try:
            rec['original_v1_GH20'] = L.Model(rs, gene_count=G).marginal(p, v, nodes=20)[0]
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            rec['original_v1_GH20_refused'] = str(exc)
        tau = v ** .5
        bound = 12.
        for _ in range(11):
            first = scalar(rs, p, G, Zn, tau, bound)
            if first['prior_tail_relative'] <= 1e-10:
                break
            bound += 8
        expanded = scalar(rs, p, G, Zn, tau, bound + 2)
        rec['scalar'] = dict(first=first, expanded=expanded,
                             refinement_difference=expanded['logmass'] - first['logmass'])
        S = expanded['logmass']
        for key in ('production', 'v18_GH20', 'v18_GH40', 'reference', 'original_v1_GH20'):
            if key in rec:
                rec[key + '_minus_scalar'] = rec[key] - S
                rec[key + '_scalar_pass'] = bool(abs(rec[key] - S) <= TOL)
        if 'v18_GH20' in rec and 'v18_GH40' in rec:
            rec['GH20_minus_GH40'] = rec['v18_GH20'] - rec['v18_GH40']
            rec['GH20_GH40_pass'] = bool(abs(rec['GH20_minus_GH40']) <= TOL)
        rec['scalar_resolved'] = bool(abs(rec['scalar']['refinement_difference']) + expanded['quad_error_logmass'] <= 1e-9)
        rec['seconds'] = time.monotonic() - t0
        rec['stats'] = m.stats
        save('point_' + pt['label'], rec)
        summary.append({k: rec.get(k) for k in ('label', 'rows', 'v', 'production_integrator',
                                                 'production_scalar_pass', 'v18_GH20_scalar_pass',
                                                 'GH20_GH40_pass', 'reference_scalar_pass',
                                                 'original_v1_GH20_scalar_pass', 'scalar_resolved',
                                                 'production_minus_scalar', 'v18_GH20_minus_scalar',
                                                 'GH20_minus_GH40', 'reference_minus_scalar')})
    save('summary', dict(points=summary, tolerance_nats=TOL))


if __name__ == '__main__':
    main()
