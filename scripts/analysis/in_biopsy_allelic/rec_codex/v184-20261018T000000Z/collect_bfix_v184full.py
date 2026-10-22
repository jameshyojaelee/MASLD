"""Collect the REC B-FIX v1.8.4 full family from one rec-v1.8.4-full-* run directory.

Based on ../v184-20261017T210000Z/collect_bfix_v184.py (sha256 80fafc4d...f120, not edited). p-values use the
unchanged v1 formulas (tests_v1.mc_result, tests_v1.p2_kappa) with failed refits as exceedances. Adds, per
test: kept start index counts, failure lists by cause, and per-draw differences from the REC v1.8.3 run and
the REC v1.8.4 class-order run (P3). Writes <run>/collect-<UTC>/summary.json. Light.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
EX = Path('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/'
          'executions/codex-rec-20260929T142434Z')
PRIOR = {'v1.8.3': EX / 'rec-v1.8.3-20261016T234346Z', 'v1.8.4-class-order': EX / 'rec-v1.8.4-20261017T210000Z'}


def draws(run, prefix):
    found = {}
    for path in sorted(run.glob('*/' + prefix + '_[0-9][0-9][0-9].json')):
        rec = json.loads(path.read_text())
        found[rec['b']] = rec
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', required=True)
    args = ap.parse_args()
    import numpy as np
    from tests_v1 import mc_result, p2_kappa
    run = Path(args.run_dir)
    B = 50
    out = {'run': str(run), 'spec': 'v1.8.3b + v1.8.4 addendum (sha256 eaa1a3b8...0978)',
           'note': 'Synthetic agreement fixture, not calibration.'}
    summary = json.loads((run / 'observed/observed_summary.json').read_text())
    out['observed'] = summary

    def kept_counts(recs):
        return dict(Counter(str(recs[b]['fit']['kept_start']) if recs[b].get('fit') else 'no fit' for b in recs))

    def rules(recs):
        return dict(Counter(str(recs[b]['fit']['start_rule']) if recs[b].get('fit') else 'no fit' for b in recs))

    def probe_gain(recs):
        """Draws where the probe start ended more than 1e-6 above both warm and cold."""
        hits = []
        for b, r in recs.items():
            st = (r.get('fit') or {}).get('starts')
            if st and st[2]['loglik'] is not None and np.isfinite(st[2]['loglik']):
                if all(st[2]['loglik'] - s['loglik'] > 1e-6 for s in st[:2] if np.isfinite(s['loglik'])):
                    hits.append(b)
        return hits

    prior = {}
    for name, d in PRIOR.items():
        prior[name] = {t: draws(d, t) for t in ('p1', 'p2', 'p3')}

    p1 = draws(run, 'p1')
    if len(p1) == B and not summary.get('observed_fit_failed'):
        exceed = [p1[b]['exceed'] for b in range(B)]
        z = [p1[b]['z'] for b in range(B)]
        old = prior['v1.8.3']['p1']
        out['P1'] = {**mc_result(exceed), 'observed_z': summary['observed_z'],
                     'failed': [b for b in range(B) if not p1[b]['converged']], 'draw_z': z,
                     'kept_start_counts': kept_counts(p1), 'start_rules': rules(p1),
                     'probe_above_warm_and_cold': probe_gain(p1),
                     'max_abs_z_diff_vs_v183': max(abs(z[b] - old[b]['z']) for b in range(B)) if len(old) == B else None,
                     'indicator_changes_vs_v183': [b for b in range(B) if len(old) == B and old[b]['exceed'] != exceed[b]]}
    else:
        out['P1'] = {'missing_draws': sorted(set(range(B)) - set(p1))}

    p2 = draws(run, 'p2')
    taus = {}
    for s in (0, 1):
        hits = sorted(run.glob('*/p2tau_s%d.json' % s))
        if hits:
            taus[s] = json.loads(hits[-1].read_text())
    if len(p2) == B and len(taus) == 2:
        kappas = [p2[b]['kappa'] for b in range(B)]
        lam = np.asarray([p2[b]['lambda'] for b in range(B)])
        kr = p2_kappa(kappas, lam, .125)
        pt = max(taus[s]['p_tau'] for s in (0, 1))
        old = prior['v1.8.3']['p2']
        out['P2'] = {'p': max(kr['p_kappa'], pt), 'p_kappa': kr['p_kappa'], 'p_tau': pt, 'kappa': kr,
                     'tau': {str(s): {k: taus[s].get(k) for k in ('source', 'lambda_fifth', 'tau0', 'LR',
                                                                  'signed_root', 'p_tau', 'reason')} for s in taus},
                     'failed_kappa_draws': [b for b in range(B) if not p2[b]['converged']],
                     'failed_tau_profiles': [s for s in taus if taus[s].get('reason') not in (None,)],
                     'kept_start_counts': kept_counts(p2), 'start_rules': rules(p2),
                     'probe_above_warm_and_cold': probe_gain(p2),
                     'max_abs_kappa_diff_vs_v183': max(abs(kappas[b] - old[b]['kappa']) for b in range(B))
                     if len(old) == B else None}
    else:
        out['P2'] = {'missing_draws': sorted(set(range(B)) - set(p2)), 'tau_sources_present': sorted(taus)}

    p3 = draws(run, 'p3')
    if len(p3) == B:
        exceed = [p3[b]['exceed'] for b in range(B)]
        lr = [p3[b]['LR'] for b in range(B)]
        old = prior['v1.8.4-class-order']['p3']
        out['P3'] = {**mc_result(exceed), 'observed_LR': summary['observed_LR'],
                     'observed_P3_failed': summary.get('observed_P3_failed'),
                     'failed': [b for b in range(B) if not p3[b]['converged']],
                     'failed_section3': [b for b in range(B) if p3[b]['section3_failed']],
                     'failed_section4': [b for b in range(B) if p3[b]['section4_failed']],
                     'draw_LR': lr, 'v_hat': [p3[b]['v_hat'] for b in range(B)],
                     'kept_member_counts': dict(Counter(str((p3[b].get('section4_pair') or {}).get('kept_member'))
                                                        for b in range(B))),
                     'missed_v0_point': [b for b in range(B) if p3[b]['missed_v0_point']],
                     'l3_star_minus_l3_max': max((p3[b]['l3_star_minus_l3'] or 0.) for b in range(B)),
                     'v_on_upper_bound': [b for b in range(B) if p3[b]['v_on_upper_bound']],
                     'kept_start_counts': kept_counts(p3), 'start_rules': rules(p3),
                     'probe_above_warm_and_cold': probe_gain(p3),
                     'seconds': [p3[b]['seconds'] for b in range(B)],
                     'LR_diff_vs_v184_class_order': [lr[b] - old[b]['LR'] for b in range(B)] if len(old) == B else None,
                     'draws_LR_changed_vs_v184_class_order_gt_1e-9':
                         [b for b in range(B) if len(old) == B and abs(lr[b] - old[b]['LR']) > 1e-9]}
    else:
        out['P3'] = {'missing_draws': sorted(set(range(B)) - set(p3))}
    target = run / ('collect-' + time.strftime('%Y%m%dT%H%M%SZ', time.gmtime()))
    target.mkdir(exist_ok=False)
    (target / 'summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print(target / 'summary.json')


if __name__ == '__main__':
    main()
