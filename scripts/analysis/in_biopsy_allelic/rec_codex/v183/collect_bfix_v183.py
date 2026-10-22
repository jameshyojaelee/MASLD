"""Collect the REC B-FIX v1.8.3 family from one rec-v1.8.3-* run directory (copy of ../collect_bfix_v18.py).

Light (reads per-draw JSON, arithmetic on 50 values per test). Writes
<run>/collect-<UTC>/summary.json. p-values use the unchanged v1 formulas
(tests_v1.mc_result, tests_v1.p2_kappa) with failed refits as exceedances.
"""
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', required=True)
    args = ap.parse_args()
    import numpy as np
    from tests_v1 import mc_result, p2_kappa
    run = Path(args.run_dir)

    def draws(prefix):
        found = {}
        for path in sorted(run.glob('*/' + prefix + '_[0-9][0-9][0-9].json')):
            rec = json.loads(path.read_text())
            found[rec['b']] = (rec, str(path))
        return found

    out = {'run': str(run), 'spec': 'v1.8.2 + v1.8.3', 'note': 'Synthetic agreement fixture, not calibration.'}
    summary = json.loads((run / 'observed/observed_summary.json').read_text())
    out['observed'] = summary
    if summary.get('observed_fit_failed'):
        out['note_observed'] = 'observed fit failed (v1.8.3 item 4); P1 not computed'
    B = 50
    p1 = draws('p1')
    if len(p1) == B and not summary.get('observed_fit_failed'):
        exceed = [p1[b][0]['exceed'] for b in range(B)]
        out['P1'] = {**mc_result(exceed), 'failed': [b for b in range(B) if not p1[b][0]['converged']],
                     'draw_z': [p1[b][0]['z'] for b in range(B)]}
    else:
        out['P1'] = {'missing_draws': sorted(set(range(B)) - set(p1))}
    p2 = draws('p2')
    taus = {}
    for s in (0, 1):
        hits = sorted(run.glob('*/p2tau_s%d.json' % s))
        if hits:
            taus[s] = json.loads(hits[-1].read_text())
    if len(p2) == B and len(taus) == 2:
        kappas = [p2[b][0]['kappa'] for b in range(B)]
        lam = np.asarray([p2[b][0]['lambda'] for b in range(B)])
        kr = p2_kappa(kappas, lam, .125)
        pt = max(taus[s]['p_tau'] for s in (0, 1))
        out['P2'] = {'p': max(kr['p_kappa'], pt), 'kappa': kr, 'tau': {str(s): {k: taus[s].get(k) for k in (
            'source', 'lambda_fifth', 'tau0', 'LR', 'signed_root', 'p_tau', 'reason')} for s in taus},
            'failed_kappa_draws': [b for b in range(B) if not p2[b][0]['converged']],
            'failed_tau_profiles': [s for s in taus if taus[s].get('reason') == 'profile fit failed']}
    else:
        out['P2'] = {'missing_draws': sorted(set(range(B)) - set(p2)), 'tau_sources_present': sorted(taus)}
    p3 = draws('p3')
    if len(p3) == B:
        exceed = [p3[b][0]['exceed'] for b in range(B)]
        out['P3'] = {**mc_result(exceed), 'observed_LR': summary['observed_LR'],
                     'failed': [b for b in range(B) if not p3[b][0]['converged']],
                     'draw_LR': [p3[b][0]['LR'] for b in range(B)],
                     'v_hat': [p3[b][0]['v_hat'] for b in range(B)],
                     'v_on_upper_bound': [b for b in range(B) if p3[b][0]['v_on_upper_bound']],
                     'seconds': [p3[b][0]['seconds'] for b in range(B)]}
    else:
        out['P3'] = {'missing_draws': sorted(set(range(B)) - set(p3))}
    target = run / ('collect-' + time.strftime('%Y%m%dT%H%M%SZ', time.gmtime()))
    target.mkdir(exist_ok=False)
    (target / 'summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print(target / 'summary.json')


if __name__ == '__main__':
    main()
