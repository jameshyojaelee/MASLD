"""Pre-launch check of draw_v184 against likelihood_v1.conditional_draw on all 50 B-FIX P3 draws.

Uses the REC v1.8.3 observed null parameters (read only). Expectations written before running:
orientation +1 rows draw the same count under both orders; orientation -1 rows can differ only when
the class uniform falls where the (H, R, A) and (H, A, R) cumulative orders disagree; every drawn
class under v1.8.4 must equal the inverse CDF over (H, R, A) computed here a second, separate way.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / 'v183'), str(HERE.parent)]
FIXTURE = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/'
               'executions/model-a-bfix-20260930T115521Z')
V183_OBS = Path('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/'
                'executions/codex-rec-20260929T142434Z/rec-v1.8.3-20261016T234346Z/observed/observed_null.json')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'pip_freeze.txt').write_text(subprocess.run([sys.executable, '-m', 'pip', 'freeze'],
                                                       capture_output=True, text=True).stdout)
    import numpy as np
    from bfix_v1 import load
    from likelihood_v1 import conditional_draw
    from likelihood_v183 import ModelV183
    from draw_v184 import conditional_draw_v184, class_weights
    design, settings = load()
    model = ModelV183.from_model(design.model())
    null_p = np.asarray(json.loads(V183_OBS.read_text())['parameters'])
    uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy')
    per_draw = []; failures = []
    for b in range(settings['B']):
        changed = {1: 0, -1: 0}; count_changed = {1: 0, -1: 0}
        for i, (row, (u1, u2)) in enumerate(zip(model.rows, uniforms[b])):
            eta, r, _, _ = model.predictors(row, null_p)
            old = conditional_draw(row, eta, r, u1, u2)
            new, g, g_v1 = conditional_draw_v184(row, eta, r, u1, u2)
            _, _, w = class_weights(row, eta, r)
            # genotype-labelled weights (H, R, A): index 1 is R for +1 and A for -1
            w_hra = w if row.orientation == 1 else np.array([w[0], w[2], w[1]])
            cum = np.cumsum(w_hra)
            want = 'HRA'[min(2, int(np.sum(cum < u1)))]
            got = ('HRA' if row.orientation == 1 else 'HAR')[g]
            if want != got:
                failures.append(('class', b, i))
            if row.orientation == 1 and (old.a != new.a or g != g_v1):
                failures.append(('plus_row_changed', b, i))
            o = int(row.orientation)
            changed[o] += int(g != g_v1); count_changed[o] += int(old.a != new.a)
        per_draw.append({'b': b, 'class_changed_minus': changed[-1], 'class_changed_plus': changed[1],
                         'count_changed_minus': count_changed[-1], 'count_changed_plus': count_changed[1]})
    res = {'rows': len(model.rows), 'minus_rows': sum(r.orientation == -1 for r in model.rows),
           'draws': settings['B'], 'failures': failures,
           'draws_with_any_class_change': sum(d['class_changed_minus'] > 0 for d in per_draw),
           'total_class_changes': sum(d['class_changed_minus'] for d in per_draw), 'per_draw': per_draw}
    (out / 'check.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps({k: v for k, v in res.items() if k != 'per_draw'}))
    if failures:
        raise SystemExit('check failed')


if __name__ == '__main__':
    main()
