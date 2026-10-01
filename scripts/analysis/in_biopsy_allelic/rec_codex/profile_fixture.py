"""Measure synthetic likelihood cost without changing the live agreement run."""
import cProfile
import hashlib
import json
import os
import pstats
import time
from pathlib import Path
import numpy as np
from bfix_v1 import load


design, settings = load()
model = design.model()
root = Path(os.environ['CODEX_REC_OUTPUT'])
out = root/('fixture_profile_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
observed = json.loads((root/'bfix_tests_21985999/observed_null.json').read_text())
parameters = np.array(observed['parameters'])
profiler = cProfile.Profile()
start = time.monotonic()
profiler.enable()
for repeat in range(10):
    ll, score, hessian = model.evaluate(parameters)
    np.testing.assert_allclose(ll, observed['loglik'], atol=1e-8, rtol=0)
profiler.disable()
elapsed = time.monotonic()-start
profiler.dump_stats(str(out/'evaluate.prof'))
with (out/'profile.txt').open('w') as fh:
    pstats.Stats(profiler, stream=fh).strip_dirs().sort_stats('cumulative').print_stats(35)
summary = {'rows': len(model.rows), 'genes': model.G, 'repeats': 10,
           'elapsed_seconds': elapsed, 'seconds_per_evaluation': elapsed/10,
           'source_sha256': hashlib.sha256(Path('likelihood_v1.py').read_bytes()).hexdigest(),
           'limit': 'Synthetic section-3 likelihood evaluation only; not full fit or real-scale runtime.'}
(out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary), flush=True)
