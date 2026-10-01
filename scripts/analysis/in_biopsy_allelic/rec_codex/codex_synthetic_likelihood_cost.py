"""Fixed-point synthetic timing; no fit, integration change or live edit.

Compare repeated fixed homozygous-class calculations with bounded memoization.
This is a private numerical-cost diagnostic, not a production optimization.
"""
import cProfile
from functools import lru_cache
import hashlib
import importlib
import importlib.util
import io
import json
import os
from pathlib import Path
import platform
import pstats
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ORIGINAL = REC / 'bfix_p3_21998519'
SOURCE_HASHES = {
    'likelihood_v1.py':'41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
    'bfix_v1.py':'324594ce783ba2464921eaee5924065f579060eded6043a48bdf4b7e392bfa56',
    'design_v1.py':'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
}
PROBE_HASH = '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    import numpy as np
    import scipy
    output = REC / ('codex_synthetic_likelihood_cost_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'; private.mkdir()
    frozen = {}
    for name, expected in SOURCE_HASHES.items():
        source = ORIGINAL / 'executed_sources' / name
        require(sha(source) == expected, 'Frozen synthetic source changed')
        frozen[source] = expected
        (private / name).write_bytes(source.read_bytes())
    probe = Path(__file__).with_name('codex_stable_bb_probe.py')
    require(sha(probe) == PROBE_HASH, 'Stable numerical representation changed')
    frozen[probe] = PROBE_HASH
    (private / probe.name).write_bytes(probe.read_bytes())
    (private / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    sys.path.insert(0, str(private))
    likelihood = importlib.import_module('likelihood_v1')
    loader = importlib.import_module('bfix_v1')
    original_receipt = json.loads((ORIGINAL / 'receipt.json').read_text())
    for name, expected in original_receipt['fixture_sha256'].items():
        path = loader.FIXTURE / name
        require(sha(path) == expected, 'Original supplied fixture changed')
        frozen[path] = expected
    observed = ORIGINAL / 'observed_null.json'
    frozen[observed] = sha(observed)
    p = np.asarray(json.loads(observed.read_text())['parameters'], dtype=float)
    design, settings = loader.load(); model = design.model()
    require(len(model.rows) == 173 and settings['B'] == 50, 'Synthetic fixture identity differs')
    spec = importlib.util.spec_from_file_location('cost_stable_probe', private / probe.name)
    stable = importlib.util.module_from_spec(spec); spec.loader.exec_module(stable)
    original_bb, original_fixed = likelihood.bb_terms, likelihood.fixed_bb

    @lru_cache(maxsize=2048)
    def cached_fixed(keys, n, mean, rho):
        values = original_fixed(np.asarray(keys), n, mean, rho)
        values.setflags(write=False)
        return values

    def fixed_wrapper(k, n, mean, rho):
        return cached_fixed(tuple(np.asarray(k).tolist()), n, mean, rho)

    def stable_wrapper(k, n, eta, linear_r, need_hessian=True):
        lp, g, h = stable.stable_bb_terms(k, n, eta, linear_r)
        return lp, g, h if need_hessian else None

    def differences(left, right):
        return [None if a is None else float(np.max(np.abs(a-b)))
                for a, b in zip(left, right)]

    records = []
    for representation, bb in [('original', original_bb), ('stable', stable_wrapper)]:
        likelihood.bb_terms = bb
        cached_fixed.cache_clear()
        for hessian in (False, True):
            likelihood.fixed_bb = original_fixed
            reference = model.evaluate(p, need_hessian=hessian)
            likelihood.fixed_bb = fixed_wrapper
            cached = model.evaluate(p, need_hessian=hessian)
            error = differences(reference, cached)
            require(all(x is None or x == 0. for x in error), 'Memoization changed a numerical result')
            # Changing u must change only the heterozygous terms; fixed-class
            # cache reuse must remain exact at both prespecified displacements.
            for u in (-0.05, 0.05):
                likelihood.fixed_bb = original_fixed
                displaced = model.evaluate(p, u=u, need_hessian=hessian)
                likelihood.fixed_bb = fixed_wrapper
                require(all(x is None or x == 0. for x in differences(
                    displaced, model.evaluate(p, u=u, need_hessian=hessian))),
                    'Displaced synthetic point changed under memoization')
            times = {}
            for name, fixed in [('uncached', original_fixed), ('cached', fixed_wrapper)]:
                likelihood.fixed_bb = fixed
                samples = []
                for _ in range(5):
                    begin = time.perf_counter(); model.evaluate(p, need_hessian=hessian)
                    samples.append(time.perf_counter()-begin)
                times[name] = samples
            records.append(dict(representation=representation, need_hessian=hessian,
                exact_memoized_identity=True, displaced_u=[-0.05, 0.05],
                seconds=times, median_speed_ratio=float(np.median(times['uncached'])/np.median(times['cached'])),
                cache_info=cached_fixed.cache_info()._asdict()))
    profiler = cProfile.Profile()
    profiler.runcall(model.evaluate, p, need_hessian=False)
    text = io.StringIO()
    pstats.Stats(profiler, stream=text).sort_stats('cumtime').print_stats(25)
    (output / 'stable_cached_fixed_point_profile.txt').write_text(text.getvalue())
    require(all(sha(path) == expected for path, expected in frozen.items()), 'Input changed during check')
    summary = dict(fixed_point_records=records, synthetic_rows=173,
        fixed_input_sha256={str(path):value for path,value in frozen.items()},
        caching='Fixed homozygous BB classes only; exact scalar arguments/count tuple; max2048 entries',
        fits_run=0, marginal_integrations_run=0, statistical_rules_changed=False,
        live_source_changed=False, scientific_calibration=False,
        limits='Single synthetic design/observed parameter point and two u displacements; warm cache timings do not predict full refit or biological-design runtime.',
        job_id=os.environ['SLURM_JOB_ID'], python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        helper_sha256=sha(Path(__file__)))
    (output / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(passed=True, fits_run=0, output=str(output))), flush=True)


if __name__ == '__main__':
    main()
