#!/usr/bin/env python3
"""Independently reconstruct the supplied P3 counts; not B3 calibration."""
import argparse
from bisect import bisect_left
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
REFERENCE = REC/'repair_p3_draw11_decimal_likelihood_22087880/executed_checker.py'
NULL = REC/'p3_parallel_20261003T130323Z/baseline/observed_null.json'
JOBS = REC/'repair_p3_gamma_draw_v3_launch_20261007T034152Z/jobs.json'
FIXED = {REFERENCE: '38b3539bf30798b1ab34d738a45936ed00c62dc88e47170ac4644b1ea9a7b99b',
         NULL: 'a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181',
         JOBS: 'a7ca972d4a4c82916ba9b5bb8ef6aede1ecad652ce95a887d24fbb6653ab8cc9'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def metadata():
    for p, h in FIXED.items():
        require(sha(p) == h, 'Fixed input differs')
    spec = importlib.util.spec_from_file_location('independent_probability_reference', REFERENCE)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    pins = base.metadata()
    pins.update({str(p): h for p, h in FIXED.items()})
    jobs = json.loads(JOBS.read_text())
    require(len(jobs) == 50 and [j['draw'] for j in jobs] == list(range(50)), 'Draw roster differs')
    for j in jobs:
        p = Path(j['output_directory'])/'draw_counts.json'
        require(p == REC/f"repair_p3_gamma_draw_{j['draw']}_{j['jobid']}"/'draw_counts.json',
                'Count receipt escaped fixed synthetic family')
        pins[str(p)] = sha(p)
    return base, pins, jobs


def cumulative(mass):
    total = sum(mass, Decimal(0))
    require(total > 0, 'Empty modeled support')
    result = []
    acc = Decimal(0)
    for value in mass:
        acc += value/total
        result.append(acc)
    result[-1] = Decimal(1)
    return result


def reconstruct(base, rows, parameters, uniforms, digits):
    with localcontext() as ctx:
        ctx.prec = digits
        p = list(map(base.dec, parameters))
        counts = [[] for _ in range(50)]
        classes = [[] for _ in range(50)]
        margin = Decimal(1)
        normalization = Decimal(0)
        for ordinal, row in enumerate(rows):
            stage = base.dec(row.stage)
            factor = 1+p[16]*stage+sum((p[17+j]*base.dec(z) for j, z in enumerate(row.z)), Decimal(0))
            eta = p[row.gene]*factor-Decimal(row.orientation)*(p[24]*stage+base.dec(row.omega))
            r = p[8+row.gene]+sum((p[25+j]*base.dec(z) for j, z in enumerate(
                (row.stage, row.background, row.duplication))), Decimal(0))
            require(abs(eta) < 200 and abs(r) < 1000 and row.n <= 5000, 'Synthetic scope exceeded')
            eps = Decimal('0.000001')
            rho = eps+(1-2*eps)*base.logistic(r)
            supports = [base.bb_support(row.n, base.logistic(eta), rho),
                        base.bb_support(row.n, base.dec(row.e), base.dec(row.phi)),
                        base.bb_support(row.n, 1-base.dec(row.e), base.dec(row.phi))]
            normalization = max(normalization, *(e for _, e in supports))
            lo = max(row.minor, (row.percent*row.n+99)//100)
            hi = row.n-lo
            require(row.n >= row.depth and lo <= hi, 'Empty integer selection set')
            masses = [s[lo:hi+1] for s, _ in supports]
            pr, pa, px = map(base.dec, (row.p_ref, row.p_alt, row.p_x))
            weights = [2*pr*pa, pr*pr+2*pr*px, pa*pa+2*pa*px]
            if row.orientation == -1:
                weights[1], weights[2] = weights[2], weights[1]
            class_cdf = cumulative([w*sum(m, Decimal(0)) for w, m in zip(weights, masses)])
            count_cdfs = [cumulative(m) for m in masses]
            for draw in range(50):
                uc, un = map(base.dec, uniforms[draw, ordinal])
                require(0 <= uc < 1 and 0 <= un < 1, 'Uniform outside unit interval')
                genotype = bisect_left(class_cdf, uc)
                cdf = count_cdfs[genotype]
                k = bisect_left(cdf, un)
                # Record distance to either boundary surrounding each selected bin.
                for c, u, index in ((class_cdf, uc, genotype), (cdf, un, k)):
                    left = c[index-1] if index else Decimal(0)
                    margin = min(margin, abs(u-left), abs(c[index]-u))
                classes[draw].append(genotype)
                counts[draw].append(lo+k)
        require(normalization < Decimal(10)**(-(digits-15)), 'Probability normalization failed')
        return dict(digits=digits, counts=counts, genotype_classes=classes,
                    minimum_uniform_to_bin_boundary=str(margin),
                    maximum_support_normalization_error=str(normalization))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    base, pins, jobs = metadata()
    if args.inspect:
        print(json.dumps(dict(pins_verified=len(pins), draws=50, no_numeric_imports=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
            'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'One numerical thread required')
    out = REC/('p3_decimal_counts50_'+os.environ['SLURM_JOB_ID'])
    out.mkdir()
    (out/'executed_checker.py').write_bytes(Path(__file__).read_bytes())
    try:
        import numpy as np
        archive = REC/'repair_p3_gamma_draw_11_22081580/executed_sources'
        sys.path.insert(0, str(archive))
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1'):
            require(Path(sys.modules[name].__file__).resolve() == archive/(name+'.py'), 'Import escaped archive')
        design, settings = loader.load()
        rows = design.model().rows
        require(len(rows) == 173 and len(design.participants) == 60, 'Fixed dimensions differ')
        null = json.loads(NULL.read_text())
        require(null['converged'] and len(null['parameters']) == 28, 'Generating null differs')
        uniforms = np.load(FIXTURE/'draws/p3_uniforms.npy', allow_pickle=False)
        require(uniforms.shape == (50, 173, 2), 'Uniform dimensions differ')
        results = []
        for digits in (80, 110):
            result = reconstruct(base, rows, null['parameters'], uniforms, digits)
            (out/f'counts_{digits}.json').write_text(json.dumps(result, indent=2)+'\n')
            results.append(result)
        require(results[0]['counts'] == results[1]['counts'] and
                results[0]['genotype_classes'] == results[1]['genotype_classes'], 'Precision changes discrete draws')
        comparisons = []
        for j in jobs:
            saved = json.loads((Path(j['output_directory'])/'draw_counts.json').read_text())
            expected = results[1]['counts'][j['draw']]
            require(saved['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in rows],
                    'Saved row order differs')
            require(hashlib.sha256(np.asarray(saved['counts'], dtype='<i8').tobytes()).hexdigest() ==
                    saved['count_vector_sha256'], 'Saved count hash differs')
            differences = [i for i, (a, b) in enumerate(zip(expected, saved['counts'])) if a != b]
            require(len(saved['counts']) == 173, 'Saved count length differs')
            comparisons.append(dict(draw=j['draw'], differing_rows=differences))
        for p, h in pins.items():
            require(sha(Path(p)) == h, 'Input changed during reconstruction')
        summary = dict(source_and_input_sha256=pins, checker_sha256=sha(Path(__file__)),
                       all50_count_vectors_match=not any(x['differing_rows'] for x in comparisons),
                       comparisons=comparisons, precision_discrete_agreement=True,
                       no_candidate_count_or_probability_function_used=True,
                       no_fit_or_calibration=True, no_human_validation=True)
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        require(summary['all50_count_vectors_match'], 'Independent counts disagree')
    except BaseException:
        (out/'failure.json').write_text(json.dumps(dict(traceback=traceback.format_exc()), indent=2)+'\n')
        raise


if __name__ == '__main__':
    main()
