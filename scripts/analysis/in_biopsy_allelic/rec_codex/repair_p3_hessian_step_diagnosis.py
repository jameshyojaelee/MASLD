#!/usr/bin/env python3
"""Retain every adjacent-step Richardson diagnostic; no accuracy acceptance.

This postdiagnostic numerical analysis uses saved synthetic fourth-order
differences. It does not select a step, change a scientific test,
rerun a probability model or certify any floating/integration error bound.
"""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
STEPS = (.001, .0003, .0001, .00003, .00001)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--input-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(args.input.resolve().is_relative_to(REC) and args.output.resolve().is_relative_to(REC),
            'Synthetic REC inputs/outputs only')
    body = args.input.read_bytes()
    require(len(body) < 64*1024**2 and hashlib.sha256(body).hexdigest() == args.input_sha,
            'Reviewed Hessian collection changed')
    collection = json.loads(body)
    groups = defaultdict(dict)
    for job in collection['jobs']:
        for record in job['candidate_comparisons']:
            expected = {'fourth_order_central': ([-2, -1, 1, 2], [1, -8, 8, -1]),
                        'fourth_order_forward_one_sided': ([0, 1, 2, 3, 4], [-25, 48, -36, 16, -3])}
            require(record['method'] in expected and record['completed'] is True
                    and (record['stencil']['offsets'], record['stencil']['coefficients'])
                    == expected[record['method']], 'Expected successful fourth-order stencil required')
            key = (job['case'], record['column'], record['score_axis'], record['complex_step'], record['bound'])
            require(record['real_step'] not in groups[key], 'Duplicate step in reference family')
            groups[key][record['real_step']] = record
    require(len(groups) == 588 and sum(map(len, groups.values())) == 2940
            and all(set(g) == set(STEPS) for g in groups.values()), 'Complete original step family required')
    results = []
    incompatible_stencils = []
    summaries = defaultdict(list)
    with localcontext() as context:
        context.prec = 80
        for key, group in sorted(groups.items()):
            require(len({r['candidate_hessian'] for r in group.values()}) == 1,
                    'Candidate curvature changed across steps')
            for coarse_h, fine_h in zip(STEPS[:-1], STEPS[1:]):
                coarse, fine = group[coarse_h], group[fine_h]
                if coarse['method'] != fine['method']:
                    incompatible_stencils.append(dict(case=key[0], column=key[1], score_axis=key[2],
                        complex_step=key[3], bound=key[4], coarse_real_step=coarse_h, fine_real_step=fine_h,
                        coarse_method=coarse['method'], fine_method=fine['method']))
                    continue
                ratio = (Decimal.from_float(coarse_h)/Decimal.from_float(fine_h))**4
                estimate = (ratio*Decimal.from_float(fine['independent_hessian'])
                            - Decimal.from_float(coarse['independent_hessian']))/(ratio-1)
                difference = estimate-Decimal.from_float(fine['candidate_hessian'])
                # These propagated estimates were never certified bounds.
                uncertainty = (ratio*Decimal.from_float(fine['abs_stencil_weighted_uncertified_error'])
                               + Decimal.from_float(coarse['abs_stencil_weighted_uncertified_error']))/(ratio-1)
                result = dict(case=key[0], column=key[1], score_axis=key[2], complex_step=key[3],
                    bound=key[4], method=fine['method'], coarse_real_step=coarse_h, fine_real_step=fine_h,
                    fourth_order_ratio=str(ratio), extrapolated_reference=str(estimate),
                    candidate_hessian=fine['candidate_hessian'], reference_minus_candidate=str(difference),
                    propagated_uncertified_reference_estimate=str(uncertainty),
                    within_nominal_1e_minus5_descriptive_only=abs(difference) <= Decimal('1e-5'),
                    higher_order_truncation_and_arithmetic_unbounded=True)
                results.append(result)
                summaries[(key[0], coarse_h, fine_h)].append(result)
    summary = [dict(case=k[0], coarse_real_step=k[1], fine_real_step=k[2], comparisons=len(values),
                    outside_nominal_budget=sum(not v['within_nominal_1e_minus5_descriptive_only'] for v in values),
                    max_absolute_difference=str(max(abs(Decimal(v['reference_minus_candidate'])) for v in values)))
               for k, values in summaries.items()]
    require(len(results)+len(incompatible_stencils) == 2352,
            'Every adjacent-step pair must be evaluated or explicitly retained as incompatible')
    packet = dict(input_sha256=args.input_sha, original_groups=588, original_records=2940,
        extrapolated_records=len(results), incompatible_stencils=incompatible_stencils,
        summaries=summary, records=results,
        postdiagnostic_analysis=True, best_step_selected=False, accuracy_accepted=False,
        certified_error_bounds=False, no_fit=True, no_calibration=True)
    with args.output.open('x') as handle:
        json.dump(packet, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(output=str(args.output), summaries=summary, accuracy_accepted=False)))


if __name__ == '__main__':
    main()
