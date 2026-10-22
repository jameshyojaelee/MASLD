#!/usr/bin/env python3
"""Compute-only GL polynomial and closed Gaussian moment checks."""
import argparse
from decimal import Decimal,localcontext
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import time


def require(ok,message):
    if not ok:
        raise ValueError(message)


def arctangent(x, tolerance):
    total, power, square = x, x, x*x
    for k in range(1,10000):
        power *= -square
        term = power/(2*k+1)
        total += term
        if abs(term) < tolerance:
            return total
    raise ValueError('Machin pi series did not converge')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module',type=Path,required=True)
    parser.add_argument('--module-sha',required=True)
    parser.add_argument('--self-sha',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    sha = lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    require(sha(Path(__file__)) == args.self_sha and sha(args.module) == args.module_sha,
            'Reviewed rule/checker source differs')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    args.out.mkdir(exist_ok=False)
    spec = importlib.util.spec_from_file_location('checked_GL',args.module)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    began = time.perf_counter()
    records = []
    for digits in (80,110):
        with localcontext() as context:
            context.prec = digits
            tolerance = Decimal(10)**(-digits+5)
            pi = 16*arctangent(Decimal(1)/5,tolerance)-4*arctangent(Decimal(1)/239,tolerance)
            variance = Decimal('0.5')
            center,conditional_cross = Decimal('0.75'),Decimal(2)/7
            mean,posterior_variance = center*variance/(1+variance),variance/(1+variance)
            second = mean*mean+posterior_variance
            third = mean**3+3*mean*posterior_variance
            mass = (2*pi/(1+variance)).sqrt()*(-center*center/(2*(1+variance))).exp()
            expected = (mass,mass*mean,mass*second,mass*(second-variance)/(2*variance**2),
                        mass*(conditional_cross+third),mass*(third-variance*mean)/(2*variance**2))
            expected_mixed = conditional_cross+2*mean*posterior_variance
            expected_r_variance = mean*posterior_variance/variance**2
            require(all(x != 0 for x in expected),'Mock has an unexercised zero moment')
            def mock(prepared,u):
                return -(u-center)**2/2,u,u*u,conditional_cross
            for rule_index in range(3):
                moments = [Decimal(0) for _ in range(6)]
                panel_ids = []
                for chunk in range(16):
                    value = module.integrate(None,variance,Decimal(0),rule_index,chunk,digits,mock)
                    panel_ids.extend(range(value['first_panel'],value['last_panel_exclusive']))
                    moments = [a+Decimal(b) for a,b in zip(moments,value['moments'])]
                require(panel_ids == list(range(value['total_panels'])),'Panel coverage differs')
                differences = [a-b for a,b in zip(moments,expected)]
                require(max(map(abs,differences)) < Decimal('1e-40'),'Closed raw moments differ')
                mixed = moments[4]/moments[0]-moments[1]*moments[2]/moments[0]**2
                r_variance = moments[5]/moments[0]-moments[1]*moments[3]/moments[0]**2
                require(abs(mixed-expected_mixed) < Decimal('1e-40'),'Mixed covariance differs')
                require(abs(r_variance-expected_r_variance) < Decimal('1e-40'),'Variance covariance differs')
                records.append(dict(digits=digits,rule_index=rule_index,
                    raw_moment_differences=list(map(str,differences)),
                    r_omega_curvature=str(mixed),r_variance_curvature=str(r_variance)))
                print(json.dumps(dict(digits=digits,rule_index=rule_index,passed=True)),flush=True)
    require(sha(Path(__file__)) == args.self_sha and sha(args.module) == args.module_sha,
            'Reviewed rule/checker source changed during execution')
    result = dict(completed=True,closed_Gaussian_and_polynomial_checks_passed=True,records=records,
        source_sha256={str(Path(__file__)):args.self_sha,str(args.module):args.module_sha},
        elapsed_seconds=time.perf_counter()-began,
        measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        no_native_integral_accuracy=True,no_fit=True,no_calibration=True)
    with (args.out/'summary.json').open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    main()
