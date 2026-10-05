#!/usr/bin/env python3
"""Reuse guarded synthetic reconstruction; check conditional preparation only.

The archived fixed-point checker is executed in a private namespace. Only labels,
the preparation-hook source and output directory prefix differ. Its numerical
case reconstruction, source/environment guards and postchecks are unchanged.
No integration, fit or calibration is invoked by this route.
"""
import argparse
import ast
import hashlib
from pathlib import Path
import sys

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SOURCE = REC/'repair_p3_fixed_point_failed0gene2_22039132/executed_sources/repair_p3_fixed_point_check.py'
SOURCE_SHA = 'd3f529ad20666c58dc8a84d01785d3f62c42f18a430a3a572a055c180b9b4a9e'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
PREPARED = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('failed0gene2','positive18'), required=True)
    parser.add_argument('--expected-checker-sha', required=True)
    parser.add_argument('--expected-prepared-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    if sha(SOURCE) != SOURCE_SHA or sha(Path(__file__)) != args.expected_checker_sha or sha(PREPARED) != args.expected_prepared_sha:
        raise ValueError('Reviewed reconstruction/check/preparation source drift')
    source = SOURCE.read_text()
    # Nonmathematical label changes are exact and individually occurrence-guarded.
    labels = {'repair_p3_fixed_point_':'repair_p3_prepared_conditional_v2_',
              "save('predictor_cache_comparison', report)":"save('prepared_conditional_comparison', report)",
              'predictor_comparison_completed=True':'prepared_conditional_comparison_completed=True'}
    for old,new in labels.items():
        if source.count(old) != 1:
            raise ValueError('Expected private label occurrence differs: '+old)
        source = source.replace(old,new)
    namespace = {'__file__':str(SOURCE),'__name__':'private_prepared_reconstruction'}
    exec(compile(ast.parse(source),str(SOURCE),'exec'),namespace)
    namespace['PREDICTOR_CACHE'] = PREPARED
    namespace['GUARDS'][Path(__file__).resolve()] = args.expected_checker_sha
    # Existing checker archives and hash-checks both new files pre/post compute.
    sys.argv = [str(SOURCE),'--case',args.case,'--predictor-points-sha',args.expected_prepared_sha,
                '--expected-helper-sha',SOURCE_SHA,'--expected-kernel-sha',KERNEL_SHA]
    if args.inspect:
        sys.argv.append('--inspect')
    namespace['main']()


if __name__ == '__main__':
    main()
