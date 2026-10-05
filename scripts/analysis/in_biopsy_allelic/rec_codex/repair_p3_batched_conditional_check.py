#!/usr/bin/env python3
"""Reuse immutable synthetic case reconstruction, isolated batched rows only."""
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
BATCHED = Path(__file__).with_name('repair_p3_batched_conditional.py')
V2 = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')
V2_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'
V2_CHECKER = Path(__file__).with_name('repair_p3_prepared_conditional_v2_check.py')
V2_CHECKER_SHA = '7443b61cea8565b6ea27320528d6a5ddcd088aff70cc0e8354e5055368013af1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=('failed0gene2','positive18'),required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-batched-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    guards = {SOURCE:SOURCE_SHA,Path(__file__).resolve():args.expected_checker_sha,
              BATCHED:args.expected_batched_sha,V2:V2_SHA,V2_CHECKER:V2_CHECKER_SHA}
    for path,digest in guards.items():
        if sha(path)!=digest:
            raise ValueError('Reviewed reconstruction/batch/dependency source drift: '+str(path))
        ast.parse(path.read_text())
    source = SOURCE.read_text()
    labels = {'repair_p3_fixed_point_':'repair_p3_batched_conditional_',
              "save('predictor_cache_comparison', report)":"save('batched_conditional_comparison', report)",
              'predictor_comparison_completed=True':'batched_conditional_comparison_completed=True'}
    for old,new in labels.items():
        if source.count(old)!=1:
            raise ValueError('Expected private label occurrence differs: '+old)
        source = source.replace(old,new)
    namespace = {'__file__':str(SOURCE),'__name__':'private_batched_conditional_reconstruction'}
    exec(compile(ast.parse(source),str(SOURCE),'exec'),namespace)
    namespace['PREDICTOR_CACHE'] = BATCHED
    # Existing guarded archive saves both new files and both immutable v2 sources.
    namespace['GUARDS'].update(guards)
    sys.argv = [str(SOURCE),'--case',args.case,'--predictor-points-sha',args.expected_batched_sha,
                '--expected-helper-sha',SOURCE_SHA,'--expected-kernel-sha',KERNEL_SHA]
    if args.inspect:
        sys.argv.append('--inspect')
    namespace['main']()


if __name__=='__main__':
    main()
