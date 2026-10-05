#!/usr/bin/env python3
"""Guarded synthetic gene-local moment comparisons; no integrals or fits."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SOURCE = REC/'repair_p3_fixed_point_failed0gene2_22039132/executed_sources/repair_p3_fixed_point_check.py'
SOURCE_SHA = 'd3f529ad20666c58dc8a84d01785d3f62c42f18a430a3a572a055c180b9b4a9e'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
MOMENTS = Path(__file__).with_name('repair_p3_local_moments.py')
BATCHED = Path(__file__).with_name('repair_p3_batched_conditional.py')
BATCHED_SHA = '184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76'
V2 = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')
V2_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'
BASES = {
    'failed0gene2':('repair_p3_fixed_point_failed0gene2_22039200',
                  'edec063451383f9a1386619186c2089d3ef8e6ee34d854e80a2894fdd28e0f9a'),
    'positive18':('repair_p3_fixed_point_positive18_22039201',
                 '4a576370be2f33ac28bf23052fe283b4a5f614bc30e983bd6a84224da857c722')}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=tuple(BASES),required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-moments-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    base_dir,base_sha = BASES[args.case]
    baseline_path = REC/base_dir/'candidate_baseline.json'
    guards = {SOURCE:SOURCE_SHA,Path(__file__).resolve():args.expected_checker_sha,
              MOMENTS:args.expected_moments_sha,BATCHED:BATCHED_SHA,V2:V2_SHA,
              baseline_path:base_sha}
    for path,digest in guards.items():
        if sha(path)!=digest:
            raise ValueError('Reviewed local moment/reconstruction dependency changed: '+str(path))
    source = SOURCE.read_text()
    labels = {'repair_p3_fixed_point_':'repair_p3_local_moments_',
              "save('predictor_cache_comparison', report)":"save('local_moment_comparison', report)",
              'predictor_comparison_completed=True':'local_moment_comparison_completed=True'}
    for old,new in labels.items():
        if source.count(old)!=1:
            raise ValueError('Immutable reconstruction output label differs')
        source=source.replace(old,new)
    namespace = {'__file__':str(SOURCE),'__name__':'private_local_moment_reconstruction'}
    exec(compile(ast.parse(source),str(SOURCE),'exec'),namespace)
    namespace['PREDICTOR_CACHE'] = MOMENTS
    namespace['GUARDS'].update(guards)
    original_loader = namespace['load_module']
    def local_loader(name,path):
        module = original_loader(name,path)
        if Path(path).name == MOMENTS.name:
            # Called after the immutable synthetic reconstruction saves its protocol.
            output = REC/('repair_p3_local_moments_'+args.case+'_'+os.environ['SLURM_JOB_ID'])
            protocol = json.loads((output/'protocol.json').read_text())
            if protocol['fixed_v'] <= 0 or protocol['fixed_gene_rows']!=24:
                raise ValueError('Expected positive fixed synthetic point')
            module.SAVED_VARIANCE = protocol['fixed_v']
            module.SAVED_LOGMASS_ESTIMATE = json.loads(baseline_path.read_text())['logmass']
            original_loader('repair_p3_local_validated_batch',Path(path).with_name(BATCHED.name))
        return module
    namespace['load_module'] = local_loader
    sys.argv = [str(SOURCE),'--case',args.case,'--predictor-points-sha',args.expected_moments_sha,
                '--expected-helper-sha',SOURCE_SHA,'--expected-kernel-sha',KERNEL_SHA]
    if args.inspect:sys.argv.append('--inspect')
    namespace['main']()


if __name__=='__main__':
    main()
