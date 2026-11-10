"""Compare all numerical network outputs before/after training-row caching.

Run on compute with explicit reference and candidate source files. The small
Cicero graph is fixed in both runs; this tests row-cache equivalence, not Cicero.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy import sparse


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Run on an allocated compute node')
    args.output.mkdir(exist_ok=False, parents=True)
    rng = np.random.default_rng(42)
    counts = sparse.csr_matrix(rng.poisson(5, size=(30, 12)))
    meta = pd.DataFrame(dict(donor_id=np.repeat([f'd{i}' for i in range(10)], 3),
                             cohort='GSE296875', lineage='hepatocyte'))
    peaks = pd.DataFrame(dict(peak_id=[f'p{i}' for i in range(12)], chrom='chr1',
                              start0=np.arange(12) * 2000, end=np.arange(12) * 2000 + 500))
    genes = pd.DataFrame([dict(gene_id=f'g{i}', chrom='chr1', tss0=2000*i,
                              strand='+' if i % 2 else '-') for i in range(1, 11)])
    genes = pd.concat([genes, pd.DataFrame([
        dict(gene_id='missing_promoter', chrom='chr1', tss0=900000, strand='+'),
        dict(gene_id='missing_chromosome', chrom='chrX', tss0=10000, strand='-')])], ignore_index=True)
    rna = sparse.csr_matrix(rng.poisson(5, size=(30, 10)))
    rna_genes = pd.DataFrame(dict(gene_id=[f'g{i}' for i in range(1, 11)]))
    graph = pd.DataFrame([dict(peak1=f'p{i}', peak2=f'p{j}', coaccess=(i+j+1)/24)
                          for i in range(12) for j in range(i+1, 12)])
    settings = dict(cicero_rscript='/usr/bin/true', bootstrap_reps=100,
                    max_genes=0, outer_folds=5, inner_folds=3)
    runs = {}
    count_arrays = [getattr(counts, field).copy() for field in ['data', 'indices', 'indptr']]
    rna_arrays = [getattr(rna, field).copy() for field in ['data', 'indices', 'indptr']]
    original_metadata, original_peaks = meta.copy(deep=True), peaks.copy(deep=True)
    for name, path in [('reference', args.reference), ('candidate', args.candidate)]:
        source_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        module = load(path, name + '_cache_check')
        module.cicero_scores = lambda *unused, **kwargs: graph
        worker_pids = {}
        if name == 'candidate' and hasattr(module, '_fit_gene_readout'):
            original_gene = module._fit_gene_readout
            original_ordered = module._ordered_gene_readouts
            def record_pid(job, context):
                result = original_gene(job, context)
                result['_fixture_pid'] = os.getpid()
                return result
            def collect_pids(jobs, context, workers):
                for result in original_ordered(jobs, context, workers):
                    fold = result['coverage']['outer_fold']
                    worker_pids.setdefault(fold, set()).add(result.pop('_fixture_pid'))
                    yield result
            module._fit_gene_readout = record_pid
            module._ordered_gene_readouts = collect_pids
        start = time.monotonic()
        local_settings = dict(settings)
        if name == 'candidate':
            local_settings['readout_workers'] = args.workers
        module.fit_group(counts, meta, peaks, genes, args.output/name, local_settings,
                         42, rna, rna_genes)
        for field, expected in zip(['data', 'indices', 'indptr'], count_arrays):
            np.testing.assert_array_equal(getattr(counts, field), expected)
        for field, expected in zip(['data', 'indices', 'indptr'], rna_arrays):
            np.testing.assert_array_equal(getattr(rna, field), expected)
        pd.testing.assert_frame_equal(meta, original_metadata, check_exact=True)
        pd.testing.assert_frame_equal(peaks, original_peaks, check_exact=True)
        if hashlib.sha256(path.read_bytes()).hexdigest() != source_sha:
            raise RuntimeError('Source changed during numerical equivalence check')
        if name == 'candidate' and any(len(pids) != args.workers for pids in worker_pids.values()):
            raise RuntimeError('Fixture worker counts by fold: '+str({k: len(v) for k,v in worker_pids.items()}))
        runs[name] = dict(seconds=time.monotonic()-start,
                          source_sha256=source_sha,
                          observed_worker_pids_by_fold={str(k):sorted(v) for k,v in worker_pids.items()})
    checked = []
    for name in ['predictions.tsv', 'links.tsv', 'comparisons.tsv', 'donor_errors.tsv',
                 'coverage.tsv', 'fold_assignments.tsv', 'method_execution_summary.tsv']:
        first = pd.read_csv(args.output/'reference'/name, sep='\t')
        second = pd.read_csv(args.output/'candidate'/name, sep='\t')
        pd.testing.assert_frame_equal(first, second, check_exact=True)
        checked.append(dict(table=name, n_rows=len(first), exact=True))
    receipt = dict(status='passed', scope='fixed_fixture_cache_equivalence',
                    numerical_outputs_exact=True, both_endpoints=True,
                    methods=['distance', 'directnet_adapted', 'cicero'],
                    n_outer_folds=5, n_inner_folds=3, runs=runs, tables=checked,
                    candidate_readout_workers=args.workers, input_arrays_and_metadata_unchanged=True,
                    cicero_scope='fixed_graph_only; native_inference_not_tested_here')
    (args.output/'verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
