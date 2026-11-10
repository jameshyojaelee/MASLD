#!/usr/bin/env python3
"""Rebuild positive descriptive links from peak identifiers on the fixed axis."""
import argparse
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

A = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--lineage', required=True)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    source = args.candidate/f'disease_links/{args.cohort}_{args.lineage}'
    output = args.candidate/f'disease_links_index_corrected/{args.cohort}_{args.lineage}'
    output.mkdir(parents=True, exist_ok=False)
    frozen = A/'results/20261024/ridge_grid_sensitivity_r1/hepatocyte/prepared/source_snapshot/network_compare.py'
    spec = importlib.util.spec_from_file_location('annotation_core', frozen)
    core = importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
    peaks = core.read_table(A/f'results/20261023/inputs/count/{args.cohort}/{args.lineage}/peaks.tsv')
    peaks[['start0', 'end']] = peaks[['start0', 'end']].astype(np.int64)
    native = pd.read_csv(source/'cicero/peaks.tsv', sep='\t')
    columns = pd.Index(peaks.peak_id).get_indexer(native.peak_id)
    assert (columns >= 0).all() and len(set(columns)) == len(columns)
    with threadpool_limits(limits=1):
        reduced = mmread(source/'cicero/counts.mtx').T.tocoo()
    aggregate = sparse.csr_matrix((reduced.data, (reduced.row, columns[reduced.col])),
                                  shape=(reduced.shape[0], len(peaks)))
    assert int(aggregate.sum()) == int(reduced.sum())
    del reduced
    (output/'cicero').mkdir()
    for path in (source/'cicero').iterdir():
        if path.name != 'links.tsv.sparse_graph.json':
            (output/'cicero'/path.name).symlink_to(path.resolve())
    for path in source.iterdir():
        if path.name not in ['cicero', 'links.tsv.gz', 'coverage.tsv', 'completion.json']:
            (output/path.name).symlink_to(path.resolve())
    graph = core.CiceroGraph.from_tsv(output/'cicero/links.tsv', peaks)
    x = core.normalize_sparse(aggregate).tocsc()
    donors = pd.read_csv(source/'cicero/donors.tsv', sep='\t', dtype=str).donor_id.to_numpy()
    genes = core.read_table(A/'results/20261023/preparation/genes.tsv')
    genes.tss0 = genes.tss0.astype(np.int64)
    geometry = core.PeakGeometry(peaks)
    fields = ['cohort','lineage','gene_id','gene_symbol','peak_id','peak_coordinate','method','rank','score',
              'score_unit','distance_to_tss_bp','n_donors','evaluation_status']
    coverage = []; n_links = 0
    with gzip.open(output/'links.tsv.gz', 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for number, gene in enumerate(genes.sort_values('gene_id').to_dict('records'), 1):
            prom, distal = geometry.query(gene, 250000, 2000)
            counts = dict(distance=0, cicero=0)
            if len(prom) and len(distal):
                promoter = core.log_columns(x, prom).sum(axis=1)
                for method in ['distance', 'cicero']:
                    selected, scores = core.selected_indices(method, distal, prom, peaks, gene, x, promoter,
                                                             donors, 42, 10, graph if method == 'cicero' else None)
                    if method == 'cicero':
                        keep = scores > 0; selected, scores = selected[keep], scores[keep]
                    counts[method] = len(selected)
                    for rank, (index, score) in enumerate(zip(selected, scores), 1):
                        peak = peaks.iloc[index]
                        writer.writerow(dict(cohort=args.cohort, lineage=args.lineage, gene_id=gene['gene_id'],
                            gene_symbol=gene['gene_name'], peak_id=peak.peak_id,
                            peak_coordinate=f'{peak.chrom}:{peak.start0}-{peak.end}', method=method, rank=rank,
                            score=float(score), score_unit='negative_distance_bp' if method == 'distance' else 'donor_weighted_coaccessibility',
                            distance_to_tss_bp=abs((int(peak.start0)+int(peak.end))/2-int(gene['tss0'])),
                            n_donors=len(set(donors)), evaluation_status='all_donor_descriptive_accessibility_link_not_edge_validated'))
                        n_links += 1
            coverage.append(dict(gene_id=gene['gene_id'], gene_symbol=gene['gene_name'], n_promoter_peaks=len(prom),
                                 n_geometric_distal_peaks=len(distal), distance_links=counts['distance'], cicero_links=counts['cicero']))
            if number % 4000 == 0:
                print(args.cohort, args.lineage, number, n_links, flush=True)
    pd.DataFrame(coverage).to_csv(output/'coverage.tsv', sep='\t', index=False)
    # Distance selection is unaffected by the Cicero index correction.
    keys = ['gene_id', 'peak_id', 'rank', 'score']
    old = pd.read_csv(source/'links.tsv.gz', sep='\t')
    new = pd.read_csv(output/'links.tsv.gz', sep='\t')
    pd.testing.assert_frame_equal(old.loc[old.method.eq('distance'), keys].reset_index(drop=True),
                                  new.loc[new.method.eq('distance'), keys].reset_index(drop=True))
    assert new.loc[new.method.eq('cicero'), 'score'].gt(0).all()
    assert new.groupby(['gene_id','method']).size().le(10).all()
    # Independently verify selected scores by peak names in the raw edge table.
    linked_genes = sorted(new.loc[new.method.eq('cicero'), 'gene_id'].unique())
    selected_genes = set(np.random.default_rng(42).choice(linked_genes, min(200, len(linked_genes)), replace=False))
    check_rows = new[new.method.eq('cicero') & new.gene_id.isin(selected_genes)]
    expected = dict(zip(zip(check_rows.gene_id, check_rows.peak_id), check_rows.score))
    observed = {key:-np.inf for key in expected}
    promoter_genes = {}
    for gene in genes[genes.gene_id.isin(selected_genes)].to_dict('records'):
        prom, _ = geometry.query(gene, 250000, 2000)
        for peak in peaks.iloc[prom].peak_id:
            promoter_genes.setdefault(peak, set()).add(gene['gene_id'])
    for chunk in pd.read_csv(source/'cicero/links.tsv', sep='\t', chunksize=250000):
        chunk = chunk[chunk.peak1.isin(promoter_genes) | chunk.peak2.isin(promoter_genes)]
        for edge in chunk.itertuples(index=False):
            for promoter_peak, target in [(edge.peak1, edge.peak2), (edge.peak2, edge.peak1)]:
                for gene_id in promoter_genes.get(promoter_peak, ()):
                    key = (gene_id, target)
                    if key in observed:
                        observed[key] = max(observed[key], edge.coaccess)
    assert expected and all(np.isclose(observed[key], value, rtol=1e-12, atol=1e-12) for key, value in expected.items())
    (output/'independent_peak_identifier_check.json').write_text(json.dumps(dict(status='passed',
        seed=42, n_genes=len(selected_genes), n_selected_edges=len(expected),
        check='scores_reconstructed_from_raw_edge_peak_names_without_sparse_graph_indices'), indent=2))
    record = json.loads((source/'completion.json').read_text())
    record.update(n_links=n_links, source=str(source.resolve()), positive_Cicero_only=True,
                  graph_index='fixed_original_peak_axis_via_peak_identifier', distance_selection_reproduced=True,
                  prediction_feature_selection_changed=False)
    (output/'completion.json').write_text(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
