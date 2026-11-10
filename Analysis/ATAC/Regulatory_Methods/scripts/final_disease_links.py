#!/usr/bin/env python3
"""Outcome-blind all-donor ATAC link maps; separate from prediction folds."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd

A = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', choices=['GSE244832','GSE281367'], required=True)
    p.add_argument('--lineage', choices=['hepatocyte','stellate'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=8)
    args = p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    frozen = A / 'results/20261024/ridge_grid_sensitivity_r1/hepatocyte/prepared/source_snapshot/network_compare.py'
    spec = importlib.util.spec_from_file_location('frozen_atac_core', frozen)
    core = importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=False)
    source = A / f'results/20261023/inputs/count/{args.cohort}/{args.lineage}'
    counts, cells, peaks = core.load_input(source,42,1000)
    config = json.loads((A / 'config/analysis.json').read_text())
    keep, exclusions = core.condition_eligibility(cells,config['network'])
    cells = cells.loc[keep].reset_index(drop=True); counts = counts[np.flatnonzero(keep)]
    counts = counts.copy(); counts.eliminate_zeros(); counts.data[:] = 1
    aggregate, meta, groups, eligibility = core.build_metacells(counts,cells,42,50,1000)
    meta.to_csv(args.output / 'aggregates.tsv',sep='\t',index=False)
    eligibility.to_csv(args.output / 'donor_eligibility.tsv',sep='\t',index=False)
    exclusions.to_csv(args.output / 'condition_exclusions.tsv',sep='\t',index=False)
    roster = []
    for i, group in enumerate(groups):
        for n in group:
            roster.append(dict(cohort=args.cohort,lineage=args.lineage,donor_id=cells.iloc[n].donor_id,
                               metacell_id=meta.iloc[i].metacell_id,cell_id=cells.iloc[n].cell_id))
    pd.DataFrame(roster).to_csv(args.output/'cell_roster.tsv.gz',sep='\t',index=False)
    genes = core.read_table(A / 'results/20261023/preparation/genes.tsv')
    genes.tss0 = genes.tss0.astype(np.int64)
    geometry = core.PeakGeometry(peaks)
    nonzero = np.flatnonzero(np.asarray(aggregate.sum(axis=0)).ravel() > 0)
    graph = core.cicero_scores(aggregate[:,nonzero],peaks.iloc[nonzero],args.output/'cicero',
        config['network']['cicero_rscript'],42,500000,distance_samples=100,donor_ids=meta.donor_id.to_numpy(),workers=args.workers)
    log = (args.output/'cicero/run.log').read_text()
    assert re.findall(r'Models with errors:\s*(\d+)',log) == ['0']
    x = core.normalize_sparse(aggregate).tocsc()
    # CiceroGraph.from_tsv preserves peaks.index, including omitted zero columns.
    donors = meta.donor_id.to_numpy()
    import gzip
    import csv
    columns = ['cohort','lineage','gene_id','gene_symbol','peak_id','peak_coordinate','method','rank','score','score_unit',
               'distance_to_tss_bp','n_donors','evaluation_status']
    coverage = []
    n_links = 0
    with gzip.open(args.output/'links.tsv.gz','wt') as handle:
        writer = csv.DictWriter(handle,fieldnames=columns,delimiter='\t'); writer.writeheader()
        for number,gene in enumerate(genes.sort_values('gene_id').to_dict('records'),1):
            prom, distal = geometry.query(gene,250000,2000)
            if not len(prom) or not len(distal):
                coverage.append(dict(gene_id=gene['gene_id'],gene_symbol=gene['gene_name'],
                    n_promoter_peaks=len(prom),n_geometric_distal_peaks=len(distal),distance_links=0,cicero_links=0))
                continue
            promoter = core.log_columns(x,prom).sum(axis=1)
            counts_by_method = {}
            for method in ['distance','cicero']:
                local_graph = graph if method == 'cicero' else None
                chosen, values = core.selected_indices(method,distal,prom,peaks,gene,x,promoter,donors,42,10,local_graph)
                if method == 'cicero':
                    positive = values > 0
                    chosen, values = chosen[positive], values[positive]
                counts_by_method[method] = len(chosen)
                for rank,(index,score) in enumerate(zip(chosen,values),1):
                    peak = peaks.iloc[index]
                    coordinate = f'{peak.chrom}:{peak.start0}-{peak.end}'
                    writer.writerow(dict(cohort=args.cohort,lineage=args.lineage,gene_id=gene['gene_id'],
                        gene_symbol=gene['gene_name'],peak_id=peak.peak_id,peak_coordinate=coordinate,method=method,
                        rank=rank,score=float(score),score_unit='negative_distance_bp' if method=='distance' else 'donor_weighted_coaccessibility',
                        distance_to_tss_bp=abs((int(peak.start0)+int(peak.end))/2-int(gene['tss0'])),
                        n_donors=len(set(donors)),evaluation_status='all_donor_descriptive_accessibility_link_not_edge_validated'))
                    n_links += 1
            coverage.append(dict(gene_id=gene['gene_id'],gene_symbol=gene['gene_name'],n_promoter_peaks=len(prom),
                n_geometric_distal_peaks=len(distal),distance_links=counts_by_method['distance'],cicero_links=counts_by_method['cicero']))
            if number % 2000 == 0:
                print(f'{args.cohort} {args.lineage} annotation gene={number} links={n_links}',flush=True)
    pd.DataFrame(coverage).to_csv(args.output/'coverage.tsv',sep='\t',index=False)
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed',cohort=args.cohort,lineage=args.lineage,
        n_donors=len(set(donors)),n_aggregates=len(meta),n_links=n_links,elapsed_seconds=time.monotonic()-started,
        hypothesis_testing=False,scientific_adoption=False),indent=2))


if __name__ == '__main__':
    main()
