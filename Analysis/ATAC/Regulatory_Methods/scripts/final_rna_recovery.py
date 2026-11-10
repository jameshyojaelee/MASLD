#!/usr/bin/env python3
"""Jointly recount sequencing runs from one source-declared RNA capture library."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np
import pandas as pd
from final_rna_identity import source_rna_libraries

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
CELLRANGER = Path('/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger')
REFERENCE = Path('/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--cpus', type=int, default=8)
    p.add_argument('--memory', type=int, default=56)
    p.add_argument('--attempt', type=int, default=1)
    args = p.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
    suffix = '' if args.attempt == 1 else f'_attempt{args.attempt}'
    output = args.candidate / 'rna_recovery' / (args.library + suffix)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    pairing = pd.read_csv(ROOT / 'data/GSE244832/metadata/donor_pairing.csv', dtype=str)
    donor = pairing.set_index('rna_sample').loc[args.library]
    source_identity = source_rna_libraries(args.candidate)[args.library]
    assert source_identity['gsm'] == donor.rna_gsm
    runs = donor.rna_srrs.split(';')
    ena = pd.read_csv(args.candidate / 'rna_source_audit/ena_run_report.tsv', sep='\t', dtype=str)
    subset = ena[ena.run_accession.isin(runs)]
    assert set(subset.run_accession) == set(runs)
    assert subset.experiment_accession.nunique() == subset.sample_accession.nunique() == 1
    assert set(subset.library_name) == {donor.rna_gsm}
    assert set(subset.sample_title.str.split(',').str[0]) == {args.library}
    fastqs = output / 'fastqs'; fastqs.mkdir()
    evidence, chemistry = [], set()
    for run in runs:
        old = ROOT / 'data/GSE244832/cellranger' / run / 'outs/molecule_info.h5'
        with h5py.File(old, 'r') as h:
            value = h['metrics_json'][()]
            metrics = json.loads(value.decode() if isinstance(value, bytes) else value)
            chemistry.add(metrics['chemistry_name'])
        directory = ROOT / 'data/GSE244832/fastq' / run
        for mate in [1,2]:
            choices = [directory / f'{run}_S1_L001_R{mate}_001.fastq.gz', directory / f'{run}_{mate}.fastq.gz']
            present = [q for q in choices if q.is_file()]
            assert present, f'Missing mate {mate}: {run}'
            source = present[0].resolve()
            with gzip.open(source, 'rt') as handle:
                record = [next(handle).strip() for _ in range(4)]
            assert record[0].startswith('@') and len(record[1]) >= (28 if mate == 1 else 25)
            target = fastqs / f'{run}_S1_L001_R{mate}_001.fastq.gz'
            target.symlink_to(source)
            evidence.append(dict(run=run, mate=mate, path=str(source), bytes=source.stat().st_size,
                                 first_read_length=len(record[1])))
    assert chemistry == {'SC3Pv3-polyA'}, f'Unexpected chemistry {chemistry}'
    # SC3Pv3 is the installed Cell Ranger option; the deposited 3-prime v3.1
    # protocol and all prior run-level molecule files support this read structure.
    help_text = subprocess.check_output([str(CELLRANGER), 'count', '--help'], text=True)
    (output / 'cellranger_count_help.txt').write_text(help_text)
    assert 'SC3Pv3' in help_text and '--include-introns' in help_text
    command = [str(CELLRANGER), 'count', '--id=' + args.library,
               '--transcriptome=' + str(REFERENCE), '--fastqs=' + str(fastqs.resolve()),
               '--sample=' + ','.join(runs),
               '--chemistry=SC3Pv3', '--include-introns=true', '--create-bam=false', '--nosecondary',
               '--localcores=' + str(args.cpus), '--localmem=' + str(args.memory),
               '--output-dir=' + str((output / 'count').resolve())]
    record = dict(library=args.library, donor_id=donor.donor_id, condition=source_identity['condition'],
                  local_crosswalk_condition=donor.condition,source_identity=source_identity,
                  cross_assay_confidence=donor.confidence, runs=runs, ena=subset.to_dict('records'),
                  fastqs=evidence, chemistry=list(chemistry), command=command,
                  reference=str(REFERENCE), reference_description=json.loads((REFERENCE / 'reference.json').read_text()),
                  cellranger_version=subprocess.check_output([str(CELLRANGER), '--version'], text=True).strip(),
                  no_run_matrix_summation=True, same_nucleus_pairing=False)
    (output / 'input_specification.json').write_text(json.dumps(record, indent=2))
    with (output / 'cellranger.log').open('w') as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True, cwd=output)
    matrix_path = output / 'count/outs/filtered_feature_bc_matrix.h5'
    with h5py.File(matrix_path, 'r') as h:
        matrix = h['matrix']; data = matrix['data'][:]
        barcodes = matrix['barcodes'].asstr()[:]
        assert len(barcodes) == len(set(barcodes))
        assert np.isfinite(data).all() and (data >= 0).all() and np.equal(data, np.floor(data)).all()
        summary = dict(status='completed', library=args.library, n_runs=len(runs), n_cells=len(barcodes),
                       n_genes=int(matrix['shape'][0]), total_umis=int(np.sum(data, dtype=np.int64)),
                       elapsed_seconds=time.monotonic()-started, same_nucleus_pairing=False, scientific_adoption=False)
    pd.read_csv(output / 'count/outs/metrics_summary.csv').to_csv(output / 'metrics.tsv', sep='\t', index=False)
    (output / 'completion.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
