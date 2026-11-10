#!/usr/bin/env python3
"""Inspect local RNA molecules and authoritative sample/run metadata on compute."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import urllib.request

import h5py
import numpy as np
import openpyxl
import pandas as pd

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
    args.output.mkdir(exist_ok=False)
    pairing = pd.read_csv(ROOT / 'data/GSE244832/metadata/donor_pairing.csv', dtype=str)
    source = ROOT / 'Analysis/SingleCell/reviews/control-label-review-20261021/sources'
    book = openpyxl.load_workbook(source / 'GSE244832_240605_18liver_metadata.xlsx', read_only=True, data_only=True)
    sheets = {sheet.title: list(sheet.values) for sheet in book}
    (args.output / 'source_workbook.json').write_text(json.dumps(sheets, indent=2, default=str))
    book.close()
    soft = (source / 'GSE244832_family.soft').read_text()
    samples = {}
    for block in soft.split('^SAMPLE = ')[1:]:
        lines = block.splitlines()
        samples[lines[0]] = [line for line in lines if line.startswith(('!Sample_title', '!Sample_relation', '!Sample_extract_protocol', '!Sample_data_processing', '!Sample_characteristics'))]
    (args.output / 'sample_evidence.json').write_text(json.dumps(samples, indent=2))
    rows, molecule_examples = [], {}
    for donor in pairing.to_dict('records'):
        for run in donor['rna_srrs'].split(';'):
            out = ROOT / 'data/GSE244832/cellranger' / run / 'outs'
            matrix, molecule = out / 'filtered_feature_bc_matrix.h5', out / 'molecule_info.h5'
            fastqs = sorted((ROOT / 'data/GSE244832/fastq' / run).glob('*.fastq.gz'))
            row = dict(cohort='GSE244832', donor_id=donor['donor_id'], specimen=donor['rna_sample'],
                       library_id=donor['rna_sample'], run=run, gsm=donor['rna_gsm'], condition=donor['condition'],
                       cross_assay_confidence=donor['confidence'], matrix=str(matrix), matrix_exists=matrix.exists(),
                       molecule=str(molecule), molecule_exists=molecule.exists(), fastq_bytes=sum(p.stat().st_size for p in fastqs))
            if matrix.exists():
                with h5py.File(matrix, 'r') as h:
                    m = h['matrix']
                    row.update(n_genes=int(m['shape'][0]), n_called_cells=int(m['shape'][1]),
                               total_umis=int(np.sum(m['data'][:], dtype=np.int64)), matrix_dtype=str(m['data'].dtype))
            if molecule.exists() and donor['rna_sample'] == 'JB_303':
                with h5py.File(molecule, 'r') as h:
                    molecule_examples[run] = {key: dict(shape=list(value.shape), dtype=str(value.dtype))
                                               if isinstance(value, h5py.Dataset) else list(value.keys()) for key, value in h.items()}
                    for key in ['library_info', 'metrics_json']:
                        if key in h:
                            value = h[key][()]
                            molecule_examples[run][key + '_value'] = value.decode() if isinstance(value, bytes) else str(value)
            rows.append(row)
    pd.DataFrame(rows).to_csv(args.output / 'rna_run_inventory.tsv', sep='\t', index=False)
    (args.output / 'molecule_examples.json').write_text(json.dumps(molecule_examples, indent=2))
    read_examples = []
    for path in sorted((ROOT / 'data/GSE244832/fastq/SRR26313200').glob('*.fastq.gz')):
        with gzip.open(path, 'rt') as handle:
            records = [[next(handle).strip() for _ in range(4)] for _ in range(20)]
        read_examples.append(dict(path=str(path), read_lengths=[len(r[1]) for r in records], first_header=records[0][0]))
    (args.output / 'fastq_structure.json').write_text(json.dumps(read_examples, indent=2))
    # Resolve accession exclusively from the deposited source.
    import re
    projects = re.findall(r'BioProject: (PRJNA\d+)', soft)
    if not projects:
        projects = re.findall(r'bioproject/(PRJNA\d+)', soft)
    assert len(set(projects)) == 1, 'Source must specify one BioProject'
    url = (f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={projects[0]}&result=read_run&format=tsv'
           '&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,library_name,library_strategy,library_layout,sample_title,fastq_bytes')
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            payload = response.read(4 * 1024 * 1024)
        (args.output / 'ena_run_report.tsv').write_bytes(payload)
        retrieval = dict(url=url, success=True)
    except Exception as error:
        retrieval = dict(url=url, success=False, error=str(error))
    (args.output / 'retrieval.json').write_text(json.dumps(retrieval, indent=2))
    summary = dict(status='completed', n_runs=len(rows), n_libraries=len(pairing),
                   matrices=sum(r['matrix_exists'] for r in rows), molecules=sum(r['molecule_exists'] for r in rows),
                   fastq_bytes=sum(r['fastq_bytes'] for r in rows), cellranger=shutil.which('cellranger'),
                   known_cellranger_exists=Path('/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger').exists(),
                   source_projects=projects, scientific_adoption=False)
    (args.output / 'completion.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
