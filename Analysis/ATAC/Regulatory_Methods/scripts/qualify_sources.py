#!/usr/bin/env python3
"""Qualify public sources, or join raw same-nucleus RNA by authoritative well/barcode."""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request


DEFAULT_ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
SOURCE_URLS = {
    'FNIH_readme': 'https://raw.githubusercontent.com/Gaulton-Lab/FNIH.Liver/main/README.md',
    'FNIH_portal': 'https://epigenome.wustl.edu/MASLD/',
    'LINGER_readme': 'https://raw.githubusercontent.com/Durenlab/LINGER/main/README.md',
    'TRIPOD_readme': 'https://raw.githubusercontent.com/yuchaojiang/TRIPOD/main/README.md',
}


def require(value, message):
    if not value:
        raise ValueError(message)


def read_table(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def fetch_small(url):
    """Read documentation only. Never follow a matrix/archive download."""
    request = urllib.request.Request(url, headers={'User-Agent': 'MASLD-source-qualification/1.0'})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = response.read(2 * 1024 * 1024 + 1)
            require(len(payload) <= 2 * 1024 * 1024, 'Documentation exceeds 2 MiB limit')
            return dict(url=url, final_url=response.url, status=response.status,
                        sha256=hashlib.sha256(payload).hexdigest(), bytes=len(payload),
                        text=payload.decode('utf-8'), observed=True)
    except (urllib.error.URLError, ValueError, UnicodeDecodeError) as error:
        return dict(url=url, observed=False, error=str(error))


def qualification_from_documents(documents):
    fn = documents.get('FNIH_readme', {})
    text = fn.get('text', '')
    portal = '\n'.join(record.get('text', '') for name, record in documents.items()
                       if name.startswith('FNIH_portal'))
    links = sorted(set(re.findall(r'''https?://[^\s<>"']+''', portal)))
    candidates = [link for link in links if re.search(r'\.h5(ad)?|\.rds|matrix|metadata|\.mtx|\.csv|\.tsv', link, re.I)]
    tracks = [link for link in links if re.search(r'\.(bw|bigwig|hic|bigbed)(?:$|[?])', link, re.I)]
    return dict(
        FNIH=dict(
            resource_role='potential human MASLD paired-multiome expansion',
            source_reported_86_livers=bool(re.search(r'86\s+livers', text)),
            processed_portal_link_verified='https://epigenome.wustl.edu/MASLD/' in text,
            source_says_raw_supplementary_forthcoming=bool(re.search(r'raw and supplementary data will be available soon', text, re.I)),
            portal_candidate_data_links=candidates,
            portal_aggregate_track_count=len(tracks),
            portal_aggregate_track_examples=tracks[:5],
            donor_demultiplexed_counts='unverified', phenotype_join='unverified',
            independence_from_GSE244832='unverified',
            qualified_for_MASH_replication=False,
            reason='Documentation/portal availability does not establish usable cell matrices, donor phenotype joins or donor independence.',
            source_owned_results='Published GRNs/QTLs remain source-owned and cannot serve as both priors and independent evaluation.'),
        other_expansion_sources=[
            dict(source='SenNet_liver', url='https://pmc.ncbi.nlm.nih.gov/articles/PMC12903365/',
                 role='potential fibrosis-context evaluation', state='not_qualified',
                 reason='Exact ATAC donor/phenotype roster and same-nucleus pairing require verification; source normal labels do not establish absence of fibrosis.'),
            dict(source='Hammond_E-MTAB-13131', url='https://link.springer.com/article/10.1038/s44321-025-00230-6',
                 role='cirrhotic regulatory context', state='not_qualified_for_MASH_replication',
                 reason='Three disease-only resection donors; no verified within-cohort healthy comparator.')],
        method_eligibility=[
            dict(method='LINGER', stage='conditional paired-multiome comparison', source=SOURCE_URLS['LINGER_readme'],
                 prerequisites=['exact same-nucleus RNA/ATAC', 'donor-grouped evaluation', 'external-prior identity and reference-overlap inventory'],
                 install_status='not_installed_by_this_analysis', installation_requires_explicit_approval=True,
                 native_method='yes only after version and inputs are qualified'),
            dict(method='TRIPOD', stage='conditional paired-multiome comparison', source=SOURCE_URLS['TRIPOD_readme'],
                 prerequisites=['exact same-nucleus RNA/ATAC', 'TF/peak/target identifiers', 'donor-aware evaluation of trio evidence'],
                 install_status='not_installed_by_this_analysis', installation_requires_explicit_approval=True)],
        source_exposure=[
            dict(cohort='GSE296875', status='previously examined in this repository',
                 permitted_role='internal held-donor paired-data evaluation', untouched_external=False,
                 MASH_contrast_qualified=False),
            dict(cohort='GSE244832', status='existing disease cohort; unpaired nuclei',
                 permitted_role='donor-level accessibility and separate-library RNA corroboration', same_nucleus=False),
            dict(cohort='GSE281367', status='existing ATAC-only cohort',
                 permitted_role='independent-donor corroboration with shared condition-blind feature ascertainment', same_nucleus=False)])


def qualify(args):
    output = Path(args.output)
    require(not output.exists(), 'Refusing existing source-qualification output')
    documents = {name: fetch_small(url) for name, url in SOURCE_URLS.items()}
    # The portal is a JS application; inspect its bounded same-origin entry script too.
    portal = documents['FNIH_portal'].get('text', '')
    scripts = re.findall(r'<script[^>]+src=["\']([^"\']+)["\']', portal)
    for index, relative in enumerate(scripts[:4]):
        url = urllib.parse.urljoin(SOURCE_URLS['FNIH_portal'], relative)
        if urllib.parse.urlparse(url).netloc == 'epigenome.wustl.edu':
            documents[f'FNIH_portal_script_{index}'] = fetch_small(url)
    output.mkdir(parents=True, exist_ok=False)
    for name, record in documents.items():
        if 'text' in record:
            (output / (name + '.txt')).write_text(record['text'])
    summary = qualification_from_documents(documents)
    summary['checked_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    summary['documents'] = {name: {k: v for k, v in record.items() if k != 'text'} for name, record in documents.items()}
    summary['large_data_downloaded'] = False
    write_json(output / 'source_qualification.json', summary)
    print(json.dumps(summary, indent=2))


def export_raw_rna(config, paired_root):
    import h5py
    import numpy as np
    from scipy import sparse
    from scipy.io import mmwrite

    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and os.environ.get('SLURM_JOB_NODELIST'),
            'Raw paired RNA extraction requires a SLURM compute allocation')
    root = Path(config.get('project_root', DEFAULT_ROOT))
    paired = config.get('paired', {})
    lock_path = Path(paired.get('source_lock', root / 'Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json'))
    lock = json.loads(lock_path.read_text())
    source_reference = lock_path.parent
    directories = [Path(paired_root)] if (Path(paired_root) / 'cells.tsv').is_file() else sorted(
        path for path in Path(paired_root).iterdir() if path.is_dir() and (path / 'cells.tsv').is_file())
    require(directories, 'No paired ATAC cell tables found')
    for directory in directories:
        require(not (directory / 'RNA_counts.mtx.gz').exists() and not (directory / 'RNA_export.json').exists(),
                'Refusing paired RNA overwrite')
    tables = {directory: read_table(directory / 'cells.tsv') for directory in directories}
    genes = {directory: read_table(directory / 'genes.tsv') for directory in directories}
    accumulators = {directory: [None] * len(cells) for directory, cells in tables.items()}
    source_evidence = []
    for well in sorted({row['well_id'] for cells in tables.values() for row in cells}):
        spec = lock['raw_rna_sources'][well]
        source = (source_reference / spec['path']).resolve()
        source_stat = source.stat()
        source_before = (source_stat.st_size, source_stat.st_mtime_ns, source_stat.st_ctime_ns)
        require(source.is_file() and source.stat().st_size == spec['bytes'] and sha(source) == spec['sha256'],
                'Raw RNA source differs: ' + well)
        with h5py.File(source, 'r') as handle:
            matrix = handle['matrix']
            decode = lambda values: [v.decode() if isinstance(v, bytes) else str(v) for v in values]
            types = decode(matrix['features/feature_type'][:])
            ids = decode(matrix['features/id'][:]); symbols = decode(matrix['features/name'][:])
            gene_positions = np.flatnonzero(np.array(types) == 'Gene Expression')
            axis = [(ids[i], symbols[i]) for i in gene_positions]
            require(len(axis) == 36601 and len(set(identity for identity, _ in axis)) == 36601, 'Raw RNA axis differs')
            axis_digest = hashlib.sha256()
            for identity, symbol in axis:
                axis_digest.update(identity.encode() + b'\0' + symbol.encode() + b'\n')
            require(axis_digest.hexdigest() == lock['raw_rna_registry']['rna_feature_id_and_symbol_sha256'], 'Raw RNA axis hash differs')
            barcodes = decode(matrix['barcodes'][:]); positions = {b: i for i, b in enumerate(barcodes)}
            require(len(positions) == len(barcodes), 'Duplicate raw barcode within well')
            feature_to_gene = np.full(len(types), -1, dtype=np.int64)
            feature_to_gene[gene_positions] = np.arange(len(gene_positions))
            pointers = matrix['indptr'][:]
            require(len(pointers) == len(barcodes) + 1 and pointers[0] == 0 and np.all(np.diff(pointers) >= 0), 'Invalid raw CSC pointers')
            selected_total = 0; selected_cells = 0
            for directory, cells in tables.items():
                require(axis == [(row['gene_id'], row['gene_symbol']) for row in genes[directory]], 'Export/raw RNA gene identities differ')
                for row_index, row in enumerate(cells):
                    if row['well_id'] != well:
                        continue
                    require(row['cell_id'] == well + '_' + row['raw_barcode'], 'Cell identity is not the exact authoritative well/barcode')
                    require(row['raw_barcode'] in positions, 'Paired nucleus absent from raw RNA: ' + row['cell_id'])
                    column = positions[row['raw_barcode']]
                    begin, end = int(pointers[column]), int(pointers[column + 1])
                    values = matrix['data'][begin:end]
                    source_indices = matrix['indices'][begin:end]
                    require(np.issubdtype(values.dtype, np.integer) and np.all(values >= 0), 'Raw RNA entries are not nonnegative integers')
                    gene_index = feature_to_gene[source_indices]
                    keep = gene_index >= 0
                    selected = sparse.csr_matrix((values[keep].astype(np.int64), gene_index[keep], np.array([0, int(keep.sum())])), shape=(1, len(axis)))
                    require(selected.has_canonical_format, 'Raw RNA sparse nucleus is not canonical')
                    require(accumulators[directory][row_index] is None, 'RNA nucleus assigned twice')
                    accumulators[directory][row_index] = selected
                    selected_total += int(values[keep].sum(dtype=np.int64)); selected_cells += 1
            source_evidence.append(dict(well=well, path=str(source), sha256=spec['sha256'],
                                        selected_cells_across_lineages=selected_cells,
                                        selected_UMIs_across_lineages=selected_total))
        source_stat = source.stat()
        require((source_stat.st_size, source_stat.st_mtime_ns, source_stat.st_ctime_ns) == source_before,
                'Raw RNA source changed during extraction: ' + well)
    for directory, chunks in accumulators.items():
        require(all(value is not None for value in chunks), 'Incomplete paired RNA coverage')
        counts = sparse.vstack(chunks, format='csr')
        require(counts.shape == (len(tables[directory]), len(genes[directory])), 'Paired RNA shape differs')
        temporary = directory / 'RNA_counts.mtx.gz.partial'
        require(not temporary.exists(), 'Refusing partial output overwrite')
        with gzip.open(temporary, 'wb') as handle:
            mmwrite(handle, counts, field='integer')
        temporary.rename(directory / 'RNA_counts.mtx.gz')
        # Sparse NPZ gives the shared evaluator a bounded native Python reader.
        require(not (directory / 'rna_counts.npz').exists(), 'Refusing RNA NPZ overwrite')
        sparse.save_npz(directory / 'rna_counts.npz', counts)
        write_json(directory / 'RNA_export.json', dict(cohort='GSE296875', cells=counts.shape[0], genes=counts.shape[1],
            same_nucleus_pairing='exact well and raw barcode from authoritative cell_identity.tsv',
            row_order='identical to cells.tsv and ATAC_counts.mtx.gz', RNA_count_unit='raw Gene Expression UMI counts',
            RNA_total_UMIs=int(counts.sum()), source_lock=str(lock_path), sources=source_evidence,
            MASH_contrast_qualified=False, untouched_external=False, orientation='cells_by_genes'))
        print(json.dumps(dict(directory=str(directory), RNA_cells=counts.shape[0], RNA_genes=counts.shape[1], RNA_total_UMIs=int(counts.sum()))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    q = sub.add_parser('qualify'); q.add_argument('--output', required=True)
    r = sub.add_parser('export-rna'); r.add_argument('--config', required=True); r.add_argument('--paired-root', required=True)
    args = parser.parse_args()
    if args.command == 'qualify':
        qualify(args)
    else:
        export_raw_rna(json.loads(Path(args.config).read_text()), args.paired_root)


if __name__ == '__main__':
    main()
