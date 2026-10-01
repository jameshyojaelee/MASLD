#!/usr/bin/env python3
"""Fetch one authorized deposited accessibility matrix and inspect it locally."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[3]
KEY = 'liver_ATAC_peakCounts_DESeqDatasets_WASP-filtered_consensusPeaks_autosomalOnly_refinedBoundaries.RData'
SIZE = 811096877
MD5 = '6f3af5bac0509909d015e8fd99cc1d5e'
META_URL = 'https://zenodo.org/api/records/15025748'


def dump(path, obj):
    with path.open('x') as f:
        json.dump(obj, f, indent=2)
        f.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--rscript', required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Acquisition/inspection must run in the authorized compute job')
    args.out.mkdir(parents=True, exist_ok=False, mode=0o700)
    private = args.out / 'private_source'
    private.mkdir(mode=0o700)
    inspector = Path(__file__).with_name('data_currin_profile_inspect.R')
    shutil.copy2(__file__, args.out / 'executed_data_currin_profile_acquire.py')
    shutil.copy2(inspector, args.out / 'executed_data_currin_profile_inspect.R')
    with urllib.request.urlopen(META_URL, timeout=90) as response:
        data = response.read(2000001)
    if len(data) > 2000000:
        raise ValueError('Metadata size exceeds allowed bounded request')
    (args.out / 'source_metadata.json').write_bytes(data)
    metadata = json.loads(data)
    entry, = [x for x in metadata['files'] if x['key'] == KEY]
    if entry['size'] != SIZE or entry['checksum'] != 'md5:' + MD5:
        raise ValueError('Deposited size/checksum changed; no acquisition performed')
    if metadata['metadata'].get('access_right') != 'open':
        raise ValueError('Source is not the inspected public open deposit')
    url = entry['links']['self']
    if not url.startswith('https://zenodo.org/api/records/15025748/files/'):
        raise ValueError('Source asset URL differs from inspected record')
    dump(args.out / 'source_terms_and_use.json', dict(
        source_url=META_URL, record_doi='10.5281/zenodo.15025748',
        metadata_bytes=len(data), metadata_sha256=hashlib.sha256(data).hexdigest(),
        access_right='open', license=metadata['metadata'].get('license'),
        terms_limit='No explicit reuse license is declared in this record metadata. Local authorized research inspection only; no donor-matrix redistribution or public rehosting.',
        use='Measured native accessibility profile assessment readiness; separate from caQTL beta and prediction drift.',
        exposure='Same Currin source family as development labels; no independent validation claim; exact foundation-model source exposure unresolved.',
        projected_peak_disk_bytes=2000000000, genotype_download=False, Wenz_outcomes=False))
    dest = private / KEY
    n, md5, sha = 0, hashlib.md5(), hashlib.sha256()
    receipt = dict(url=url, expected_bytes=SIZE, expected_md5=MD5, status='incomplete')
    try:
        with urllib.request.urlopen(url, timeout=90) as response, dest.open('xb') as out:
            for block in iter(lambda: response.read(4 << 20), b''):
                n += len(block)
                if n > SIZE:
                    raise ValueError('Download exceeds exact expected asset size')
                out.write(block)
                md5.update(block)
                sha.update(block)
        if n != SIZE or md5.hexdigest() != MD5:
            raise ValueError('Asset byte-size or MD5 mismatch')
        receipt['status'] = 'verified_complete'
    finally:
        receipt.update(bytes_received=n, md5=md5.hexdigest(), sha256=sha.hexdigest())
        dump(args.out / 'asset_download_receipt.json', receipt)
    peaks = ROOT / 'GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/supplementalData1_liver_ATAC_peaks.bed.gz'
    dump(args.out / 'inspection_input_receipt.json', dict(
        peak_metadata_path=str(peaks), peak_metadata_sha256=hashlib.sha256(peaks.read_bytes()).hexdigest(),
        executed_R_sha256=hashlib.sha256(inspector.read_bytes()).hexdigest(),
        executed_python_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    subprocess.run([args.rscript, str(args.out / 'executed_data_currin_profile_inspect.R'), str(dest), str(peaks), str(args.out)], check=True)
    total = sum(p.stat().st_size for p in args.out.rglob('*') if p.is_file())
    if total > 2000000000:
        raise ValueError('Output exceeds authorized 2GB disk cap; retained for inspection')
    dump(args.out / 'completion.json', dict(status='acquisition_and_aggregate_inspection_complete',
        output_bytes=total, job_id=os.environ['SLURM_JOB_ID'], profile_accuracy_not_run=True))


if __name__ == '__main__':
    main()
