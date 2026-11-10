"""Inspect public paired-liver inputs before any prediction fitting."""
import concurrent.futures
import gzip
import json
import os
from pathlib import Path
import urllib.request
import fsspec
import h5py

assert os.environ.get('SLURM_JOB_ID')
root = Path('Analysis/ATAC/Regulatory_Methods/results/20261024/external_liver_qualification_r1')
data = json.loads((root / 'sennet_public_liver.json').read_text())
datasets = [a for a in data if 'Multiome' in a['dataset_type']]

def inspect(a):
    record = dict(uuid=a['uuid'], sennet_id=a['sennet_id'],
                  donors=[s['sennet_id'] for s in a['sources']],
                  public_files=[f['rel_path'] for f in a['files'] if 'hubmap_ui/' not in f['rel_path']])
    url = 'https://assets.api.sennetconsortium.org/' + a['uuid'] + '/mudata_raw.h5mu'
    try:
        with fsspec.open(url, 'rb', block_size=2**20).open() as handle:
            with h5py.File(handle, 'r') as h:
                record['modalities'] = {key: dict(shape=list(map(int, h['mod'][key]['X'].attrs['shape'])),
                    layers=list(h['mod'][key]['layers'])) for key in h['mod']}
                record['paired_barcode_order'] = all(
                    list(h['mod'][key]['obs']['_index'][:]) == list(h['mod']['rna']['obs']['_index'][:])
                    for key in h['mod'])
    except Exception as error:
        record['error'] = str(error)
    print(record, flush=True)
    return record

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    records = list(pool.map(inspect, datasets))
(root / 'sennet_input_audit.json').write_text(json.dumps(records, indent=2)+'\n')

sets, details = [], []
for accession, label in [('GSM6997754','Ctrl'),('GSM6997756','Ctrl3'),('GSM6997758','Fontan1'),
                         ('GSM6997760','Fontan2'),('GSM6997762','Fontan3'),('GSM6997764','Fontan4')]:
    name = f'{accession}_ATAC_{label}_features.txt.gz'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6997nnn/{accession}/suppl/{name}'
    payload = urllib.request.urlopen(url, timeout=60).read()
    (root / name).write_bytes(payload)
    peaks = set(gzip.decompress(payload).decode().splitlines())
    sets.append(peaks)
    details.append(dict(accession=accession, source=url, unique_peaks=len(peaks)))
result = dict(samples=details, exact_shared_peaks=len(set.intersection(*sets)),
              union_peaks=len(set.union(*sets)), identical_peak_universes=all(x==sets[0] for x in sets))
(root / 'fontan_peak_audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(result, flush=True)
