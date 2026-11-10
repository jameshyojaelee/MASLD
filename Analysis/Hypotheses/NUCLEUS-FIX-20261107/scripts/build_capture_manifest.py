"""Resolve physical captures from deposited experiment identities."""
import csv
import json
from pathlib import Path

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
TASK = ROOT / 'Analysis/Hypotheses/NUCLEUS-FIX-20261107'
OLD = ROOT / 'Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata'
REC = ROOT / 'Analysis/ATAC/Regulatory_Methods/results/20261025/final_masld_r1'
read = lambda p: list(csv.DictReader(p.open(), delimiter='\t'))
libraries = read(OLD / 'library_metadata_source_F.tsv')
source = {}
for ds in ('GSE185477', 'GSE202379', 'GSE174748', 'GSE189600', 'Liver_Atlas', 'GSE244832'):
    path = REC / 'rna_source_audit/ena_run_report.tsv' if ds == 'GSE244832' else TASK / 'checks' / (ds + '_ena.tsv')
    for row in read(path):
        source[(ds, row['run_accession'])] = row['experiment_accession']
groups = {}
for row in libraries:
    ds, run = row['dataset'], row['sample']
    cap = run if ds == 'GSE136103' else source[(ds, run)]
    groups.setdefault((ds, cap), []).append(row)
completed = {}
for marker in (REC / 'rna_recovery').glob('*/completion.json'):
    entry = json.loads(marker.read_text())
    if entry['status'] == 'completed': completed[entry['library']] = marker.parent
pair244 = list(csv.DictReader((ROOT / 'data/GSE244832/metadata/donor_pairing.csv').open()))
pair185 = list(csv.DictReader((ROOT / 'data/GSE185477/metadata/donor_pairing.csv').open()))
rows = []
adt = {r['run'] for r in json.loads((TASK / 'checks/liver_assay_identity.json').read_text()) if r['assay'] == 'ADT'}
for (ds, cap), group in sorted(groups.items()):
    excluded = sorted(r['sample'] for r in group if r['sample'] in adt)
    runs = sorted(r['sample'] for r in group if r['sample'] not in adt)
    assert runs, (ds, cap, 'no RNA library')
    assert len({r['donor'] for r in group}) == 1, (ds, cap, 'conflicting donor keys')
    for field in ('condition_harmonized', 'source_diagnosis', 'preparation_method'):
        assert len({r[field] for r in group}) == 1, (ds, cap, field)
    if len(runs) > 1:
        if ds == 'GSE244832':
            match, = [r for r in pair244 if set(r['rna_srrs'].split(';')) == set(runs)]
            directory = completed[match['rna_sample']]
        elif ds == 'GSE185477':
            match, = [r for r in pair185 if set(r['rna_srrs'].split(';')) == set(runs)]
            directory = TASK / 'recovery' / match['donor_id']
        else:
            directory = TASK / 'recovery' / cap
        matrix = directory / 'count/outs/filtered_feature_bc_matrix.h5'
    else:
        base = ROOT / ('Liver_Atlas/cellranger/human' if ds == 'Liver_Atlas' else f'data/{ds}/cellranger')
        if ds == 'GSE189600': base /= 'human'
        matrix = base / runs[0] / 'outs/filtered_feature_bc_matrix.h5'
    rows.append(dict(dataset=ds, capture=cap, representative_run=runs[0], runs=';'.join(runs),
                     excluded_runs=';'.join(excluded), donor=group[0]['donor'], n_runs=len(runs), matrix=str(matrix)))
with (TASK / 'capture_manifest.tsv').open('x') as f:
    writer = csv.DictWriter(f, rows[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
print('captures:', len(rows), 'run records:', sum(r['n_runs'] for r in rows),
      'repeated captures:', sum(r['n_runs'] > 1 for r in rows))
