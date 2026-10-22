"""Apply primary-source diagnoses to metadata copies; never rewrite old releases."""
import csv
from pathlib import Path
import shutil
import hashlib
import re

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014'
OLD = ROOT / 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/metadata'
HEALTHY = {'H06', 'H07', 'H10', 'H22'}
STEATOSIS = {'H04', 'H30', 'H32'}
ORDINAL = {'Healthy': '0', 'Steatosis': '1', 'Steatohepatitis': '2', 'Cirrhosis': '3', 'Unresolved': 'NA'}


def read(path, delimiter='\t'):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main():
    evidence = read(OUT / 'source/Liver_Atlas_TableS6_donor_join.tsv')
    labels = {d['donor']: ('Healthy' if d['source_sample'] in HEALTHY else
              'Steatosis' if d['source_sample'] in STEATOSIS else
              'Cirrhosis' if d['source_sample'] == 'H21' else 'Unresolved') for d in evidence}
    titles = {d['characteristics: shortFileName']: d['title'] for d in read(ROOT / 'Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv')}
    run_donors = {d['Run']: 'Liver_Atlas_' + re.search(r'H\d+', titles[d['shortfilename']]).group() for d in read(ROOT / 'Liver_Atlas/metadata/SraRunTable.csv', ',') if re.search(r'H\d+', titles.get(d['shortfilename'], ''))}
    inputs = [OLD / 'donor_metadata_stage_corrected.tsv', OLD / 'donor_metadata_stage_corrected_donor.tsv', ROOT / 'Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/current_source_diagnosis_donor_roster.tsv']
    for path in inputs:
        original = OUT / 'originals' / path.name
        if not original.exists():
            shutil.copy2(path, original)
        rows = read(original)
        keys = list(rows[0])
        for key in ('disease_stage_source', 'source_diagnosis'):
            if key not in keys:
                keys.append(key)
        for row in rows:
            unit = row.get('sample', row.get('donor'))
            if row['dataset'] == 'GSE174748' and unit in {'SRR14602057', 'SRR14602058'}:
                label = 'Cirrhosis'
                row['source_diagnosis'] = 'NAFLD cirrhotic; Filliol Nature 2022 Methods snRNA-seq analysis'
                row['disease_stage_source'] = 'doi:10.1038/s41586-022-05289-6'
            elif row['dataset'] == 'Liver_Atlas':
                donor = run_donors.get(unit, unit)
                assert donor in labels, (unit, donor)
                label = labels[donor]
                row['exclude_stage_analysis'] = 'TRUE'
                row['source_diagnosis'] = 'Guilliams Table S6; descriptive only; see exact donor companion'
                row['disease_stage_source'] = 'doi:10.1016/j.cell.2021.12.018 Table S6'
            else:
                continue
            row['disease_stage_coarse'] = label
            for field in ('disease_stage_numeric', 'stage_ordinal'):
                if field in row:
                    row[field] = ORDINAL[label]
        with (OUT / 'metadata' / path.name).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=keys, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    with (OUT / 'source/SHA256SUMS').open('w') as handle:
        for path in sorted((OUT / 'source').iterdir()):
            if path.name != 'SHA256SUMS':
                handle.write(hashlib.sha256(path.read_bytes()).hexdigest() + '  ' + path.name + '\n')


if __name__ == '__main__':
    main()
