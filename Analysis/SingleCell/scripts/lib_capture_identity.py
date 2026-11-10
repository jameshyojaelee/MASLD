"""Select one molecule-count matrix per source-verified physical capture."""
import csv
from pathlib import Path


def select_captures(samples, manifest):
    with Path(manifest).open() as handle:
        captures = list(csv.DictReader(handle, delimiter='\t'))
    human = {s['srr']: s for s in samples if s['species'] == 'human'}
    result = [s for s in samples if s['species'] != 'human']
    covered = set()
    for row in captures:
        runs = row['runs'].split(';')
        excluded = set(filter(None, row.get('excluded_runs', '').split(';')))
        source_runs = set(runs) | excluded
        if covered.intersection(source_runs):
            raise ValueError('A run belongs to more than one physical capture')
        covered.update(source_runs)
        # Some cohorts have separate importers (e.g. public GSE136103 matrices).
        # Qualify captures present in this producer without inventing new inputs.
        if not set(runs).intersection(human):
            continue
        if not set(runs).issubset(human):
            raise ValueError(f"Capture is missing source runs: {row['capture']}")
        representative = row['representative_run']
        sample = human[representative].copy()
        if any(human[r]['dataset'] != row['dataset'] for r in runs):
            raise ValueError('Capture/source dataset mismatch')
        matrix = Path(row['matrix'])
        if not matrix.is_file():
            raise FileNotFoundError(f'Joint count matrix is not complete: {matrix}')
        sample.update(h5_path=str(matrix), capture=row['capture'], source_runs=row['runs'],
                      count_state='one_joint_matrix_per_physical_capture')
        result.append(sample)
    unknown = set(human) - covered
    if unknown:
        raise ValueError(f'Human runs need source capture qualification before integration: {sorted(unknown)}')
    return result
