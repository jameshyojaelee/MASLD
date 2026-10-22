"""Reconcile the fixed 66-pair LDSC family, preserving historical successful estimates."""
import csv
import hashlib
import json
import math
import pathlib
import re
import sys

run = pathlib.Path(sys.argv[1])
base = pathlib.Path('GWAS/ldsc/rg_matrix_pdff_corrected_20261014/rg_12trait_66pair_candidate.tsv')
rows = list(csv.DictReader(base.open(), delimiter='\t'))
filtering = json.loads((run / 'pair_filtering.json').read_text())
assert len(rows) == 66 and len(filtering) == 27
assert json.loads((run / 'input_sha256_before.json').read_text()) == json.loads((run / 'input_sha256_after.json').read_text())
new = {}
for item in filtering:
    key = tuple(sorted((item['trait1'], item['trait2'])))
    name = item['trait1'] + '__' + item['trait2']
    log = run / (name + '.log')
    text = log.read_text()
    h2 = re.findall(r'Total Observed scale h2: ([^ ]+) \(([^)]+)\)', text)
    exceptions = [line for line in text.splitlines() if any(token in line for token in ('KeyError:', 'Traceback', 'ValueError:', 'FloatingPointError:'))]
    estimates = [line.split() for line in text.splitlines() if len(line.split()) >= 6 and line.split()[0].endswith('.sumstats.gz') and line.split()[1].endswith('.sumstats.gz')]
    assert estimates, 'Missing native terminal result: ' + name
    result = estimates[-1]
    positive_h2 = len(h2) == 2 and all(float(value) > 0 for value, se in h2)
    finite = all(value not in ('NA', 'nan', 'NaN') and math.isfinite(float(value)) for value in (result[2], result[3], result[5]))
    if positive_h2 and finite:
        assert not exceptions and item['returncode'] == 0, (name, exceptions)
        state = 'estimable'
        reason = 'Native estimate; positive estimated heritability in both traits'
    elif len(h2) == 2 and not positive_h2:
        state = 'untestable_nonpositive_h2'
        reason = 'Native LDSC genetic correlation undefined with nonpositive estimated heritability'
    else:
        raise AssertionError('Unresolved computational failure: ' + name + ': ' + repr(exceptions))
    new[key] = dict(rg=result[2], se=result[3], nominal_p=result[5], state=state, source=str(log), disposition_reason=reason, first_h2=h2[0][0], first_h2_se=h2[0][1], second_h2=h2[1][0], second_h2_se=h2[1][1])
for row in rows:
    key = tuple(sorted((row['trait1'], row['trait2'])))
    if row['state'] == 'historical_missing':
        row.update(new[key])
    else:
        source = pathlib.Path(row['source'])
        if source.suffix == '.tsv':
            native = next(item for item in csv.DictReader(source.open(), delimiter='\t') if item['trait'] in key)
            first_h2, first_se, second_h2, second_se = native['pdff_h2'], native['pdff_h2_se'], native['other_h2'], native['other_h2_se']
        else:
            text = source.read_text()
            first_h2, first_se = re.findall(r'Total Observed scale h2: ([^ ]+) \(([^)]+)\)', text)[0]
            native = next(line.split() for line in text.splitlines() if len(line.split()) >= 8 and all(line.split()[i].endswith('.sumstats.gz') for i in (0, 1)) and tuple(sorted(pathlib.Path(line.split()[i]).name.replace('.sumstats.gz', '') for i in (0, 1))) == key)
            second_h2, second_se = native[6], native[7]
        assert float(first_h2) > 0 and float(second_h2) > 0, 'Existing estimate has nonpositive native h2: ' + repr(key)
        row.update(disposition_reason='Existing successful native estimate retained; positive estimated heritability in both traits', first_h2=first_h2, first_h2_se=first_se, second_h2=second_h2, second_h2_se=second_se)

for row in rows:
    row['estimate_qualification'] = 'Unconstrained rg estimate outside [-1,1]; retain SE and do not clip' if row['state'] == 'estimable' and abs(float(row['rg'])) > 1 else ''
fields = list(rows[0])
assert all(row['state'] != 'historical_missing' for row in rows)
with (run / 'rg_12trait_66pair_terminal.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
(run / 'matrix_source_sha256.json').write_text(json.dumps({str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted({pathlib.Path(row['source']) for row in rows})}, indent=2))
counts = {state: sum(row['state'] == state for row in rows) for state in sorted({row['state'] for row in rows})}
(run / 'CHECKED.txt').write_text('Fixed 12-trait,66-pair family: ' + json.dumps(counts) + '. All 11 corrected PDFF pairs retained; no nominal P adjustment invented.\n')
print(json.dumps(counts))
