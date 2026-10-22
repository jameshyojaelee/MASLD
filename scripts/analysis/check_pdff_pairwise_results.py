"""Verify eight native LDSC fits and the exact installed allele-eligibility rule."""
import ast
import csv
import hashlib
import itertools
import json
import pathlib
import re
import sys

out = pathlib.Path(sys.argv[1])
source = pathlib.Path('/gpfs/commons/home/jameslee/micromamba/envs/ldsc/lib/python2.7/site-packages/ldscore/sumstats.py')
native = source.read_text()
namespace = {'it': itertools}
exec(native[native.index('_N_CHR ='):native.index('def _splitp')], namespace)
producer = pathlib.Path('scripts/analysis/pdff_pairwise_compatible.py')
tree = ast.parse(producer.read_text())
function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'compatible')
ours = {'COMPLEMENT': str.maketrans('ACGT', 'TGCA')}
exec(compile(ast.Module(body=[function], type_ignores=[]), str(producer), 'exec'), ours)
checks = 0
for alleles in itertools.product('ACGT', repeat=4):
    first = dict(A1=alleles[0], A2=alleles[1])
    second = dict(A1=alleles[2], A2=alleles[3])
    assert ours['compatible'](first, second) == (''.join(alleles) in namespace['MATCH_ALLELES'])
    checks += 1
records = json.loads((out / 'pair_filtering.json').read_text())
assert len(records) == 8
rows = []
for record in records:
    log = (out / (record['trait'] + '.log')).read_text()
    assert 'ERROR computing rg' not in log and 'Traceback' not in log, record['trait']
    match = re.search(r'Genetic Correlation: ([^ ]+) \(([^)]+)\)', log)
    assert match, record['trait']
    h2 = re.findall(r'Total Observed scale h2: ([^ ]+) \(([^)]+)\)', log)
    assert len(h2) == 2
    rg, se = map(float, match.groups())
    assert se > 0 and all(float(value[0]) > 0 for value in h2)
    p = float(re.search(r'^P: ([^\n]+)', log, re.MULTILINE).group(1))
    rows.append(dict(trait=record['trait'], rg=rg, se=se, nominal_p=p,
                     pdff_h2=float(h2[0][0]), pdff_h2_se=float(h2[0][1]),
                     other_h2=float(h2[1][0]), other_h2_se=float(h2[1][1]),
                     compatible_ids=record['compatible'], excluded_ids=record['excluded']))
with (out / 'pairwise_estimates.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
assert json.loads((out / 'input_sha256_before.json').read_text()) == json.loads((out / 'input_sha256_after.json').read_text())
hashes = {}
for path in [source] + sorted(pathlib.Path('GWAS/ldsc_ref/eur_w_ld_chr').glob('*')):
    if path.is_file():
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        hashes[str(path)] = digest.hexdigest()
(out / 'ldsc_reference_sha256.json').write_text(json.dumps(hashes, indent=2))
(out / 'CHECKED.txt').write_text(f'{checks} exhaustive allele-rule combinations match installed LDSC; eight positive-heritability estimates; no LDSC errors; original sumstats hashes unchanged. Reference hashes recorded after fitting.\n')
print(json.dumps(rows, indent=2))
