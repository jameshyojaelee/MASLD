"""Independently check existing r2 outputs; never change their contents."""
import collections
import csv
import hashlib
import json
import pathlib
import sys

release = pathlib.Path('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z')
output = pathlib.Path(sys.argv[1])
output.mkdir(parents=True, exist_ok=False)


def sha256(path):
    sha = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


checks = []
for line in (release / 'output_sha256.tsv').read_text().splitlines():
    expected, path = line.split(maxsplit=1)
    checks.append(dict(path=path, expected=expected, observed=sha256(path)))
assert all(x['expected'] == x['observed'] for x in checks), 'r2 named outputs changed'
studies = collections.Counter()
sources = collections.Counter()
positive_genes = set()
seen = set()
duplicates = []
with (release / 'susie_coloc_all_gwas.csv').open() as handle:
    for row in csv.DictReader(handle):
        key = row['gwas_name'], row['chr'], row['ensembl']
        if key in seen:
            duplicates.append(key)
        seen.add(key)
        studies[row['gwas_name']] += 1
        sources[row['eligibility_source']] += 1
        try:
            if float(row['PP.H4.susie']) > 0.5:
                positive_genes.add(row['ensembl'])
        except ValueError:
            pass
assert not duplicates, 'Duplicate study/chromosome/gene rows'
assert len(studies) == 50, 'Incomplete study portfolio'
assert len(positive_genes) == 708
assert sources['rerun_shared_posterior_check'] == 1307
result = dict(studies=dict(studies), master_rows=len(seen),
              full_portfolio_primary_genes=len(positive_genes),
              eligibility_sources=dict(sources), duplicate_rows=len(duplicates),
              checked_output_hashes=checks)
(output / 'review.json').write_text(json.dumps(result, indent=2))
(output / 'COMPLETE.txt').write_text('All 50 studies present; unique study/chromosome/gene keys; 708 full-portfolio primary genes; 1,307 corrected replacements; recorded output hashes matched.\n')
print(json.dumps({k: v for k, v in result.items() if k not in ('studies', 'checked_output_hashes')}))
