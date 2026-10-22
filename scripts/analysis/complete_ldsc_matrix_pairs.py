"""Remove incompatible rsID allele joins, then fit native pairwise LDSC."""
import csv
import gzip
import hashlib
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path.cwd()
OUT = pathlib.Path(sys.argv[1])
PYTHON = '/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin/python'
LDSC = '/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin/ldsc.py'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def digest(path):
    sha = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def read(path):
    with gzip.open(path, 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        rows = {}
        duplicates = set()
        for row in reader:
            key = row['SNP']
            if key in rows:
                duplicates.add(key)
            rows[key] = row
        for key in duplicates:
            del rows[key]
    return fields, rows, len(duplicates)


def compatible(first, second):
    a = first['A1'] + first['A2']
    b = second['A1'] + second['A2']
    valid = lambda pair: len(pair) == 2 and set(pair) <= set('ACGT') and pair[0] != pair[1] and pair[0] != pair[1].translate(COMPLEMENT)
    return valid(a) and valid(b) and a in (b, b[::-1], b.translate(COMPLEMENT), b[::-1].translate(COMPLEMENT))


def write(path, fields, rows):
    with gzip.open(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


matrix = pathlib.Path('GWAS/ldsc/rg_matrix_pdff_corrected_20261014/rg_12trait_66pair_candidate.tsv')
pairs = [row for row in csv.DictReader(matrix.open(), delimiter='\t') if row['state'] == 'historical_missing']
assert len(pairs) == 27
inputs = sorted({ROOT / 'GWAS/ldsc/munged' / (row[key] + '.sumstats.gz') for row in pairs for key in ('trait1', 'trait2')})
before = {str(path): digest(path) for path in inputs}
(OUT / 'input_sha256_before.json').write_text(json.dumps(before, indent=2))
subprocess.run([PYTHON, '-m', 'pip', 'freeze'], stdout=(OUT / 'pip_freeze.txt').open('w'), check=True)
report = []
for row in pairs:
    first_trait, second_trait = row['trait1'], row['trait2']
    name = first_trait + '__' + second_trait
    first_source = ROOT / 'GWAS/ldsc/munged' / (first_trait + '.sumstats.gz')
    second_source = ROOT / 'GWAS/ldsc/munged' / (second_trait + '.sumstats.gz')
    first_fields, first_rows, first_duplicates = read(first_source)
    second_fields, second_rows, second_duplicates = read(second_source)
    shared = sorted(first_rows.keys() & second_rows.keys())
    keep = [key for key in shared if compatible(first_rows[key], second_rows[key])]
    excluded = [key for key in shared if not compatible(first_rows[key], second_rows[key])]
    with (OUT / (name + '_excluded.tsv')).open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('SNP', 'first_A1', 'first_A2', 'second_A1', 'second_A2'))
        for key in excluded:
            writer.writerow((key, first_rows[key]['A1'], first_rows[key]['A2'], second_rows[key]['A1'], second_rows[key]['A2']))
    first = OUT / (name + '_first.sumstats.gz')
    second = OUT / (name + '_second.sumstats.gz')
    write(first, first_fields, (first_rows[key] for key in keep))
    write(second, second_fields, (second_rows[key] for key in keep))
    command = [PYTHON, LDSC, '--rg', str(first) + ',' + str(second),
               '--ref-ld-chr', 'GWAS/ldsc_ref/eur_w_ld_chr/',
               '--w-ld-chr', 'GWAS/ldsc_ref/eur_w_ld_chr/', '--out', str(OUT / name)]
    (OUT / (name + '_command.json')).write_text(json.dumps(command, indent=2))
    with (OUT / (name + '.stdout')).open('w') as handle:
        status = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT).returncode
    report.append(dict(trait1=first_trait, trait2=second_trait, shared=len(shared), compatible=len(keep), excluded=len(excluded), first_duplicate_ids=first_duplicates, second_duplicate_ids=second_duplicates, returncode=status))
    (OUT / 'pair_filtering.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report[-1]), flush=True)
    del first_rows, second_rows
(OUT / 'input_sha256_after.json').write_text(json.dumps({str(path): digest(path) for path in inputs}, indent=2))
assert before == json.loads((OUT / 'input_sha256_after.json').read_text()), 'Inputs changed during execution'
(OUT / 'COMPLETE.txt').write_text('Twenty-seven pairwise LDSC commands executed; input hashes unchanged. Inspect each native heritability and terminal outcome.\n')
