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
TRAITS = ('mvp_nafld', 'mvp_cirr', 'finngen_cirr', 'ukbb_alt', 'mvp_alt',
          'ukbb_ast', 'mvp_ast', 'finngen_obesity')
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


inputs = [ROOT / 'GWAS/ldsc/munged' / (trait + '.sumstats.gz')
          for trait in ('pdff',) + TRAITS]
before = {str(path): digest(path) for path in inputs}
assert hashlib.md5(inputs[0].read_bytes()).hexdigest() == '0d2e4e00efedd6207ef759068d377b23', 'PDFF input is not the corrected file'
(OUT / 'input_sha256_before.json').write_text(json.dumps(before, indent=2))
subprocess.run([PYTHON, '-m', 'pip', 'freeze'], stdout=(OUT / 'pip_freeze.txt').open('w'), check=True)
pdff_fields, pdff, pdff_duplicates = read(inputs[0])
report = []
for trait, source in zip(TRAITS, inputs[1:]):
    fields, other, duplicates = read(source)
    shared = sorted(pdff.keys() & other.keys())
    keep = [key for key in shared if compatible(pdff[key], other[key])]
    excluded = [key for key in shared if not compatible(pdff[key], other[key])]
    with (OUT / (trait + '_excluded.tsv')).open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('SNP', 'pdff_A1', 'pdff_A2', 'other_A1', 'other_A2'))
        for key in excluded:
            writer.writerow((key, pdff[key]['A1'], pdff[key]['A2'], other[key]['A1'], other[key]['A2']))
    first = OUT / (trait + '_pdff.sumstats.gz')
    second = OUT / (trait + '.sumstats.gz')
    write(first, pdff_fields, (pdff[key] for key in keep))
    write(second, fields, (other[key] for key in keep))
    command = [PYTHON, LDSC, '--rg', str(first) + ',' + str(second),
               '--ref-ld-chr', 'GWAS/ldsc_ref/eur_w_ld_chr/',
               '--w-ld-chr', 'GWAS/ldsc_ref/eur_w_ld_chr/', '--out', str(OUT / trait)]
    (OUT / (trait + '_command.json')).write_text(json.dumps(command, indent=2))
    with (OUT / (trait + '.stdout')).open('w') as handle:
        status = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT).returncode
    report.append(dict(trait=trait, shared=len(shared), compatible=len(keep),
                       excluded=len(excluded), pdff_duplicate_ids=pdff_duplicates,
                       other_duplicate_ids=duplicates, returncode=status))
    print(json.dumps(report[-1]), flush=True)
    del other
(OUT / 'pair_filtering.json').write_text(json.dumps(report, indent=2))
after = {str(path): digest(path) for path in inputs}
(OUT / 'input_sha256_after.json').write_text(json.dumps(after, indent=2))
assert before == after, 'Inputs changed during execution'
assert all(row['returncode'] == 0 for row in report), 'At least one pair failed; inspect logs'
(OUT / 'COMPLETE.txt').write_text('Eight pairwise LDSC commands completed; input hashes unchanged. Inspect heritability and warnings before interpretation.\n')
