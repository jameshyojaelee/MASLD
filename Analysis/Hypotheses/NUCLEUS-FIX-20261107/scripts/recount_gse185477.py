"""Jointly count the runs of one source-verified physical capture."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
TASK = ROOT / 'Analysis/Hypotheses/NUCLEUS-FIX-20261107'
assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
donor = sys.argv[1]
pairing = list(csv.DictReader((ROOT / 'data/GSE185477/metadata/donor_pairing.csv').open()))
row, = [r for r in pairing if r['donor_id'] == donor]
runs = row['rna_srrs'].split(';')
ena = list(csv.DictReader((TASK / 'checks/GSE185477_ena.tsv').open(), delimiter='\t'))
evidence = [r for r in ena if r['run_accession'] in runs]
assert {r['run_accession'] for r in evidence} == set(runs)
assert len({r['experiment_accession'] for r in evidence}) == 1
assert len({r['sample_accession'] for r in evidence}) == 1
assert len({r['sample_title'] for r in evidence}) == 1
out = TASK / 'recovery' / donor
out.mkdir(parents=True, exist_ok=False)
fastqs = out / 'fastqs'; fastqs.mkdir()
for run in runs:
    for mate in (1, 2):
        source = ROOT / f'data/GSE185477/fastq/{run}/{run}_S1_L001_R{mate}_001.fastq.gz'
        assert source.is_file(), source
        (fastqs / source.name).symlink_to(source)
cmd = ['/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger', 'count',
       '--id=' + donor, '--transcriptome=/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A',
       '--fastqs=' + str(fastqs), '--sample=' + ','.join(runs),
       '--create-bam=false', '--nosecondary', '--chemistry=auto',
       '--localcores=8', '--localmem=56', '--output-dir=' + str(out / 'count')]
(out / 'source_and_command.json').write_text(json.dumps(dict(source=evidence, command=cmd), indent=2))
subprocess.run(cmd, check=True, cwd=out)
(out / 'completion.json').write_text(json.dumps(dict(library=donor, status='completed', n_runs=len(runs)), indent=2))
